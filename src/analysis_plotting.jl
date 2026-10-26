#
# General head section of any file within this module 
#

#
# Using statements
#
using ProgressBars
using Plots

#
# Include statements
#
include("analysis_helper.jl")

#
# Export statements
#
export plot_RSA_run
export plot_single_molecule
export animate_RSA_run
export plot_count_area_histograms
export plot_count_area_convergence
export plot_single_run_convergence
export plot_effective_gap_size
export gif
export savefig
export write_RSA_structures

#
# Definition of global variables
#

# Names of the metrics of a data set which are accepted by the analysis functions
const statistics_property_names = ["means", "standarddeviations", "minvalues", "minvalue_ids", "maxvalues", "maxvalue_ids", "avgvalues", "avgvalue_ids"]
const statistics_property_labels = ["Mean value", "Standard deviation", "Minimum value", "Run with the minimum value", "Maximum value", "Run with the maximum value", "Most average value", "Most average run"]

#
# Specific section for this file
#

# A function to get the coordinates of a displaced moleculue based on a status entry
function displaced_molecule_status(status, status_element, molecules, grids, lattice)

    # Get the information
    molecule_id, grid_id, point_id, rotation_id = status[:, status_element]

    # Get the coordinates
    molec_coords = molecules[molecule_id].coordinates_rotated[rotation_id]

    # Get the elements
    molec_elements = molecules[molecule_id].elements_sorted

    # Get the gridpoint
    gridpoint_coords = @view grids[grid_id].points[:, point_id]

    # Move structure to gridpoint
    displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords, [0.0, 0.0, 0.0], lattice.dimension)

    # Return coordinates and elements
    return molec_coords, displaced_coordinates, molec_elements, gridpoint_coords, molecule_id, grid_id, point_id

end

# A function to plot the points of all grids
function plot_grid_points(Ngrids, grids, lattice; pixel_per_angstrom = 10.0, silent = true, withmargins = false)
    
    # Define the image resolution
    x_axis_size = lattice.transcellvectors[1,1] + lattice.transcellvectors[1,2]
    y_axis_size = lattice.transcellvectors[2,1] + lattice.transcellvectors[2,2]
    x_axis_resolution = x_axis_size * pixel_per_angstrom
    y_axis_resolution = y_axis_size * pixel_per_angstrom
    
    if silent != true
        println("Resolution: " * string(x_axis_resolution) * " x " * string(y_axis_resolution))
    end

    # Define color of grid points
    grid_palette = palette(:darktest, Ngrids)

    # Define the white space around the simulation cell
    plot_margin, plot_framestyle, plot_ticks = simulation_cell_attributes(withmargins)

    # Plot the grids
    final_plot = 0
    for grid_id in 1:Ngrids
        if grid_id == 1
            final_plot = scatter(grids[grid_id].points[1,:],grids[grid_id].points[2,:], markersize=pixel_per_angstrom/10, legend=false, showaxis=false, grid=false, size=(x_axis_resolution, y_axis_resolution), xlims=(0, x_axis_size), ylims=(0, y_axis_size), color = grid_palette[grid_id], widen = false, margin = plot_margin, framestyle = plot_framestyle, ticks = plot_ticks)
        else
            scatter!(grids[grid_id].points[1,:],grids[grid_id].points[2,:], markersize=pixel_per_angstrom/10, xlims=(0, x_axis_size), ylims=(0, y_axis_size), color = grid_palette[grid_id])
        end
    end

    # Return the grid plot
    return final_plot

end


# A function to plot all molecules on their selected gridpoints
"""

    plot_RSA_run(status, Ngrids, grids, Nmolecules, molecules, lattice)
    plot_RSA_run(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, boundary_cells = 1, silent=true, gridplot = nothing, withmargins = false)

Create an image of the surface covered by adsorbates.

# Input
- `status`: A status field of a rsa\\_run\\_results\\_struct object.
- `Ngrids`: Integer number of present grid types.
- `grids`: A grid_struct object.
- `Nmolecules`: Integer number of present molecule types.
- `molecules`: A molecule_struct object.
- `lattice`: A lattice_struct object.

# Optional input
- `pixel_per_angstrom`: Resolution of the image controlled by the number of pixels for a distance of 1 angstrom as floating point number.
- `boundary_cells`: Integer number of boundary cells used to plot periodic boundary conditions.
- `silent`: Bool flag to request additional output statements.
- `gridplot`: A prepared plot of the grid points to prevent the repeated generation of an identical background.
- `withmargins`: Flag to add a default white space around the simulation cell. By default no white space is added.

# Return values
- A plots object of the covered surface.
"""
function plot_RSA_run(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, boundary_cells = 1, silent=true, gridplot = nothing, withmargins = false)

    # Define the axes size
    x_axis_size = lattice.transcellvectors[1,1] + lattice.transcellvectors[1,2]
    y_axis_size = lattice.transcellvectors[2,1] + lattice.transcellvectors[2,2]
    
    # A prepared plot of the grid points can be handed over to skip the generation of an identical background
    if gridplot === nothing
        final_plot = plot_grid_points(Ngrids, grids, lattice; pixel_per_angstrom = pixel_per_angstrom, silent = silent, withmargins = withmargins)
    else
        final_plot = deepcopy(gridplot)
    end 
    
    # Define color of molecule
    molecule_palette = palette(:darktest, Nmolecules)

    # Define the marker size of every molecule type
    markersize_per_molecule = Vector{Vector{Float64}}(undef, Nmolecules)
    for molecule_id in 1:Nmolecules
        markersize_per_molecule[molecule_id] = Float64.(atomic_information[molecules[molecule_id].elements_sorted[:], 3]) * pixel_per_angstrom
    end

    # Create the vectors collecting the atoms of all adsorbates of every molecule type
    x_coordinates = Vector{Vector{Float64}}(undef, Nmolecules)
    y_coordinates = Vector{Vector{Float64}}(undef, Nmolecules)
    markersizes = Vector{Vector{Float64}}(undef, Nmolecules)
    for molecule_id in 1:Nmolecules
        x_coordinates[molecule_id] = Vector{Float64}(undef, 0)
        y_coordinates[molecule_id] = Vector{Float64}(undef, 0)
        markersizes[molecule_id] = Vector{Float64}(undef, 0)
    end

    # Collect every molecule
    for element_id in axes(status, 2)

        # Get the displaced coordinates and elements
        molec_coords, displaced_coordinates, molec_elements, gridpoint_coords, molecule_id, grid_id, point_id = displaced_molecule_status(status, element_id, molecules, grids, lattice)
 
        # Add the adsorbate
        append!(x_coordinates[molecule_id], displaced_coordinates[1,:])
        append!(y_coordinates[molecule_id], displaced_coordinates[2,:])
        append!(markersizes[molecule_id], markersize_per_molecule[molecule_id])

        # Check whether the adsorbate was within the boundary cells and add all needed periodic images
        transx, transy = grids[grid_id].mapping[2:3,point_id]
        translations = boundary_translation_vectors(transx, transy, boundary_cells, lattice)
        for translation_id in eachindex(translations)
            gridpoint_coords_moved = gridpoint_coords + translations[translation_id]
            displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
            append!(x_coordinates[molecule_id], displaced_coordinates[1,:])
            append!(y_coordinates[molecule_id], displaced_coordinates[2,:])
            append!(markersizes[molecule_id], markersize_per_molecule[molecule_id])
         end

    end

    # Plot all adsorbates of every molecule type
    for molecule_id in 1:Nmolecules
        if isempty(x_coordinates[molecule_id])
            continue
        end
        scatter!(final_plot, x_coordinates[molecule_id], y_coordinates[molecule_id], markersize = markersizes[molecule_id], markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
    end

    # Return final plot
    return final_plot

end

# A function to create an animation of an rsa run
"""

    animate_RSA_run(stepinfo, Ngrids, grids, Nmolecules, molecules, lattice)
    animate_RSA_run(stepinfo, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, boundary_cells = 1, startstep = 1, laststep = 0, withmargins = false)

Create an animation of a RSA simulation.

# Input
- `stepinfo`: A stepinfo field of a rsa_run_results_struct.
- `Ngrids`: Number of present grid types.
- `grids`: A grid_struct object.
- `Nmolecules`: Number of present molecule types.
- `molecules`: A molecule_struct object.
- `lattice`: A lattice_struct object.

# Optional input
- `pixel_per_angstrom`: Resolution of the image controlled by the number of pixels for a distance of 1 angstrom.
- `boundary_cells`: Number of boundary cells used to plot periodic boundary conditions.
- `startstep`: First RSA step shown by the animation.
- `laststep`: Last RSA step shown by the animation. A value of zero (default) requests all steps of the given stepinfo.
- `withmargins`: Flag to add a default white space around the simulation cell. By default no white space is added.

# Return values
- A plots object containing the animation of the RSA simulation.

# Hints
- Generation of large animations is extremely slow.
- Only reasonable to use for the a few thousand steps.
"""
function animate_RSA_run(stepinfo, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, boundary_cells = 1, startstep = 1, laststep = 0, withmargins = false)

    # Throw a warning in case the resolution is getting to large
    if pixel_per_angstrom > 10
        println("Hint to the user: Your pixel_per_angstrom is larger then 10. Expect the generation of the animation to slow down.")
    end
    if pixel_per_angstrom > 100
        println("Hint to the user: Your pixel_per_angstrom has an extremely large value! The generation of the animation will be slow and you will get a large animation file. Consider reducing this value.")
    end

    # Get the number of performed RSA steps
    Nsteps = size(stepinfo,2)

    # Define the range of steps to be animated
    # A laststep of zero requests all steps of the given stepinfo
    if laststep ≤ 0 || laststep > Nsteps
        laststep = Nsteps
    end
    if startstep < 1
        startstep = 1
    end
    if startstep > laststep
        println("The first step of the animation (" * string(startstep) * ") is larger than the last step (" * string(laststep) * ").")
        error("Animation Range Error")
    end

    # Create the animation object
    anim = Animation()

    # Create the plot of the grid points once and use it as the background of every frame
    gridplot = plot_grid_points(Ngrids, grids, lattice; pixel_per_angstrom = pixel_per_angstrom, withmargins = withmargins)

    # Preallocate matrices
    realsize = 0
    status = Matrix{Int64}(undef, 4, Nsteps)

    # Update the status matrix for all steps in front of the first frame
    for step_id in 1:startstep-1
        realsize = update_status_by_stepinfo!(status, realsize, stepinfo, step_id)
    end

    # Create the frames
    for frame_id in ProgressBar(startstep:laststep)

        # Get the information for this step and update the status matrix
        realsize = update_status_by_stepinfo!(status, realsize, stepinfo, frame_id)

        # Create the frame
        substatus = @view status[1:4,1:realsize]
        newframe = plot_RSA_run(substatus, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = pixel_per_angstrom, boundary_cells = 1, silent=true, gridplot = gridplot, withmargins = withmargins)

        # Add the frame to the animation
        frame(anim, newframe)

    end

    # Return the animation
    return anim
    
end

# A function to plot a single molecule
"""

    plot_single_molecule(molecule_id, Nmolecules, molecules)
    plot_single_molecule(molecule_id, Nmolecules, molecules; pixel_per_angstrom = 10.0)

Create an image of a single adsorbate.

# Input
- `molecule_id`: Number of the molecule in the molecules vector to be used.
- `Nmolecules`: Number of present molecule types.
- `molecules`: A molecule_struct object.

# Optional input
- `pixel_per_angstrom`: Resolution of the image controlled by the number of pixels for a distance of 1 angstrom.

# Return values
- A plots object of the adsorbate.
"""
function plot_single_molecule(molecule_id, Nmolecules, molecules; pixel_per_angstrom = 10.0)

    # To be comparable to other plots the factor of pixel_per_angstrom is increased here
    pixel_per_angstrom *= 10.0

    # Get the coordinates
    molec_coords = molecules[molecule_id].coordinates

    # Get the elements
    molec_elements = molecules[molecule_id].elements

    # Calculate the centroid
    centroid = calculate_centroid(molec_coords)

    # Move the centroid to the origin of the coordinate system
    origin = [0.0, 0.0, 0.0]
    moved_molecule = move_structure_to_point(molec_coords, origin, centroid, molecules[molecule_id].dimension)
    #println(moved_molecule)

    # Define the image resolution
    markersize_vector = atomic_information[molec_elements[:], 3] * pixel_per_angstrom 

    # Get the largest distance of the molecule to the origin
    maxdistance, distances = get_largest_vdW_distance_to_point(molecules[molecule_id].Natoms, molec_elements, moved_molecule, origin, "2D")
    #println("Maxdistance: " * string(maxdistance))

    # Define the size of the axes
    axis_size = 2.0* maxdistance
    #println("Axis size: " * string(axis_size))

    # Define color of molecule
    molecule_palette = palette(:darktest, Nmolecules)  

    # Add to plot
    final_plot = scatter(moved_molecule[1,:],moved_molecule[2,:], markersize=markersize_vector, size=(axis_size * pixel_per_angstrom,axis_size * pixel_per_angstrom), xlims=(-0.5*axis_size,0.5*axis_size), ylims=(-0.5*axis_size,0.5*axis_size), markerstrokewidth=0, legend=false, showaxis=false,grid=false, color = molecule_palette[molecule_id])

    # Return final plot
    return final_plot

end

# A function to write all structures of a single run with their atomic coordinates to a file of xyz format
function write_RSA_structures(run_id, rsa_results, Nmolecules, molecules, file_path)

    # Get the number of adsorbates per molecule
    adsorbate_count = zeros(Int64, Nmolecules)
    for molecule_id in 1:Nmolecules
        adsorbate_count[molecule_id] = size(findall(x -> x == molecule_id, @view rsa_results[run_id].status[1,:]),1)
    end
    
    # Get the total number of atoms
    Total_NAtoms = 0
    for molecule_id in 1:Nmolecules
        Total_NAtoms += adsorbate_count[molecule_id] * molecules[molecule_id].Natoms
    end

    # Get the number of adsorbates to print
    Nadsorbates = size(rsa_results[run_id].status,2)

    # Open file
    io = open(file_path, "w")

    # Write header
    println(io, string(Total_NAtoms))
    println(io,"Atomic coordinates of RSA simulation: " * string(run_id))

    # Write every molecule
    for adsorbate_id in 1:Nadsorbates

        # Get the information
        molecule_id, grid_id, point_id, rotation_id = rsa_results[run_id].status[:,adsorbate_id]

        # Get the coordinates
        molec_coords = molecules[molecule_id].coordinates_rotated[rotation_id]

        # Get the elements
        molec_elements = molecules[molecule_id].elements_sorted

        # Get the gridpoint
        gridpoint_coords = @view grids[grid_id].points[:, point_id]

        # Move structure to gridpoint
        displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords, [0.0, 0.0, 0.0], molecules[molecule_id].dimension)
 
        # Write every atom coordinates
        for atom_id in axes(displaced_coordinates, 2)
            println(io, string(atomic_information[molec_elements[atom_id],2]) * "   " * join(string.(displaced_coordinates[:, atom_id]),"   "))
        end

    end

    # Close the file
    close(io)

end

# A function to plot a histogram based on a given vector of numbers
function plot_histogram(data_vector; labelx = "Value", labely = "Count", stepsize = 1.0, resolution = 600, distribution = "none", plotonly = true, threshold = Vector{Float64}(undef, 0), normalized = false)

    # Derive all metrics of the data set
    if distribution == "truncated"
        mean, standarddeviation, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id = calculate_data_set_metrics(data_vector; zerogaussian = true)
    elseif distribution == "gaussian"
        mean, standarddeviation, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id = calculate_data_set_metrics(data_vector; zerogaussian = false)
    else
        mean, standarddeviation = 0.0, 0.0
        minvalue, minvalue_id = findmin(data_vector)
        maxvalue, maxvalue_id = findmax(data_vector)
        avgvalue, avgvalue_id = 0.0, 0
    end

    # Define the bins of the histrogram and the limits of the x axis
    if abs(maxvalue - minvalue) < stepsize
        bins = range(minvalue - stepsize / 2, maxvalue + stepsize / 2, length = 2)
        lowerlimit = minvalue - 3 * stepsize / 2
        upperlimit = maxvalue + 3 * stepsize / 2
    else
        bins = range(minvalue, maxvalue, step = stepsize)
        lowerlimit = minvalue
        upperlimit = maxvalue
    end

    # Scale the gaussian distribution (to the area of the visible bins)
    if normalized == true
        plots_normalization = :pdf
        gaussian_scaling = 1.0
    else
        plots_normalization = :none
        gaussian_scaling = gaussian_scaling_factor(data_vector, bins)
    end

    # Plot the histogram
    histo = histogram(data_vector, xlabel = labelx, ylabel = labely, normalize = plots_normalization, xlims=(lowerlimit, upperlimit), legend=false, bins = bins, dpi = resolution)

    # Add the normal distribution
    if distribution == "gaussian" && standarddeviation > 1.0e-10
        x = range(lowerlimit, upperlimit, step = stepsize/10)
        y = @. gaussian_scaling * 1/(standarddeviation * sqrt(2*π)) * exp(-0.5 * (x - mean)^2 / standarddeviation^2) 
        plot!(x,y, width = 4, lc = "red")
    elseif distribution == "truncated" && standarddeviation > 1.0e-10
        x = range(lowerlimit, upperlimit, step = stepsize/10)
        y = @. gaussian_scaling * 2 * 1/(standarddeviation * sqrt(2*π)) * exp(-0.5 * (x - mean)^2 / standarddeviation^2) 
        plot!(x,y, width = 4, lc = "red")
    end

    # Add a vertical line if requested
    for ele in threshold
        vline!([ele], linestyle = :dash, linecolor = :red, linewidth = resolution/200)
    end

    # Return all results
    if plotonly == true
        return histo
    else
        return histo, mean, standarddeviation, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id
    end
    
end

# A function to create histogram plots for adsorbate count and surface coverage
"""

    plot_count_area_histograms(Nruns, rsa_results, Nmolecules, molecules, lattice)
    plot_count_area_histograms(Nruns, rsa_results, Nmolecules, molecules, lattice; status = true, plotonly = true, count = 1.0, area = 1.0, normalized = false)

Create histograms counting the number of adsorbed molecules and the covered area. For every molecule type each histogram is generated. In addition, a final set of histograms is generated for all molecule types combined.

# Input
- `Nruns`: Total number of RSA simulations.
- `rsa_results`: A rsa_run_results_struct object.
- `Nmolecules`: Number of present molecule types.
- `molecules`: A molecule_struct object.
- `lattice`: A lattice_struct object.

# Optional input
- `status`: Flag forcing the recalculation of the status based on the stepinfo field.
- `plotonly`: Flag to request additional metrics.
- `count`: Bin size for histogram showing molecule counts.
- `area`: Bin size for histograms showing covered area.
- `normalized`: Flag to normalize the histogram with Plots internals.

# Return values
- `plotonly = true (default)`: A vector of histograms is returned. Count and covered surface is contained pairwise for every molecule type while the second last element contains the total adsorbate count and the last element the total covered area.
- `plotonly = false`: In addition to the histogram vector, vectors storing the mean values, standarddeviation, min and max values, as well as the simulation closest to the mean value are returned. In the following order: histograms, means, standarddeviations, minvalues, minvalue_ids, maxvalues, maxvalue_ids, averagevalues, averagevalue_ids.
"""
function plot_count_area_histograms(Nruns, rsa_results, Nmolecules, molecules, lattice; status = true, plotonly = true, count = 1.0, area = 1.0, normalized = false)

    # Get the adsorbate count and the covered area (in %) of every run
    adsorbate_count_per_run, adsorbate_area_per_run = calculate_count_area_per_run(Nruns, rsa_results, Nmolecules, molecules, lattice; status = status) 

    # Generate total adsorbate count per run
    total_adsorbate_count = sum(adsorbate_count_per_run, dims = 1)

    # Generate total area per run
    total_area = sum(adsorbate_area_per_run, dims = 1)

    # Generate the vectors storing all information and plots
    if plotonly == true
        histos = Vector{Any}(undef, (Nmolecules * 2) + 2)
    else
        histos = Vector{Any}(undef, (Nmolecules * 2) + 2)
        means = Vector{Float64}(undef, (Nmolecules * 2) + 2)
        standarddeviations = Vector{Float64}(undef, (Nmolecules * 2) + 2)
        minvalues = Vector{Float64}(undef, (Nmolecules * 2) + 2)
        minvalue_ids = Vector{Int64}(undef, (Nmolecules * 2) + 2)
        maxvalues = Vector{Float64}(undef, (Nmolecules * 2) + 2)
        maxvalue_ids = Vector{Int64}(undef, (Nmolecules * 2) + 2)
        avgvalues = Vector{Float64}(undef, (Nmolecules * 2) + 2)
        avgvalue_ids = Vector{Int64}(undef, (Nmolecules * 2) + 2)
    end

    # Plot the histograms
    if plotonly == true
        
        plot_id = 0
        for molecule_id in 1:Nmolecules
            plot_id += 1
            histos[plot_id] = plot_histogram(adsorbate_count_per_run[molecule_id,:]; labelx = "Adsorbate count - molecule " * string(molecule_id), labely = histogram_label(normalized), stepsize = count, resolution = 600, plotonly = true, distribution = "gaussian", normalized = normalized)
            plot_id += 1
            histos[plot_id] = plot_histogram(adsorbate_area_per_run[molecule_id,:]; labelx = "Covered area in % - molecule " * string(molecule_id), labely = histogram_label(normalized), stepsize = area, resolution = 600, plotonly = true, distribution = "gaussian", normalized = normalized)
        end
        
        plot_id += 1
        histos[plot_id] = plot_histogram(total_adsorbate_count[:]; labelx = "Adsorbate count", labely = histogram_label(normalized), stepsize = count, resolution = 600, plotonly = true, distribution = "gaussian", normalized = normalized)
        plot_id += 1
        histos[plot_id] = plot_histogram(total_area[:]; labelx = "Covered area in %", labely = histogram_label(normalized), stepsize = area, resolution = 600, plotonly = true, distribution = "gaussian", normalized = normalized)

    else

        plot_id = 0
        for molecule_id in 1:Nmolecules
            plot_id += 1
            histos[plot_id], means[plot_id], standarddeviations[plot_id], minvalues[plot_id], minvalue_ids[plot_id], maxvalues[plot_id], maxvalue_ids[plot_id], avgvalues[plot_id], avgvalue_ids[plot_id] = plot_histogram(adsorbate_count_per_run[molecule_id,:]; labelx = "Adsorbate count - molecule " * string(molecule_id), labely = histogram_label(normalized), stepsize = count, resolution = 600, plotonly = false, distribution = "gaussian", normalized = normalized)
            plot_id += 1
            histos[plot_id], means[plot_id], standarddeviations[plot_id], minvalues[plot_id], minvalue_ids[plot_id], maxvalues[plot_id], maxvalue_ids[plot_id], avgvalues[plot_id], avgvalue_ids[plot_id] = plot_histogram(adsorbate_area_per_run[molecule_id,:]; labelx = "Covered area in % - molecule " * string(molecule_id), labely = histogram_label(normalized), stepsize = area, resolution = 600, plotonly = false, distribution = "gaussian", normalized = normalized)
        end
        
        plot_id += 1
        histos[plot_id], means[plot_id], standarddeviations[plot_id], minvalues[plot_id], minvalue_ids[plot_id], maxvalues[plot_id], maxvalue_ids[plot_id], avgvalues[plot_id], avgvalue_ids[plot_id] = plot_histogram(total_adsorbate_count[:]; labelx = "Adsorbate count", labely = histogram_label(normalized), stepsize = count, resolution = 600, plotonly = false, distribution = "gaussian", normalized = normalized)
        plot_id += 1
        histos[plot_id], means[plot_id], standarddeviations[plot_id], minvalues[plot_id], minvalue_ids[plot_id], maxvalues[plot_id], maxvalue_ids[plot_id], avgvalues[plot_id], avgvalue_ids[plot_id] = plot_histogram(total_area[:]; labelx = "Covered area in %", labely = histogram_label(normalized), stepsize = area, resolution = 600, plotonly = false, distribution = "gaussian", normalized = normalized)

    end

    # Return results
    if plotonly == true
        return histos
    else
        return histos, means, standarddeviations, minvalues, minvalue_ids, maxvalues, maxvalue_ids, avgvalues, avgvalue_ids
    end

end

# A function to plot the effective gap size
"""

    plot_effective_gap_size(status, Ngrids, grids, Nmolecules, molecules, lattice)
    plot_effective_gap_size(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, gapsonly = false, withstroke = true, plotonly = true, normalized = false, withmargins = false)

Create an image of the effective gap sizes as well as the histogram showing the frequency of all gap sizes.

# Input
- `status`: A status field of a rsa_run_results_struct.
- `Ngrids`: Number of present grid types.
- `grids`: A grid_struct object.
- `Nmolecules`: Number of present molecule types.
- `molecules`: A molecule_struct object.
- `lattice`: A lattice_struct object.

# Optional input
- `pixel_per_angstrom`: Resolution of the image controlled by the number of pixels for a distance of 1 angstrom.
- `gapsonly`: Flag to request a visualization of only the gaps (removing all adsorbates).
- `withstroke`: Flag to add a stroke to the visualization of the gap sizes.
- `plotonly`: Flag to request additional metrics.
- `normalized`: Flag to normalize the histogram with Plots internals.
- `withmargins`: Flag to add a default white space around the simulation cell. By default no white space is added.
- `threshold`: Flag to add threshold values to be marked in the histogram. Values specified in a vector.

# Return values
- `plotonly = true (default)`: Returns the histogram showing the frequency of all gap sizes and a plots object for the visualization of the gaps in the following order: histogram, plot.
- `plotonly = false`: In addition to the default case, a vector containing the obtained effective gap sizes as well as a vector of the corresponding free grid point are returned. Information are returned in the following order: histogram, plot, gap sizes, free grid points.
"""
function plot_effective_gap_size(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, gapsonly = false, withstroke = true, plotonly = true, stepsize = 0.2, threshold = Vector{Float64}(undef, 0), normalized = false, withmargins = false)

    # Get the shells of every grid point
    neighbour_shell_list = create_neighbour_shell_lists(Ngrids, grids, Nmolecules, molecules, lattice)

    # Get the effective gap sizes
    effective_gap_sizes, free_grid_points = calculate_effective_gap_size(status, neighbour_shell_list, Ngrids, grids, molecules, lattice; dim = 2)
    
    # Create the histogram of the effective gap sizes
    gaps_histogram = plot_histogram(effective_gap_sizes; labelx = "Effective Gap Size in Å", labely = histogram_label(normalized), stepsize = stepsize, threshold = threshold, normalized = normalized)

    # Get a plot of all adsorbates
    if gapsonly == false
        plot_gaps = plot_RSA_run(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = pixel_per_angstrom, boundary_cells = 1, silent=true, withmargins = withmargins)
    end

    # Collect coordinates of free grid points
    Nradii = size(effective_gap_sizes,1)
    coordinates = Matrix{Float64}(undef,2,Nradii)
    for element_id in axes(free_grid_points, 2)

        # Get the free grid point
        grid_id, point_id = free_grid_points[:,element_id]

        # Get the coordinates of the grid point
        coordinates[:,element_id] = grids[grid_id].points[1:2,point_id]       

    end

    # Add to the plot
    if gapsonly == false
        if withstroke == true
            scatter!(coordinates[1,:],coordinates[2,:], markersize = effective_gap_sizes * pixel_per_angstrom, markerstrokewidth = 10 / pixel_per_angstrom, color = :grey)
        else
            scatter!(coordinates[1,:],coordinates[2,:], markersize = effective_gap_sizes * pixel_per_angstrom, markerstrokewidth = 0, color = :grey)
        end
    else
        # Define the image resolution
        x_axis_size = lattice.transcellvectors[1,1] + lattice.transcellvectors[1,2]
        y_axis_size = lattice.transcellvectors[2,1] + lattice.transcellvectors[2,2]
        x_axis_resolution = x_axis_size * pixel_per_angstrom
        y_axis_resolution = y_axis_size * pixel_per_angstrom
    
        # Define color of grid points
        grid_palette = palette(:darktest, Ngrids)

        # Define the white space around the simulation cell
        plot_margin, plot_framestyle, plot_ticks = simulation_cell_attributes(withmargins)

        # Plot the grids
        #for grid_id in 1:Ngrids
        #    if grid_id == 1
        #        plot_gaps = scatter(grids[grid_id].points[1,:],grids[grid_id].points[2,:], markersize=pixel_per_angstrom/10, legend=false, showaxis=false, grid=false, size=(x_axis_resolution, y_axis_resolution), xlims=(0, x_axis_size), ylims=(0, y_axis_size), color = grid_palette[grid_id], widen = false)
        #    else
        #        scatter!(grids[grid_id].points[1,:],grids[grid_id].points[2,:], markersize=pixel_per_angstrom/10, color = grid_palette[grid_id])
        #    end
        #end

        # Add the gaps
        if withstroke == true
            plot_gaps = scatter(coordinates[1,:],coordinates[2,:], markersize = effective_gap_sizes * pixel_per_angstrom, markerstrokewidth = 10 / pixel_per_angstrom, legend=false, showaxis=false, grid=false, size=(x_axis_resolution, y_axis_resolution), xlims=(0, x_axis_size), ylims=(0, y_axis_size), widen = false, color = :grey, margin = plot_margin, framestyle = plot_framestyle, ticks = plot_ticks)
        else
            plot_gaps = scatter(coordinates[1,:],coordinates[2,:], markersize = effective_gap_sizes * pixel_per_angstrom, markerstrokewidth = 0, legend=false, showaxis=false, grid=false, size=(x_axis_resolution, y_axis_resolution), xlims=(0, x_axis_size), ylims=(0, y_axis_size), widen = false, color = :grey, margin = plot_margin, framestyle = plot_framestyle, ticks = plot_ticks)
        end

    end

    # Return results
    if plotonly == true
        return gaps_histogram, plot_gaps
    else
        return gaps_histogram, plot_gaps, effective_gap_sizes, free_grid_points
    end

end

function get_all_free_grid_points(status, molecules, Ngrids, grids, lattice, neighbour_shell_list)

    # Create the default vectors
    free_point = Vector{Vector{Bool}}(undef, Ngrids)
    for grid_id in 1:Ngrids
        free_point[grid_id] = [true for i in 1:grids[grid_id].Npoints]
    end

    # Set all points covered by an adsorbate to false
    for element_id in axes(status, 2)
        
        # Get all information
        molec_coords, displaced_coordinates, molec_elements, gridpoint_coords, molecule_id, grid_A_id, point_A_id = displaced_molecule_status(status, element_id, molecules, grids, lattice)
        self_unique_point, self_transx, self_transy = grids[grid_A_id].mapping[1:3, point_A_id]

        # Set the point itself to false
        free_point[grid_A_id][point_A_id] = false

        # Loop over grids and points of the first shell
        for grid_B_id in 1:Ngrids
            points_B = neighbour_shell_list[grid_A_id][self_unique_point][grid_B_id][1]
            for point_B_id in eachindex(points_B)
                
                # Map the shell point to the actual point
                # Here: Mapping has to use "+" as I shift to the final position of the adsorbate
                mapping_unique_point, mapping_transx, mapping_transy = grids[grid_B_id].mapping[1:3, points_B[point_B_id]]
                new_transx = mapping_transx + self_transx
                new_transy = mapping_transy + self_transy
                absolute_point_id = map_translation_to_gridpoint(mapping_unique_point, grids[grid_B_id].Nuniquepoints, new_transx, lattice.Ncellx, new_transy, lattice.Ncelly)

                # Get the coordinates of this point
                point_coords = @view grids[grid_B_id].points[:,absolute_point_id]

                # Check for overlap of the adsorbate and the grid point
                covered = point_covered_by_vdW_radii_2D(molec_elements, displaced_coordinates, point_coords, lattice)
                if covered == true
                    free_point[grid_B_id][absolute_point_id] = false
                end

            end
        end

    end

    # Return results
    return free_point

end

function find_first_shell_with_adsorbate(status, Ngrids, grids, lattice, neighbour_shell_list, grid_id, point_id)

    # Map the current point to its unique point
    self_unique_point, self_transx, self_transy = grids[grid_id].mapping[1:3, point_id]

    # Get the number of shells per grid
    Nshells = size.(neighbour_shell_list[grid_id][self_unique_point], 1)
    Nmax_shells = findmax(Nshells)[1] 

    # Find the closest shell with an adsorbate
    shell_occupied = 0
    for shell_id in 1:Nmax_shells
        
        # Loop over the grids
        for grid_B_id in 1:Ngrids

            # Skip this grid in case the maximum number of shells was reached
            if shell_id > Nshells[grid_B_id]
                continue
            end

            # Get all points of this shell
            shell_points = neighbour_shell_list[grid_id][self_unique_point][grid_B_id][shell_id]

            # Check whether any point is occupied
            for shell_point_id in eachindex(shell_points)

                # Get the mapping of this point
                mapping_unique_point, mapping_transx, mapping_transy = grids[grid_B_id].mapping[1:3, shell_points[shell_point_id]]

                # Get the actual grid point
                # Here: The mapping has to use a "+" as I shift to the final adsorbate position
                new_transx = mapping_transx + self_transx
                new_transy = mapping_transy + self_transy
                absolute_point_id = map_translation_to_gridpoint(mapping_unique_point, grids[grid_B_id].Nuniquepoints, new_transx, lattice.Ncellx, new_transy, lattice.Ncelly)

                # Is this point occupied by an adsorbate
                occupied = present_column(status, [grid_B_id, absolute_point_id], 2:3)
                if occupied == true
                    #println("Shell occupied")
                    #println(grid_B_id)
                    #println(shell_id)
                    #println(shell_point_id)
                    #println(absolute_point_id)
                    shell_occupied = shell_id
                    break
                end

            end

            # Break this loop in case an occupied shell was found
            if shell_occupied != 0
                break
            end

        end

        # Break this loop in case an occupied shell was found
        if shell_occupied != 0
            break
        end

    end

    # Return results
    return shell_occupied

end

function get_distance_closest_adsorbate(status, Ngrids, grids, lattice, molecules, neighbour_shell_list, grid_id, point_id, shell_occupied; dim = 2)

    # Debug
    #println("Searching distance to adsorbate for grid/point:")
    #println(grid_id)
    #println(point_id)
    
    # Map the current point to its unique point
    self_unique_point, self_transx, self_transy = grids[grid_id].mapping[1:3, point_id]
    
    # Get the coordinates of this point
    free_point_coords = grids[grid_id].points[:,point_id]

    # Get the number of shells per grid
    Nshells = size.(neighbour_shell_list[grid_id][self_unique_point], 1)

    # Get the distance to the closest adsorbate 
    mindistance = Inf
    for shell_id in shell_occupied:shell_occupied + 1

        #println("Shell ID: " * string(shell_id))

        # Loop over the grids
        for grid_B_id in 1:Ngrids

            #println("Grid B ID: " * string(grid_B_id))

            # Skip this grid in case the maximum number of shells was reached
            if shell_id > Nshells[grid_B_id]
                continue
            end

            # Get all points of this shell
            shell_points = neighbour_shell_list[grid_id][self_unique_point][grid_B_id][shell_id]

            # Find the occupied points
            for shell_point_id in eachindex(shell_points)

                #println("Shell point ID: " * string(shell_point_id))

                # Get the mapping of this point
                mapping_unique_point, mapping_transx, mapping_transy = grids[grid_B_id].mapping[1:3, shell_points[shell_point_id]]

                # Get the actual grid point
                # Here: The mapping has to use a "+" as I shift to the final adsorbate structure
                new_transx = mapping_transx + self_transx
                new_transy = mapping_transy + self_transy
                absolute_point_id = map_translation_to_gridpoint(mapping_unique_point, grids[grid_B_id].Nuniquepoints, new_transx, lattice.Ncellx, new_transy, lattice.Ncelly)

                #println("Absolute point ID: " * string(absolute_point_id))

                # Is this point occupied by an adsorbate
                status_id = findfirst_subset(status, [grid_B_id, absolute_point_id], 2:3)
                if status_id != 0
                    
                    # Generate the structure
                    molec_coords, displaced_coordinates, molec_elements, gridpoint_coords, molecule_id, grid_A_id, point_A_id = displaced_molecule_status(status, status_id, molecules, grids, lattice)
                    #println("Status ID: " * string(status_id))
                    
                    # Get the distance to the grid point
                    distance = get_smallest_vdW_distance_to_point(molec_elements, displaced_coordinates, free_point_coords, lattice; dim = dim, distanceonly = true, pbc = true)
                    #println("Distance: " * string(distance))

                    # Keep the smallest distance
                    if distance < mindistance
                        mindistance = distance
                    end

                end

            end

        end

    end

    # Return result
    return mindistance

end

# A function to calculate the effective gap size for a given RSA simulation
function calculate_effective_gap_size(status, neighbour_shell_list, Ngrids, grids, molecules, lattice; dim = 2)
    
    # Create empty arrays to store results
    effective_gap_sizes = Vector{Float64}(undef, 0)
    free_grid_points = Matrix{Int64}(undef, 2, 0)

    # Step 1: Generate a list of free gridpoints
    free_point = get_all_free_grid_points(status, molecules, Ngrids, grids, lattice, neighbour_shell_list)


    # Step 2: Loop over all free gridpoints to find the closest adsorbate; the shells are used to limit distance calculations   
    # Loop over all grid points
    for grid_id in 1:Ngrids
        for point_id in 1:grids[grid_id].Npoints

            # Is the point free
            if free_point[grid_id][point_id] == false
                continue
            end

            # To simplify plotting: Add this point to a list
            free_grid_points = hcat(free_grid_points, [grid_id, point_id])

            # Find the first shell containing an adsorbate
            occupied_shell = find_first_shell_with_adsorbate(status, Ngrids, grids, lattice, neighbour_shell_list, grid_id, point_id)

            # At this point the closest occupied shell "shell_occupied" is known
            # Use this and the next shell to find all occupied points - one of these points has the closest distance to the free grid point
            mindistance = get_distance_closest_adsorbate(status, Ngrids, grids, lattice, molecules, neighbour_shell_list, grid_id, point_id, occupied_shell; dim = dim)
            
            # Add the minimal distance to the list
            push!(effective_gap_sizes, mindistance)

        end
    end

    # Return the distance list
    return effective_gap_sizes, free_grid_points

end

# A function to plot the convergence of the metrics of the count and area histograms with the number of RSA runs
"""

    plot_count_area_convergence(Nruns, rsa_results, Nmolecules, molecules, lattice)
    plot_count_area_convergence(Nruns, rsa_results, Nmolecules, molecules, lattice; status = true, plotonly = true, properties = ["means", "standarddeviations"], stride = 1, reference = true, resolution = 600)

Create plots showing the convergence of the metrics with the number of RSA simulations.  


# Input
- `Nruns`: Total number of RSA simulations.
- `rsa_results`: A rsa_run_results_struct object.
- `Nmolecules`: Number of present molecule types.
- `molecules`: A molecule_struct object.
- `lattice`: A lattice_struct object.

# Optional input
- `status`: Flag forcing the recalculation of the status based on the stepinfo field.
- `plotonly`: Flag to request the plotted data in addition to the plots.
- `properties`: Vector of the metrics to be plotted. Accepted values are "means", "standarddeviations", "minvalues", "minvalue\\_ids", "maxvalues", "maxvalue\\_ids", "avgvalues", and "avgvalue\\_ids".
- `stride`: Number of runs added between two evaluations of the metrics. The final evaluation always includes all runs.
- `reference`: Flag to add the value obtained with all RSA runs as a dashed horizontal line.
- `resolution`: Resolution of the images controlled by the dpi value.

# Return values
- `plotonly = true (default)`: A matrix of plots. First index indicates the plotted molecule (area or count) while second index follows the order of the requested properties.
- `plotonly = false`: In addition to the matrix of plots, the vector of the evaluated numbers of runs and the plotted values are returned. The values are given as a vector over the data sets, with every element being a matrix of the evaluated numbers of runs in rows and the requested properties in columns. Information are returned in the following order: plots, numbers of runs, values.
"""
function plot_count_area_convergence(Nruns, rsa_results, Nmolecules, molecules, lattice; status = true, plotonly = true, properties = ["means", "standarddeviations"], stride = 1, reference = true, resolution = 600)

    # Get the adsorbate count and the covered area (in %) of every run
    adsorbate_count_per_run, adsorbate_area_per_run = calculate_count_area_per_run(Nruns, rsa_results, Nmolecules, molecules, lattice; status = status)

    # Collect all data sets to be evaluated
    data_sets, data_labels = collect_count_area_data_sets(Nmolecules, adsorbate_count_per_run, adsorbate_area_per_run)
    Nsets = size(data_sets, 1)

    # Map the requested properties to their position within the metrics of a data set
    Nproperties = size(properties, 1)
    property_ids = Vector{Int64}(undef, Nproperties)
    for property_id in 1:Nproperties
        property_ids[property_id] = statistics_property_id(properties[property_id])
    end

    # Define the numbers of runs to be evaluated
    # The last evaluation always includes all runs
    Nruns_values = collect(stride:stride:Nruns)
    if isempty(Nruns_values) || last(Nruns_values) != Nruns
        push!(Nruns_values, Nruns)
    end
    Npoints = size(Nruns_values, 1)

    # Evaluate the metrics for an increasing number of runs
    convergence_data = Vector{Matrix{Float64}}(undef, Nsets)
    for set_id in 1:Nsets

        # Generate the matrix storing the requested properties of this data set
        convergence_data[set_id] = Matrix{Float64}(undef, Npoints, Nproperties)

        # Loop over all numbers of runs
        for point_id in 1:Npoints

            # Derive all metrics of the first N runs
            metrics = calculate_data_set_metrics(@view data_sets[set_id][1:Nruns_values[point_id]])

            # Store the requested properties
            for property_id in 1:Nproperties
                convergence_data[set_id][point_id, property_id] = Float64(metrics[property_ids[property_id]])
            end

        end

    end

    # Plot the convergence: One plot for every combination of data set and property
    convergence_plots = Matrix{Any}(undef, Nsets, Nproperties)
    for set_id in 1:Nsets
        for property_id in 1:Nproperties

            # Plot the property against the number of included runs
            convergence_plots[set_id, property_id] = plot(Nruns_values, convergence_data[set_id][:, property_id], xlabel = "Number of RSA runs", ylabel = statistics_property_labels[property_ids[property_id]] * "\n" * data_labels[set_id], legend = false, width = 2, dpi = resolution)

            # Add the value obtained with all runs as a reference
            if reference == true
                hline!(convergence_plots[set_id, property_id], [convergence_data[set_id][Npoints, property_id]], linestyle = :dash, linecolor = :red, linewidth = resolution/300)
            end

        end
    end

    # Return results
    if plotonly == true
        return convergence_plots
    else
        return convergence_plots, Nruns_values, convergence_data
    end

end

# A function to plot the convergence of properties of a single RSA run
"""

    plot_single_run_convergence(stepinfo, Nmolecules, molecules, lattice)
    plot_single_run_convergence(stepinfo, Nmolecules, molecules, lattice; errorrange = 0.0, startstep = 1, laststep = 0, plotonly = true, resolution = 600)

Create a plot showing the convergence of adsorbate count and covered area of a single RSA simulation with the number of performed RSA steps.

# Input
- `stepinfo`: A stepinfo field of a rsa_run_results_struct.
- `Nmolecules`: Number of present molecule types.
- `molecules`: A molecule_struct object.
- `lattice`: A lattice_struct object.

# Optional input
- `errorrange`: Range in % which is added to and subtracted from the final value of every property. The resulting range is highlighted within the plot. A value of zero (default) requests no range.
- `startstep`: First RSA step shown by the plot.
- `laststep`: Last RSA step shown by the plot. A value of zero (default) requests all steps of the given stepinfo.
- `plotonly`: Flag to request the plotted data in addition to the plot.
- `resolution`: Resolution of the image controlled by the dpi value.

# Return values
- `plotonly = true (default)`: A vector containing the plots for each adsorbate and all adsorbates together.
- `plotonly = false`: In addition to the plots, a vector containing the ploted data as well as a vector stating at which step convergence is reached are returned.
"""
function plot_single_run_convergence(stepinfo, Nmolecules, molecules, lattice; errorrange = 0.0, startstep = 1, laststep = 0, plotonly = true, resolution = 600)

    # Get the number of performed RSA steps
    Nsteps = size(stepinfo, 2)

    # Define the range of steps to be plotted
    # A laststep of zero (default) requests all steps of the given stepinfo
    if laststep ≤ 0 || laststep > Nsteps
        laststep = Nsteps
    end
    if startstep < 1
        startstep = 1
    end
    if startstep > laststep
        println("The first plotted step (" * string(startstep) * ") is larger than the last plotted step (" * string(laststep) * ").")
        error("Analysis Range Error")
    end

    # Get the count and area per step
    count_per_step, area_per_step = calculate_count_area_per_step(Nsteps, stepinfo, Nmolecules, molecules, lattice)

    # Collect all data sets to be plotted
    data_sets, data_labels = collect_count_area_data_sets(Nmolecules, count_per_step, area_per_step; capital = true)
    Nsets = size(data_sets, 1)

    # Plot the convergence: One plot for every data set
    convergence_plots = Vector{Any}(undef, Nsets)
    convergence_steps = zeros(Int64, Nsets)
    for set_id in 1:Nsets

        # The bare plot
        convergence_plot = plot(xlabel = "RSA step", ylabel = data_labels[set_id], legend = false, dpi = resolution)

        # Add an error range if requested
        if errorrange > 0.0
            
            # Define the range based on the final value of the property
            finalvalue = data_sets[set_id][Nsteps]
            lowervalue = finalvalue * (1.0 - errorrange/100)
            uppervalue = finalvalue * (1.0 + errorrange/100)

            # Highlight the range over the complete range of plotted steps
            plot!(convergence_plot, [startstep, laststep], [lowervalue, lowervalue], fillrange = [uppervalue, uppervalue], fillalpha = 0.2, fillcolor = :grey, linewidth = 0)

            # Get the first step from which the property stays within the range
            convergence_steps[set_id] = find_convergence_step(data_sets[set_id], lowervalue, uppervalue)

        end

        # Plot the data
        plot!(convergence_plot, [startstep:laststep], data_sets[set_id][startstep:laststep], width = 2)
        convergence_plots[set_id] = convergence_plot

    end

    # Return results
    if plotonly == true
        return convergence_plots
    else
        return convergence_plots, data_sets, convergence_steps
    end

end