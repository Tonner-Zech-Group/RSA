#
# General head section of any file within this module 
#

#
# Using statements
#
using ProgressBars


#
# Include statements
#

#
# Export statements
#
export plot_RSA_run
export plot_single_molecule
export animate_RSA_run
export plot_count_area_histograms
export plot_effective_gap_size
export gif
export savefig

#
# Definition of global variables
#

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

# A function to plot all molecules on their selected gridpoints
"""

    plot_RSA_run(status, Ngrids, grids, Nmolecules, molecules, lattice)
    plot_RSA_run(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, boundary_cells = 1, silent=true)

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

# Return values
- A plots object of the covered surface.
"""
function plot_RSA_run(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, boundary_cells = 1, silent=true)

    # Define the image resolution
    x_axis_size = lattice.transvectors[1,1] + lattice.transvectors[2,1]
    y_axis_size = lattice.transvectors[1,2] + lattice.transvectors[2,2]
    x_axis_resolution = x_axis_size * pixel_per_angstrom
    y_axis_resolution = y_axis_size * pixel_per_angstrom
    
    if silent != true
        println("Resolution: " * string(x_axis_resolution) * " x " * string(y_axis_resolution))
    end

    # Define color of grid points
    grid_palette = palette(:darktest, Ngrids)

    # Plot the grids
    final_plot = 0
    for grid_id in 1:Ngrids
        if grid_id == 1
            final_plot = scatter(grids[grid_id].points[1,:],grids[grid_id].points[2,:], markersize=pixel_per_angstrom/10, legend=false, showaxis=false, grid=false, size=(x_axis_resolution, y_axis_resolution), xlims=(0, x_axis_size), ylims=(0, y_axis_size), color = grid_palette[grid_id], widen = false)
        else
            scatter!(grids[grid_id].points[1,:],grids[grid_id].points[2,:], markersize=pixel_per_angstrom/10, xlims=(0, x_axis_size), ylims=(0, y_axis_size), color = grid_palette[grid_id])
        end
    end
    
    #println("Number of grids: " * string(Ngrids))   
    
    # Define color of molecule
    molecule_palette = palette(:darktest, Nmolecules)    

    # Plot every molecule
    for element_id in axes(status, 2)

        # Get the displaced coordinates and elements
        molec_coords, displaced_coordinates, molec_elements, gridpoint_coords, molecule_id, grid_id, point_id = displaced_molecule_status(status, element_id, molecules, grids, lattice)
 
        # Define the marker size
        markersize_vector = atomic_information[molec_elements[:], 3] * pixel_per_angstrom 

        # Add to plot
        scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])

        # Check whether the adsorbate was within the boundary cells
        transx, transy = grids[grid_id].mapping[2:3,point_id]
        # Increase by one to adapt the scale to "1 to Ncells"
        transx += 1
        transy += 1
        if transx ≤ boundary_cells
            if transy ≤ boundary_cells
                gridpoint_coords_moved = gridpoint_coords + lattice.transvectors[1,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
                gridpoint_coords_moved = gridpoint_coords + lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
                gridpoint_coords_moved = gridpoint_coords + lattice.transvectors[1,:] + lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
            elseif (lattice.Ncelly - boundary_cells) < transy
                gridpoint_coords_moved = gridpoint_coords + lattice.transvectors[1,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
                gridpoint_coords_moved = gridpoint_coords - lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
                gridpoint_coords_moved = gridpoint_coords + lattice.transvectors[1,:] - lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
            else
                gridpoint_coords_moved = gridpoint_coords + lattice.transvectors[1,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
            end
        elseif (lattice.Ncellx - boundary_cells) < transx
            if transy ≤ boundary_cells
                gridpoint_coords_moved = gridpoint_coords - lattice.transvectors[1,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
                gridpoint_coords_moved = gridpoint_coords + lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
                gridpoint_coords_moved = gridpoint_coords - lattice.transvectors[1,:] + lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id]) 
            elseif (lattice.Ncelly - boundary_cells) < transy
                gridpoint_coords_moved = gridpoint_coords - lattice.transvectors[1,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
                gridpoint_coords_moved = gridpoint_coords - lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
                gridpoint_coords_moved = gridpoint_coords - lattice.transvectors[1,:] - lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
            else
                gridpoint_coords_moved = gridpoint_coords - lattice.transvectors[1,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
            end
        else
            if transy ≤ boundary_cells
                gridpoint_coords_moved = gridpoint_coords + lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
            elseif (lattice.Ncelly - boundary_cells) < transy
                gridpoint_coords_moved = gridpoint_coords - lattice.transvectors[2,:]
                displaced_coordinates = move_structure_to_point(molec_coords, gridpoint_coords_moved, [0.0, 0.0, 0.0], lattice.dimension)
                scatter!(displaced_coordinates[1,:],displaced_coordinates[2,:], markersize = markersize_vector, markerstrokewidth = 0, xlims = (0, x_axis_size), ylims = (0, y_axis_size), widen = false, color = molecule_palette[molecule_id])
            end
        end

    end

    # Return final plot
    return final_plot

end

# A function to create an animation of an rsa run
"""

    animate_RSA_run(stepinfo, Ngrids, grids, Nmolecules, molecules, lattice)
    animate_RSA_run(stepinfo, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, boundary_cells = 1)

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

# Return values
- A plots object containing the animation of the RSA simulation.

# Hints
- Generation of large animations is extremely slow.
- Only reasonable to use for the initial ~1000 steps (with stepinfo[:,1:1000]).
"""
function animate_RSA_run(stepinfo, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, boundary_cells = 1)

    # Throw a warning in case the resolution is getting to large
    if pixel_per_angstrom > 10
        println("Hint to the user: Your pixel_per_angstrom is larger then 10. Expect the generation of the animation to slow down.")
    end
    if pixel_per_angstrom > 100
        println("Hint to the user: Your pixel_per_angstrom has an extremely large value! The generation of the animation will be slow and you will get a large animation file. Consider reducing this value.")
    end

    # Create the animation object
    anim = Animation()

    # Preallocate matrices
    realsize = 0
    Nframes = size(stepinfo,2)
    status = Matrix{Int64}(undef, 4, Nframes)

    # Create the frames
    for frame_id in ProgressBar(1:Nframes)

        # Get the information for this step
        selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2 = @view stepinfo[6:12, frame_id]
        
        # Update the status matrix
        if selected_event_type == 1
            realsize += 1
            status[1:4,realsize] = [selected_molecule, selected_grid_type, selected_grid_point, selected_event]
        elseif selected_event_type == 2
            change_column = findfirst_column(status, [selected_molecule, selected_grid_type, selected_grid_point], 3)
            status[4,change_column] = selected_event
        elseif selected_event_type == 3
            change_column = findfirst_column(status, [selected_molecule, selected_grid_type, selected_grid_point], 3)
            status[2:4,change_column] = [selected_subevent, selected_event, selected_event_2]
        elseif selected_event_type == 4
            change_column = findfirst_column(status, [selected_molecule, selected_grid_type, selected_grid_point], 3)
            status[1:4,change_column] = [selected_subevent, selected_grid_type, selected_grid_point, selected_event]
        end

        # Create the frame
        substatus = @view status[1:4,1:realsize]
        newframe = plot_RSA_run(substatus, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, boundary_cells = 1, silent=true)

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

# A function to reduce the run information into a final status matrix
function reduce_rsa_run_info(stepinfo)
    
    # Generate the empty matrix
    maxsize = size(stepinfo, 2)
    realsize = 0
    reduced_info = Matrix{Int64}(undef, 4, maxsize)

    # Update the matrix based on every performed rsa step
    for info_id in axes(stepinfo, 2)

        # Get the information for this step
        selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2 = @view stepinfo[6:12, info_id]

        # Update the reduced_info matrix
        if selected_event_type == 1
            realsize += 1
            reduced_info[1:4,realsize] = [selected_molecule, selected_grid_type, selected_grid_point, selected_event]
        elseif selected_event_type == 2
            change_column = findfirst_column(reduced_info, [selected_molecule, selected_grid_type, selected_grid_point], 3)
            reduced_info[4,change_column] = selected_event
        elseif selected_event_type == 3
            change_column = findfirst_column(reduced_info, [selected_molecule, selected_grid_type, selected_grid_point], 3)
            reduced_info[2:4,change_column] = [selected_subevent, selected_event, selected_event_2]
        elseif selected_event_type == 4
            change_column = findfirst_column(reduced_info, [selected_molecule, selected_grid_type, selected_grid_point], 3)
            reduced_info[1:4,change_column] = [selected_subevent, selected_grid_type, selected_grid_point, selected_event]
        end

    end

    # Return the result
    return reduced_info[:,1:realsize]

end

# A function to reduce the run information of all runs into a finals status matrix (using a rsa_run_results_struct object)
function reduce_rsa_allrun_info!(Nruns, rsa_results)

    # Loop over all runs
    for run_id in 1:Nruns

        # Get the status matrix
        rsa_results[run_id].status = reduce_rsa_run_info(rsa_results[run_id].stepinfo)

    end

end

# A function to count the number of adsorbates per run
function count_adsorbates_per_rsa_run(Nruns, Nmolecules, rsa_results)
    
    # Allocate the count vector
    run_adsorbate_count = Matrix{Int64}(undef, Nmolecules, Nruns)

    # Loop over every RSA run
    for run_id in 1:Nruns

        # Generate the adsorbate count vector
         adsorbate_count = zeros(Int64, Nmolecules)

        # Count the occurence of any molecule in the status matrix
        for molecule_id in 1:Nmolecules
            adsorbate_count[molecule_id] = size(findall(x -> x == molecule_id, @view rsa_results[run_id].status[1,:]),1)
        end

        # Add to the final matrix
        run_adsorbate_count[:, run_id] = adsorbate_count

    end

    # Return results
    return run_adsorbate_count

end

# A function to calculate the surface area for every molecule of a molecules object
function calculate_surface_area_molecules(Nmolecules, molecules; resolution = 0.01)

    # Create the vector
    molecules_area = Vector{Float64}(undef, Nmolecules)

    # Loop over all molecules
    for molecule_id in 1:Nmolecules
        molecules_area[molecule_id] = calculate_surface_area(molecules[molecule_id].elements, molecules[molecule_id].coordinates, resolution = resolution)
    end

    # Return results
    return molecules_area

end

# Function to calculate the covered surface based on the molecule
function calculate_surface_area(molecule_elements, molecule_coords; resolution = 0.01)
    
    # Move the molecule to the center of the coordinate system
    origin = [0.0, 0.0, 0.0]
    centroid = calculate_centroid(molecule_coords)
    dimension = size(molecule_coords,1)
    moved_coords = move_structure_to_point(molecule_coords, origin, centroid, dimension)

    # Get the radius of the molecule
    Natoms = size(molecule_coords,2)
    maxradius, radii = get_largest_vdW_distance_to_point(Natoms, molecule_elements, moved_coords, origin, "2D") 

    # Loop over a grid and count the gridpoints covered by the molecule
    # We add 2% to the radius to be on the save side
    count = 0
    for x_value in range(-1.02 * maxradius, 1.02 * maxradius, step = resolution)
        for y_value in range(-1.02 * maxradius, 1.02 * maxradius, step = resolution)
            point = [x_value, y_value]
            covered = point_covered_by_vdW_radii_2D(molecule_elements, moved_coords, point)
            if covered == true
                count += 1
            end

        end
    end

    # Convert count to area
    area = Float64(count) * resolution^2

    # Return result
    return area

end

# Function to derive gaussian distribution over a given vector of numbers
# Returns mean and variance of the distribution as float
function gaussian_distribution(data_vector; zerogaussian = false)
    
    # Number of values
    Nvalues = size(data_vector,1)

    # Calculate the mean value
    if zerogaussian == false
        mean = sum(data_vector) / Nvalues
    else
        mean = 0.0
    end

    # Calculate the variance
    tmp_values = Float64.(deepcopy(data_vector))
    tmp_values .-= mean
    variance = sqrt(sum(abs2, tmp_values) / Nvalues)

    return mean, variance

end

# A function to plot a histogram based on a given vector of numbers
function plot_histogram(data_vector; labelx = "Value", labely = "Frequentness", stepsize = 1.0, resolution = 600, distribution = "none", plotonly = true, threshold = 0.0)

    # Derive properties of a normal distribution
    if distribution == "truncated"
        mean, variance = gaussian_distribution(data_vector, zerogaussian = true)
    elseif distribution == "gaussian"
        mean, variance = gaussian_distribution(data_vector, zerogaussian = false)
    else
        mean, variance = 0.0, 0.0
    end


    # Get more data
    #Nvalues = size(data_vector, 1)
    minvalue, minvalue_id = findmin(data_vector)
    maxvalue, maxvalue_id = findmax(data_vector)
    if distribution == "truncated" || distribution == "gaussian"
        avgvalue, avgvalue_id = findmin(abs.(data_vector .- mean))
        avgvalue = data_vector[avgvalue_id]
    else
        avgvalue, avgvalue_id = 0.0, 0
    end

    # Plot the histogram
    histo = histogram(data_vector, xlabel = labelx, ylabel = labely, normalize = :pdf, xlims=(minvalue, maxvalue), legend=false, bins = range(minvalue, maxvalue, step = stepsize), dpi = resolution)

    # Add the normal distribution
    if distribution == "gaussian"
        x = range(minvalue, maxvalue, step = stepsize/10)
        y = @. 1/(variance * sqrt(2*π)) * exp(-0.5 * (x - mean)^2 / variance^2) 
        plot!(x,y, width = 4, lc = "red")
    elseif distribution == "truncated"
        x = range(minvalue, maxvalue, step = stepsize/10)
        y = @. 2 * 1/(variance * sqrt(2*π)) * exp(-0.5 * (x - mean)^2 / variance^2) 
        plot!(x,y, width = 4, lc = "red")
    end

    # Add a vertical line if requested
    if threshold > 0.0
        vline!([threshold], linestyle = :dash, linecolor = :red, linewidth = resolution/200)
    end

    # Return all results
    if plotonly == true
        return histo
    else
        return histo, mean, variance, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id
    end
    
end

# A function to create histogram plots for adsorbate count and surface coverage
"""

    plot_count_area_histograms(Nruns, rsa_results, Nmolecules, molecules, lattice)
    plot_count_area_histograms(Nruns, rsa_results, Nmolecules, molecules, lattice; status = true, plotonly = true, count = 1.0, area = 1.0)

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

# Return values
- `plotonly = true (default)`: A vector of histograms is returned. Count and covered surface is contained pairwise for every molecule type while the second last element contains the total adsorbate count and the last element the total covered area.
- `plotonly = false`: In addition to the histogram vector, vectors storing the mean values, variance, min and max values, as well as the simulation closest to the mean value are returned. In the following order: histograms, means, variances, minvalues, minvalue_ids, maxvalues, maxvalue_ids, averagevalues, averagevalue_ids.
"""
function plot_count_area_histograms(Nruns, rsa_results, Nmolecules, molecules, lattice; status = true, plotonly = true, count = 1.0, area = 1.0)

    # If the status is not present for all RSA runs recalculate the status
    if status == false
        reduce_rsa_allrun_info!(Nruns, rsa_results)
    end

    # Count the number of adsorbates per run
    adsorbate_count_per_run = count_adsorbates_per_rsa_run(Nruns, Nmolecules, rsa_results)

    # Calculate the are per molecule
    molecules_area = calculate_surface_area_molecules(Nmolecules, molecules)

    # Calculate the surface area
    surface_area = norm(cross(lattice.transvectors[1,:], lattice.transvectors[2,:]))

    # Generate the covered area per molecule per run (in %)
    adsorbate_area_per_run = Matrix{Float64}(undef, Nmolecules, Nruns)
    for run_id in 1:Nruns
        adsorbate_area_per_run[:, run_id] = adsorbate_count_per_run[:, run_id] .* molecules_area / surface_area * 100
    end

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
        variances = Vector{Float64}(undef, (Nmolecules * 2) + 2)
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
            histos[plot_id] = plot_histogram(adsorbate_count_per_run[molecule_id,:]; labelx = "Adsorbate count - molecule " * string(molecule_id), labely = "Normalized Frequentness", stepsize = count, resolution = 600, plotonly = true, distribution = "gaussian")
            plot_id += 1
            histos[plot_id] = plot_histogram(adsorbate_area_per_run[molecule_id,:]; labelx = "Covered area in % - molecule " * string(molecule_id), labely = "Normalized Frequentness", stepsize = area, resolution = 600, plotonly = true, distribution = "gaussian")
        end
        
        plot_id += 1
        histos[plot_id] = plot_histogram(total_adsorbate_count[:]; labelx = "Adsorbate count", labely = "Normalized Frequentness", stepsize = count, resolution = 600, plotonly = true, distribution = "gaussian")
        plot_id += 1
        histos[plot_id] = plot_histogram(total_area[:]; labelx = "Covered area in %", labely = "Normalized Frequentness", stepsize = area, resolution = 600, plotonly = true, distribution = "gaussian")

    else

        plot_id = 0
        for molecule_id in 1:Nmolecules
            plot_id += 1
            histos[plot_id], means[plot_id], variances[plot_id], minvalues[plot_id], minvalue_ids[plot_id], maxvalues[plot_id], maxvalue_ids[plot_id], avgvalues[plot_id], avgvalue_ids[plot_id] = plot_histogram(adsorbate_count_per_run[molecule_id,:]; labelx = "Adsorbate count - molecule " * string(molecule_id), labely = "Normalized Frequentness", stepsize = count, resolution = 600, plotonly = false, distribution = "gaussian")
            plot_id += 1
            histos[plot_id], means[plot_id], variances[plot_id], minvalues[plot_id], minvalue_ids[plot_id], maxvalues[plot_id], maxvalue_ids[plot_id], avgvalues[plot_id], avgvalue_ids[plot_id] = plot_histogram(adsorbate_area_per_run[molecule_id,:]; labelx = "Covered area in % - molecule " * string(molecule_id), labely = "Normalized Frequentness", stepsize = area, resolution = 600, plotonly = false, distribution = "gaussian")
        end
        
        plot_id += 1
        histos[plot_id], means[plot_id], variances[plot_id], minvalues[plot_id], minvalue_ids[plot_id], maxvalues[plot_id], maxvalue_ids[plot_id], avgvalues[plot_id], avgvalue_ids[plot_id] = plot_histogram(total_adsorbate_count[:]; labelx = "Adsorbate count", labely = "Normalized Frequentness", stepsize = count, resolution = 600, plotonly = false, distribution = "gaussian")
        plot_id += 1
        histos[plot_id], means[plot_id], variances[plot_id], minvalues[plot_id], minvalue_ids[plot_id], maxvalues[plot_id], maxvalue_ids[plot_id], avgvalues[plot_id], avgvalue_ids[plot_id] = plot_histogram(total_area[:]; labelx = "Covered area in %", labely = "Normalized Frequentness", stepsize = area, resolution = 600, plotonly = false, distribution = "gaussian")

    end

    # Return results
    if plotonly == true
        return histos
    else
        return histos, means, variances, minvalues, minvalue_ids, maxvalues, maxvalue_ids, avgvalues, avgvalue_ids
    end

end

# A function generating shell lists for every unique grid point
# The tolerance is multiplied with the largest molecule radius to define the shell size
function create_neighbour_shell_lists(Ngrids, grids, Nmolecules, molecules, lattice; tolerance = 1.05)

    # Create the vector of vectors to store the shells
    # Grid --> Unique point --> Grid --> Shells --> Points
    neighbour_shell_list = Vector{Vector{Vector{Vector{Vector{Int64}}}}}(undef, Ngrids)
    for grid_A_id in 1:Ngrids
        neighbour_shell_list[grid_A_id] = Vector{Vector{Vector{Vector{Int64}}}}(undef, grids[grid_A_id].Nuniquepoints)
        for unique_id in 1:grids[grid_A_id].Nuniquepoints
            neighbour_shell_list[grid_A_id][unique_id] = Vector{Vector{Vector{Int64}}}(undef, Ngrids)
        end
    end

    # Get the largest radius of any molecule
    radius = 0.0
    for molecule_id in 1:Nmolecules
        if molecules[molecule_id].maxradius > radius
            radius = molecules[molecule_id].maxradius
        end
    end

    # Add the tolerance to this radius
    radius *= tolerance

    # Loop over all grids
    for grid_A_id in 1:Ngrids

        # Loop over all unique points
        for unique_id in 1:grids[grid_A_id].Nuniquepoints

            # Get the current unique point
            upoint = @view grids[grid_A_id].uniquepoints[:, unique_id]

            # Loop over all grids
            for grid_B_id in 1:Ngrids

                # Calculate the distance vectors between the unique point and all points of this grid
                distance_vectors =  grids[grid_B_id].points .- upoint

                # Correct for PBC
                corrected_distance_vectors = apply_pbc_to_coordinates(distance_vectors, lattice.transvectors, lattice.inversevectors)

                # Get the distances
                distances = vec(sqrt.(sum(abs2, corrected_distance_vectors, dims = 1)))

                # Sort the distances for increasing distances
                sorted_index = sortperm(distances)

                # Create the sorted distances
                sorted_distances = distances[sorted_index]

                # Define empty vector to store shells
                shells_vector = Vector{Vector{Int64}}(undef, 0)

                # Loop over the sorted distances and add them to shells based on the distance
                # In case both grids are identical: Skip the first entry as this is always the point itself (with a distance of zero)
                if grid_A_id == grid_B_id
                    start_id = 2
                else
                    start_id = 1
                end
                shell = 1
                for distance_id in 2:grids[grid_B_id].Npoints
                    if sorted_distances[distance_id] > shell * radius
                        push!(shells_vector, sorted_index[start_id:distance_id-1])
                        start_id = distance_id
                        shell += 1
                    end
                end
                
                # Add the last block to the shell vector
                push!(shells_vector, sorted_index[start_id:grids[grid_B_id].Npoints])

                # Add the shell vector to the final vector
                neighbour_shell_list[grid_A_id][unique_id][grid_B_id] = shells_vector

            end

        end

    end

    # Return the result
    return neighbour_shell_list

end

# A function to plot the effective gap size
"""

    plot_effective_gap_size(status, Ngrids, grids, Nmolecules, molecules, lattice)
    plot_effective_gap_size(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, gapsonly = false, withstroke = true, plotonly = true)

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


# Return values
- `plotonly = true (default)`: Returns the histogram showing the frequency of all gap sizes and a plots object for the visualization of the gaps in the following order: histogram, plot.
- `plotonly = false`: In addition to the default case, a vector containing the obtained effective gap sizes as well as a vector of the corresponding free grid point are returned. Information are returned in the following order: histogram, plot, gap sizes, free grid points.
"""
function plot_effective_gap_size(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = 10.0, gapsonly = false, withstroke = true, plotonly = true, stepsize = 0.2, threshold = 0.0)

    # Get the shells of every grid point
    neighbour_shell_list = create_neighbour_shell_lists(Ngrids, grids, Nmolecules, molecules, lattice)

    # Get the effective gap sizes
    effective_gap_sizes, free_grid_points = calculate_effective_gap_size(status, neighbour_shell_list, Ngrids, grids, molecules, lattice; dim = 2)
    
    # Create the histogram of the effective gap sizes
    gaps_histogram = plot_histogram(effective_gap_sizes; labelx = "Effective Gap Size in Å", labely = "Normalized Frequentness", stepsize = stepsize, threshold = threshold)

    # Get a plot of all adsorbates
    if gapsonly == false
        plot_gaps = plot_RSA_run(status, Ngrids, grids, Nmolecules, molecules, lattice; pixel_per_angstrom = pixel_per_angstrom, boundary_cells = 1, silent=true)
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
        x_axis_size = lattice.transvectors[1,1] + lattice.transvectors[2,1]
        y_axis_size = lattice.transvectors[1,2] + lattice.transvectors[2,2]
        x_axis_resolution = x_axis_size * pixel_per_angstrom
        y_axis_resolution = y_axis_size * pixel_per_angstrom
    
        # Define color of grid points
        grid_palette = palette(:darktest, Ngrids)

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
            plot_gaps = scatter(coordinates[1,:],coordinates[2,:], markersize = effective_gap_sizes * pixel_per_angstrom, markerstrokewidth = 10 / pixel_per_angstrom, legend=false, showaxis=false, grid=false, size=(x_axis_resolution, y_axis_resolution), xlims=(0, x_axis_size), ylims=(0, y_axis_size), widen = false, color = :grey)
        else
            plot_gaps = scatter(coordinates[1,:],coordinates[2,:], markersize = effective_gap_sizes * pixel_per_angstrom, markerstrokewidth = 0, legend=false, showaxis=false, grid=false, size=(x_axis_resolution, y_axis_resolution), xlims=(0, x_axis_size), ylims=(0, y_axis_size), widen = false, color = :grey)
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
