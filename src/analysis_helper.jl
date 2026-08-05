#
# General head section of any file within this module 
#

#
# Using statements
#

#
# Include statements
#

#
# Definition of global variables
#

#
# Specific section for this file
#

# A helper function to define the y-axes label of histograms
function histogram_label(normalized)

    if normalized == true
        return "Normalized count"
    else
        return "Count"
    end

end

# A function to calculate the scaling factor for gaussian distributions so that they fit to a histogram
# The area of the visible bins (count x width) has to fit the gaussian distribution (usually normalized to 1)
function gaussian_scaling_factor(data_vector, bins)

    # Count all values within the range of the bins
    Nbinned = 0
    for value in data_vector
        if value ≥ first(bins) && value ≤ last(bins)
            Nbinned += 1
        end
    end

    # Return the scaling factor: number of shown values times the bin width
    return Float64(Nbinned) * Float64(step(bins))

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

# A function to derive the number of adsorbates and the covered surface area (in %) of every RSA run
# Both matrices are returned with the molecule types in rows and the RSA runs in columns
function calculate_count_area_per_run(Nruns, rsa_results, Nmolecules, molecules, lattice; status = true)

    # If the status is not present for all RSA runs recalculate the status
    if status == false
        reduce_rsa_allrun_info!(Nruns, rsa_results)
    end

    # Count the number of adsorbates per run
    adsorbate_count_per_run = count_adsorbates_per_rsa_run(Nruns, Nmolecules, rsa_results)

    # Calculate the are per molecule
    molecules_area = calculate_surface_area_molecules(Nmolecules, molecules)

    # Calculate the surface area
    surface_area = norm(cross(lattice.transcellvectors[:,1], lattice.transcellvectors[:,2]))

    # Generate the covered area per molecule per run (in %)
    adsorbate_area_per_run = Matrix{Float64}(undef, Nmolecules, Nruns)
    for run_id in 1:Nruns
        adsorbate_area_per_run[:, run_id] = adsorbate_count_per_run[:, run_id] .* molecules_area / surface_area * 100
    end

    # Return results
    return adsorbate_count_per_run, adsorbate_area_per_run

end

# A function to derive the number of adsorbates and the covered surface area (in %) of every step of a single RSA run
# Both matrices are returned with the molecule types in rows and the RSA runs in columns
function calculate_count_area_per_step(Nsteps, stepinfo, Nmolecules, molecules, lattice)

    # Count the number of adsorbates of every molecule type after every step
    count_per_step = count_adsorbates_per_step(stepinfo, Nmolecules)

    # Calculate the are per molecule
    molecules_area = calculate_surface_area_molecules(Nmolecules, molecules)

    # Calculate the surface area
    surface_area = norm(cross(lattice.transcellvectors[:,1], lattice.transcellvectors[:,2]))

    # Generate the covered area per molecule after every step (in %)
    area_per_step = Matrix{Float64}(undef, Nmolecules, Nsteps)
    for step_id in 1:Nsteps
        area_per_step[:, step_id] = count_per_step[:, step_id] .* molecules_area / surface_area * 100
    end

    # Return results
    count_per_step, area_per_step

end

# Function to derive gaussian distribution over a given vector of numbers
# Returns mean and standard deviation of the distribution as float
function gaussian_distribution(data_vector; zerogaussian = false)
    
    # Number of values
    Nvalues = size(data_vector,1)

    # Calculate the mean value
    if zerogaussian == false
        mean = sum(data_vector) / Nvalues
    else
        mean = 0.0
    end

    # Calculate the standard deviation
    tmp_values = Float64.(deepcopy(data_vector))
    tmp_values .-= mean
    standarddeviation = sqrt(sum(abs2, tmp_values) / Nvalues)

    return mean, standarddeviation

end

# A function to derive all metrics of a data set (mean, standarddeviation, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id)
function calculate_data_set_metrics(data_vector; zerogaussian = false)

    # Derive properties of a normal distribution
    mean, standarddeviation = gaussian_distribution(data_vector, zerogaussian = zerogaussian)

    # Get the smallest and the largest value
    minvalue, minvalue_id = findmin(data_vector)
    maxvalue, maxvalue_id = findmax(data_vector)

    # Get the value closest to the mean value
    avgvalue, avgvalue_id = findmin(abs.(data_vector .- mean))
    avgvalue = data_vector[avgvalue_id]

    # Return all metrics
    return mean, standarddeviation, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id

end

# A function to map a requested metric of a data set to its position within the metrics
function statistics_property_id(property)

    # Find the requested property
    property_id = findfirst(x -> x == property, statistics_property_names)

    # Stop in case the property is unknown
    if property_id === nothing
        println("Unknown property: " * string(property))
        println("Accepted properties: " * join(statistics_property_names, ", "))
        error("Analysis Property Error")
    end

    # Return the position of the property
    return property_id

end

# A function to collect all data sets evaluated by the count and area histograms
# Count and area pairs for every molecule type followed by the total count and the total covered area
function collect_count_area_data_sets(Nmolecules, adsorbate_count_per_run, adsorbate_area_per_run; capital = false)

    # Generate the vectors storing the data sets and their labels
    data_sets = Vector{Vector{Float64}}(undef, (Nmolecules * 2) + 2)
    data_labels = Vector{String}(undef, (Nmolecules * 2) + 2)

    # Add the count and the covered area of every molecule type
    set_id = 0
    for molecule_id in 1:Nmolecules
        set_id += 1
        data_sets[set_id] = Float64.(adsorbate_count_per_run[molecule_id,:])
        if capital == false
            data_labels[set_id] = "adsorbate count - molecule " * string(molecule_id)
        else
            data_labels[set_id] = "Adsorbate count - molecule " * string(molecule_id)
        end
        set_id += 1
        data_sets[set_id] = adsorbate_area_per_run[molecule_id,:]
        if capital == false
            data_labels[set_id] = "covered area in % - molecule " * string(molecule_id)
        else
            data_labels[set_id] = "Covered area in % - molecule " * string(molecule_id)
        end
    end

    # Add the total count and the total covered area
    set_id += 1
    data_sets[set_id] = Float64.(vec(sum(adsorbate_count_per_run, dims = 1)))
    if capital == false
        data_labels[set_id] = "adsorbate count"
    else
        data_labels[set_id] = "Adsorbate count"
    end
    set_id += 1
    data_sets[set_id] = vec(sum(adsorbate_area_per_run, dims = 1))
    if capital == false
        data_labels[set_id] = "covered area in %"
    else
        data_labels[set_id] = "Covered area in %"
    end

    # Return results
    return data_sets, data_labels

end

# A function to update a status matrix based on a single step of a stepinfo matrix
# Returns the updated number of adsorbates stored within the status matrix
function update_status_by_stepinfo!(status, realsize, stepinfo, step_id)

    # Get the information for this step
    selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2 = @view stepinfo[6:12, step_id]
        
    # Update the status matrix
    if selected_event_type == 1
        realsize += 1
        status[1:4,realsize] = [selected_molecule, selected_grid_type, selected_grid_point, selected_event]
    elseif selected_event_type == 2
        change_column = findfirst_column(@view(status[:,1:realsize]), [selected_molecule, selected_grid_type, selected_grid_point], 3)
        status[4,change_column] = selected_event
    elseif selected_event_type == 3
        change_column = findfirst_column(@view(status[:,1:realsize]), [selected_molecule, selected_grid_type, selected_grid_point], 3)
        status[2:4,change_column] = [selected_subevent, selected_event, selected_event_2]
    elseif selected_event_type == 4
        change_column = findfirst_column(@view(status[:,1:realsize]), [selected_molecule, selected_grid_type, selected_grid_point], 3)
        status[1:4,change_column] = [selected_subevent, selected_grid_type, selected_grid_point, selected_event]
    end

    # Return the number of adsorbates
    return realsize

end

# A function to reduce the run information into a final status matrix
function reduce_rsa_run_info(stepinfo)
    
    # Generate the empty matrix
    maxsize = size(stepinfo, 2)
    realsize = 0
    reduced_info = Matrix{Int64}(undef, 4, maxsize)

    # Update the matrix based on every performed rsa step
    for info_id in axes(stepinfo, 2)

        realsize = update_status_by_stepinfo!(reduced_info, realsize, stepinfo, info_id)

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

# A function to derive the plot attributes controlling the white space around the simulation cell
# Returns the margin, the framestyle, and the ticks used for the plot
function simulation_cell_attributes(withmargins)

    if withmargins == true
        return 1 * Plots.mm, :axes, :auto
    else
        return -2 * Plots.mm, :none, false
    end

end

# A function to count the number of adsorbates of every molecule type after every step of a single RSA run
function count_adsorbates_per_step(stepinfo, Nmolecules)

    # Generate the matrix storing the counts of every step
    Nsteps = size(stepinfo, 2)
    count_per_step = Matrix{Int64}(undef, Nmolecules, Nsteps)

    # Count the adsorbates based on every performed rsa step
    current_count = zeros(Int64, Nmolecules)
    for step_id in axes(stepinfo, 2)

        # Get the information for this step
        selected_molecule, selected_event_type, selected_subevent = @view stepinfo[8:10, step_id]

        # Update the current counts
        # An adsorption adds a new adsorbate of the selected molecule type
        # A conversion replaces an adsorbate by an adsorbate of the molecule type given by the subevent
        # Rotations and diffusions do not change the number of adsorbates
        if selected_event_type == 1
            current_count[selected_molecule] += 1
        elseif selected_event_type == 4
            current_count[selected_molecule] -= 1
            current_count[selected_subevent] += 1
        end

        # Store the counts of this step
        count_per_step[:, step_id] = current_count

    end

    # Return results
    return count_per_step

end

# A function to find the first step from which a value stays within a given range
# A returned step of zero indicates that the value does not stay within the range
function find_convergence_step(values, lowervalue, uppervalue)

    # Find the last value outside of the range
    last_outside = 0
    for point_id in eachindex(values)
        if values[point_id] < lowervalue || values[point_id] > uppervalue
            last_outside = point_id
        end
    end

    # The property stays within the range starting with the following step
    if last_outside == size(values, 1)
        return 0
    else
        return last_outside + 1
    end

end