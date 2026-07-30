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

# A function to derive all metrics of a data set (mean, variance, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id)
function calculate_data_set_metrics(data_vector; zerogaussian = false)

    # Derive properties of a normal distribution
    mean, variance = gaussian_distribution(data_vector, zerogaussian = zerogaussian)

    # Get the smallest and the largest value
    minvalue, minvalue_id = findmin(data_vector)
    maxvalue, maxvalue_id = findmax(data_vector)

    # Get the value closest to the mean value
    avgvalue, avgvalue_id = findmin(abs.(data_vector .- mean))
    avgvalue = data_vector[avgvalue_id]

    # Return all metrics
    return mean, variance, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id

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
function collect_count_area_data_sets(Nmolecules, adsorbate_count_per_run, adsorbate_area_per_run)

    # Generate the vectors storing the data sets and their labels
    data_sets = Vector{Vector{Float64}}(undef, (Nmolecules * 2) + 2)
    data_labels = Vector{String}(undef, (Nmolecules * 2) + 2)

    # Add the count and the covered area of every molecule type
    set_id = 0
    for molecule_id in 1:Nmolecules
        set_id += 1
        data_sets[set_id] = Float64.(adsorbate_count_per_run[molecule_id,:])
        data_labels[set_id] = "adsorbate count - molecule " * string(molecule_id)
        set_id += 1
        data_sets[set_id] = adsorbate_area_per_run[molecule_id,:]
        data_labels[set_id] = "covered area in % - molecule " * string(molecule_id)
    end

    # Add the total count and the total covered area
    set_id += 1
    data_sets[set_id] = Float64.(vec(sum(adsorbate_count_per_run, dims = 1)))
    data_labels[set_id] = "adsorbate count"
    set_id += 1
    data_sets[set_id] = vec(sum(adsorbate_area_per_run, dims = 1))
    data_labels[set_id] = "covered area in %"

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
        change_column = findfirst_column(status, [selected_molecule, selected_grid_type, selected_grid_point], 3)
        status[4,change_column] = selected_event
    elseif selected_event_type == 3
        change_column = findfirst_column(status, [selected_molecule, selected_grid_type, selected_grid_point], 3)
        status[2:4,change_column] = [selected_subevent, selected_event, selected_event_2]
    elseif selected_event_type == 4
        change_column = findfirst_column(status, [selected_molecule, selected_grid_type, selected_grid_point], 3)
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