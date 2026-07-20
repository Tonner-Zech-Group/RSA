# Collection of all functions related to calculation steps

#
# General head section of any file within this module 
#

#
# Using statements
#


#
# Export statements
#
 
#
# Definition of global variables
#

#
# Specific section for this file
#

# Function to move a structure defined by a point set to a given point in space
function move_structure_to_point(coordinates, target_point, starting_point, dimension)

    # Calculate the difference to the point
    difference = target_point - starting_point

    # Move by the difference
    displaced_coordinates = coordinates[1:dimension,:] .+ difference

    # Return moved structure
    return displaced_coordinates

end

# Function to calculate the largest distance of a vdw-radius to a given point 
# returns the largest value as well as a vector with each distance
# TODO: To be more general, add the possibility for pbc to this function
function get_largest_vdW_distance_to_point(Natoms, molecule_elements, molecule_coords, point, case="3D")

    # Default
    distance = 0.0
    distances = Vector{Float64}(undef, Natoms)

    # Check all atoms
    for atom_id in axes(molecule_coords,2)

        # Get the difference
        if case == "3D"
            difference = norm(molecule_coords[:,atom_id] - point) + (atomic_information[molecule_elements[atom_id], 3])
        elseif case == "2D"
            difference = norm(molecule_coords[1:2,atom_id] - point[1:2]) + (atomic_information[molecule_elements[atom_id], 3])
        end
        
        # Compare to current max distance
        if difference > distance
            distance = difference
        end
        distances[atom_id] = difference

    end
    
    # Return
    return distance, distances

end

# Function to calculate the smallest distance of a vdw-radius to a given point 
function get_smallest_vdW_distance_to_point(molecule_elements, molecule_coords, point, lattice; dim=3, distanceonly = true, pbc = true)

    # Get the distanes
    distance_vectors =  molecule_coords .- point

    # Correct for periodic boundary conditions
    if pbc == true
        distance_vectors = apply_pbc_to_coordinates(distance_vectors, lattice.transvectors, lattice.inversevectors)
    end

    # Calculate the distances
    distances = vec(sqrt.(sum(abs2, distance_vectors[1:dim,:], dims = 1)))
    
    # Correct by the vdw radii
    distances .-= Float64.(atomic_information[molecule_elements,3])

    # Find the smallest distance
    mindistance, mindistance_id = findmin(distances)

    # Return the results
    if distanceonly == true
        return mindistance
    else
        return mindistance, distances
    end 
    
end

# Function to generate all rotations and store them
# The initial coordinates are centered on the fixpoint (and thereby changed!)
function create_all_rotations_centered(coordinates, Nrotations, rotationvalues, fixpoint, dimension)
    # Defaults
    origin = zeros(dimension)
    AllRotations = Vector{Matrix{Float64}}(undef, Nrotations)

    # Move initial structure to origin
    displaced_coordinates = move_structure_to_point(coordinates, origin, fixpoint, dimension)

    # Create rotation and store them
    Threads.@threads for rot_angle_id in axes(rotationvalues,1)
        # Copy the original coordinates
        molec_1 = deepcopy(displaced_coordinates)

        # Rotate and store
        AllRotations[rot_angle_id] = rotate_structure(molec_1, rotationvalues[rot_angle_id], dimension)
    end

    # Return result
    return AllRotations
end

# Function to check whether a given point is covered by vdW radii
function point_covered_by_vdW_radii_2D(molecule_elements, molecule_coords, point)
    
    # Default
    covered = false

    # Check all atoms
    for ele in axes(molecule_coords, 2)
        
        distance = norm(molecule_coords[1:2,ele] - point[1:2])
        if distance < atomic_information[molecule_elements[ele], 3]
            covered = true
            break
        end

    end

    # Return
    return covered

end 

function point_covered_by_vdW_radii_2D(molecule_elements, molecule_coords, point, lattice; dim = 2)
    
    # Default
    covered = false

    # Check all atoms
    for ele in axes(molecule_coords, 2)
        
        distance_vector = molecule_coords[:,ele] - point
        distance_vector = apply_pbc_to_coordinates(distance_vector, lattice.transvectors, lattice.inversevectors)
        distance = norm(distance_vector[1:dim])

        if distance < atomic_information[molecule_elements[ele], 3]
            covered = true
            break
        end
        
    end

    # Return
    return covered

end

# Function to calculate the centroid of a given set of points
function calculate_centroid(molecule_coords)
    
    # Get number of atoms
    Natoms = size(molecule_coords, 2)
    centroid = sum(molecule_coords, dims=2) / Natoms

    # Return result
    return centroid

end

# Function to calculate the centroid of a limited set of points
function caclulate_partial_centroid(molecule_coords, atom_list)
    
    # Get number of atoms
    Natoms = size(atom_list, 1)
    dimension = size(molecule_coords, 1)

    # Sum over all atoms in the atom list
    centroid = zeros(dimension)

    for ele in atom_list
        centroid += molecule_coords[:,ele]
    end

    # Normalize
    centroid = centroid / Natoms

    # Return result as matrix
    return reshape(centroid, dimension, 1)

end

# Function to clock-wise rotate a set of given points around the z-axis by a given angle in degree
function rotate_structure(coordinates, rotation_angle, dimension)

    # Get rotation matrix, convert degree to rad
    rotation_matrix = calculate_rotation_matrix(rotation_angle*pi/180, dimension)

    # Rotate
    rotated_coordinates = rotation_matrix * coordinates

    # Return
    return rotated_coordinates

end

# Function to calculate a rotation matrix for a rotation around the z-axis by a given angle in rad
function calculate_rotation_matrix(rotation_angle, dimension)

    # Define Matrix
    if dimension == 3
        rotation_matrix = Matrix{Float64}(undef, 3, 3)
        rotation_matrix = [cos(rotation_angle) -1*sin(rotation_angle) 0;
                sin(rotation_angle) cos(rotation_angle) 0;
                0 0 1]
    elseif dimension == 2
        rotation_matrix = Matrix{Float64}(undef, 2, 2)
        rotation_matrix = [cos(rotation_angle) -1*sin(rotation_angle);
                sin(rotation_angle) cos(rotation_angle)]
    end

    # Return
    return rotation_matrix

end

# A function to calculate a release and block matrix based on entries only available in the first or second matrix, respectively.
function delta_matrix(matrix_A, matrix_B)

    # Get the intersect of matrix A and B
    Bcols = Set(eachcol(matrix_B))
    intersect_matrix = matrix_A[:, @views [i for i in axes(matrix_A,2) if matrix_A[:,i] ∈ Bcols]]
    Icols = Set(eachcol(intersect_matrix))
    
    # Get the release and block matrix
    release_matrix = matrix_A[:, @views [i for i in axes(matrix_A,2) if matrix_A[:,i] ∉ Icols]]
    block_matrix = matrix_B[:, @views [i for i in axes(matrix_B,2) if matrix_B[:,i] ∉ Icols]]

    # Return the results
    return release_matrix, block_matrix

end

# Function to map a complete list of affected unique points to the actual points
function map_affected_grid_points_to_actual_points(Affected_Points_Rotations, molecule_id, grid_id, unique_point_id, rotation_id, grids, lattice, reference_transx, reference_transy)

    # Generate a copy of the unique affected points
    mapped_points = deepcopy(Affected_Points_Rotations[molecule_id][grid_id][unique_point_id][rotation_id])

    # Map all unique values to the actual values
    for point_id in axes(mapped_points,2)
        mapping_unique_point, mapping_transx, mapping_transy = grids[mapped_points[2,point_id]].mapping[1:3, mapped_points[3,point_id]]
        new_transx = mapping_transx + reference_transx
        new_transy = mapping_transy + reference_transy
        affected_point = map_translation_to_gridpoint(mapping_unique_point, grids[mapped_points[2,point_id]].Nuniquepoints, new_transx, lattice.Ncellx, new_transy, lattice.Ncelly)
        mapped_points[3,point_id] = affected_point
    end

    # Return the results
    return mapped_points

end

# A function to find the first column in a matrix 
# dim gives the number of elements to be considered in every column
function findfirst_column(matrix, column, dim)
    for column_id in axes(matrix,2)
        if column == @view matrix[1:dim,column_id]
            return column_id
        end
    end
end

# A function to find the first column in a matrix 
# subset gives the elements to be considered in every column
function findfirst_subset(matrix::Matrix{Int64}, column::Vector{Int64}, subset::UnitRange{Int64})
    for column_id in axes(matrix,2)
        if column == @view matrix[subset,column_id]
            return column_id
        end
    end
    return 0
end

# A function to find the last column in a matrix 
# subset gives the elements to be considered in every column
# Keep for debugging purposes, but not used in the current version of the code
function findlast_subset(matrix::Matrix{Int64}, column::Vector{Int64}, subset::UnitRange{Int64})
    for column_id in reverse(axes(matrix,2))
        if column == @view matrix[subset,column_id]
            return column_id
        end
    end
    return 0
end

# A function to check whether a subset is present in a matrix
# subset gives the elements to be considered in every column
function present_column(matrix::Matrix{Int64}, column::Vector{Int64}, subset::UnitRange{Int64})
    for column_id in axes(matrix,2)
        if column == @view matrix[subset,column_id]
            return true
        end
    end
    return false
end

# A function to fill a prealocated matrix. The matrix is automatically extended if necessary.
function fill_preallocated_status_matrix!(current_size, stepsize, matrix, step, value1, value2, value3, value4, value5, value6, value7, value8, value9, value10, value11, value12)
    if step > current_size
        matrix_extension = Matrix{Int64}(undef, 12, stepsize)
        matrix = hcat(matrix, matrix_extension)
        current_size += stepsize
        matrix[:, step] .= value1, value2, value3, value4, value5, value6, value7, value8, value9, value10, value11, value12
    else
        matrix[:, step] .= value1, value2, value3, value4, value5, value6, value7, value8, value9, value10, value11, value12
    end

    # Return changed values
    return matrix, current_size
end

