# Collection of all functions related to periodic boundary conditions

#
# General head section of any file within this module 
#

#
# Include statements
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

# Function to replicate a given lattice based on translations in x,y,z directions
function replicate_lattice_with_translation(lattice, transx, transy)
    
    # Replicate
    fulllattice = copy(lattice)
    fulllattice[:,1] *= transx + 1
    fulllattice[:,2] *= transy + 1
    
    # Get inverse lattice
    invfulllattice = inv(fulllattice)

    # Return results
    return fulllattice, invfulllattice

end

# Function to replicate a set of grid points based on translations in x,y,z directions
function replicate_gridpoints_with_translation(lattice, gridpoints, transx, transy)
    
    # Define empty matrix
    gridmatrix = Array{Float64}(undef,0)

    # Loop over translation
    for i in 0:transx
        for j in 0:transy
            
            # Store in temp matrix
            tmp_matrix = copy(gridpoints)
            for point in axes(tmp_matrix, 2)
                tmp_matrix[:,point] += i * lattice[:,1] + j * lattice[:,2]
            end

            #println(tmp_matrix)

            # The first matrix is just stored, others are added to the matrix
            if i == 0 && j == 0
                gridmatrix = tmp_matrix
            else
                gridmatrix = hcat(gridmatrix, tmp_matrix)
            end
        end
    end 

    # Total number of gridpoints
    NGridpoints = size(gridmatrix, 2)

    # Generate mapping (gridpoint on unit cell gridpoint)
    # Format:
    # First row: unique grid point 
    # Second row: trans x value
    # Third row: trans y value
    # Meaning of column: i-th grid point
    NGridpoints_unitcell = size(gridpoints,2)
    NTranslations = (transx + 1) * (transy + 1)
    GridpointsMapping = Array{Int64}(undef, 3, NGridpoints)
    unitcell_point_id = [n for n = 1:NGridpoints_unitcell]

    GridpointsMapping[1,:] = repeat(unitcell_point_id, outer = NTranslations)
    iteration = 0
    for i in 0:transx
        for j in 0:transy

            iteration += 1
            range_start = (iteration - 1) * NGridpoints_unitcell + 1
            range_end = iteration * NGridpoints_unitcell
            GridpointsMapping[2,range_start:range_end] .= i
            GridpointsMapping[3,range_start:range_end] .= j

        end
    end

    # Return result
    return gridmatrix, NGridpoints, GridpointsMapping

end

"""

    apply_pbc_to_coordinates(coordinates, fulllattice, invfulllattice)

Function to move a list of points (given as vector or matrix) to the unit cell by applying pbc.
Importantly: The ranged used here is [-0.5, 0.5] for the fractional coordinates. Do not use this function for displacement.
Should only be used for a distance in case the unit cell is orthogonal. Better use [`apply_minimum_image_convention`](@ref) instead.

# Input
- `coordinates` in cartesian coordinates

# Return values
Returns same structure as in `coordinates` in cartesian coordinates.
"""
function apply_pbc_to_coordinates!(coordinates, fulllattice, invfulllattice)

    # Convert to fractional coordinates
    fractional_coordinates = invfulllattice * coordinates

    # Apply pbc
    for element in eachindex(fractional_coordinates)
        while fractional_coordinates[element] < -0.5
            fractional_coordinates[element] += 1.0
        end
        while fractional_coordinates[element] > 0.5
            fractional_coordinates[element] -= 1.0
        end
    end

    # Convert back to cartesian coordinates
    coordinates = fulllattice * fractional_coordinates

    # Return
    return coordinates

end

"""

    apply_pbc_to_fractional_coordinates(coordinates, fulllattice)

Function to move a list of points (given as vector or matrix) to the unit cell by applying pbc.
Importantly: The ranged used here is [-0.5, 0.5] for the fractional coordinates. Do not use this function for displacement.
Should only be used for a distance in case the unit cell is orthogonal. Better use [`apply_minimum_image_convention`](@ref) instead.

# Input
- `coordinates` in fractional coordinates

# Return values
Returns same structure as in `coordinates` in cartesian coordinates.
"""
function apply_pbc_to_fractional_coordinates!(coordinates, fulllattice)

    # Apply pbc
    for element in eachindex(coordinates)
        while coordinates[element] < -0.5
            coordinates[element] += 1.0
        end
        while coordinates[element] > 0.5
            coordinates[element] -= 1.0
        end
    end

    # Convert back tp cartesian coordinates
    cartesian_coordinates = fulllattice * coordinates

    # Return
    return cartesian_coordinates

end

"""

    closest_image_within_plane(vector, fulllattice)

Function to find the closest periodic image of a vector which is already wrapped into the unit cell.

# Input
- `vector` in cartesian coordinates

# Return values
Returns the shortest vector in cartesian coordinates.
"""
function closest_image_within_plane(vector, fulllattice)

    # For an orthogonal lattice the wrapped vector is always the closest periodic image
    # Uncomment in case the usage of minimum image convention is getting too slow otherwise run always the exact test
    #if abs(dot(fulllattice[:,1], fulllattice[:,2])) < 1e-10
    #    return vector
    #end

    # Check the neighbouring images within the surface plane
    closest_vector = vector
    closest_distance = norm(vector)
    for shift_x in -1:1
        for shift_y in -1:1
            candidate_vector = vector + shift_x * fulllattice[:,1] + shift_y * fulllattice[:,2]
            candidate_distance = norm(candidate_vector)
            if candidate_distance < closest_distance
                closest_vector = candidate_vector
                closest_distance = candidate_distance
            end
        end
    end

    # Return
    return closest_vector

end

"""

    minimum_image_vector(vector, fulllattice, invfulllattice)

Function to reduce a single distance vector to its closest periodic image.
The ranged used here is [-0.5, 0.5] for the fractional coordinates.

# Input
- `vector` in cartesian coordinates

# Return values
Returns the shortest vector in cartesian coordinates.
"""
function minimum_image_vector(vector, fulllattice, invfulllattice)

    # Wrap the vector into the unit cell
    fractional_vector = invfulllattice * vector
    fractional_vector .-= round.(fractional_vector)

    # Return the closest periodic image
    return closest_image_within_plane(fulllattice * fractional_vector, fulllattice)

end

"""

    minimum_image_vector(vector, fulllattice, invfulllattice)

Function to reduce a single distance vector to its closest periodic image.
The ranged used here is [-0.5, 0.5] for the fractional coordinates.

# Input
- `vector` in fractional coordinates

# Return values
Returns the shortest vector in cartesian coordinates.
"""
function minimum_image_from_fractional_coordinates(coordinates, fulllattice)

    # Wrap the fractional coordinates into the unit cell
    fractional_vector = coordinates .- round.(coordinates)

    # Return the closest periodic image
    return closest_image_within_plane(fulllattice * fractional_vector, fulllattice)

end

"""

    apply_minimum_image_convention(coordinates, fulllattice, invfulllattice)

Function to reduce distance vectors to their closest periodic images.
Accepts a single distance vector or a matrix storing one distance vector per column. 
Importantly: The ranged used here is [-0.5, 0.5] for the fractional coordinates. Do not use this function for displacement.

# Input
- `coordinates` in cartesian coordinates

# Return values
Returns same structure as in `coordinates` in cartesian coordinates.
"""
function apply_minimum_image_convention(coordinates, fulllattice, invfulllattice)

    # Reduce a single distance vector
    if ndims(coordinates) == 1
        return minimum_image_vector(coordinates, fulllattice, invfulllattice)
    end

    # Reduce every column of a matrix
    reduced_coordinates = Matrix{Float64}(undef, size(coordinates, 1), size(coordinates, 2))
    for column_id in axes(coordinates, 2)
        reduced_coordinates[:, column_id] = minimum_image_vector(coordinates[:, column_id], fulllattice, invfulllattice)
    end

    # Return
    return reduced_coordinates

end

# A function to map a unitcell gridpoint with its x & y translation to a final gridpoint
# This function automatically applies PBC to the translations
function map_translation_to_gridpoint(unique_point, Nunique_points, transx, Ncellx, transy, Ncelly)

    # Apply PBC to the translation values
    while transx < 0
        transx += Ncellx
    end
    while transx > (Ncellx - 1)
        transx -= Ncellx
    end
    while transy < 0
        transy += Ncelly
    end
    while transy > (Ncelly - 1)
        transy -= Ncelly
    end
    
    # Return the number of the gridpoint
    return unique_point + transy * Nunique_points + transx * Ncelly * Nunique_points

end

# A function to derive the translation vectors needed to plot the periodic images of an adsorbate
# Adsorbates are only shifted in case they are located within the boundary cells of the lattice
# The returned vector is empty for all adsorbates outside of the boundary cells
function boundary_translation_vectors(transx, transy, boundary_cells, lattice)

    # Create the vector storing the translations
    translations = Vector{Vector{Float64}}(undef, 0)

    # Increase the translations by one to adapt the scale to "1 to Ncells"
    transx += 1
    transy += 1

    # Derive the direction of the shift along the first lattice vector
    shift_x = 0
    if transx ≤ boundary_cells
        shift_x = 1
    elseif (lattice.Ncellx - boundary_cells) < transx
        shift_x = -1
    end

    # Derive the direction of the shift along the second lattice vector
    shift_y = 0
    if transy ≤ boundary_cells
        shift_y = 1
    elseif (lattice.Ncelly - boundary_cells) < transy
        shift_y = -1
    end

    # Collect all needed translations
    # An adsorbate within a corner of the lattice needs three periodic images
    if shift_x != 0
        push!(translations, shift_x * lattice.transcellvectors[:,1])
    end
    if shift_y != 0
        push!(translations, shift_y * lattice.transcellvectors[:,2])
    end
    if shift_x != 0 && shift_y != 0
        push!(translations, shift_x * lattice.transcellvectors[:,1] + shift_y * lattice.transcellvectors[:,2])
    end

    # Return the translations
    return translations

end