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

"""

    is_reduced_lattice_basis(fulllattice; tolerance = 1e-3)

Function to check whether the two in-plane vectors of the full lattice form a reduced (Lagrange-Gauss) basis,
i.e. |a1·a2| ≤ min(|a1|², |a2|²) / 2. Only for a reduced basis is the closest periodic image of a wrapped
distance vector guaranteed to be among the neighbouring images checked in [`closest_image_within_plane`](@ref).
Orthogonal cells and hexagonal cells with an angle of 60° or 120° and equal translations are reduced.
The hexagonal case lies exactly on the boundary, so a relative `tolerance` absorbs the rounding of the input vectors.

# Input
- `fulllattice` full lattice (lattice vectors as columns, including the translations)

# Return values
Returns true in case the basis is reduced.
"""
function is_reduced_lattice_basis(fulllattice; tolerance = 1e-3)

    # In-plane lattice vectors
    vector_1 = @view fulllattice[:,1]
    vector_2 = @view fulllattice[:,2]

    # Reduction criterion
    shortest_squared = min(dot(vector_1, vector_1), dot(vector_2, vector_2))
    return abs(dot(vector_1, vector_2)) <= 0.5 * shortest_squared * (1.0 + tolerance)

end

"""

    reduce_lattice_basis(fulllattice, invfulllattice)

Function to generate a reduced (Lagrange-Gauss) basis of the full lattice for the minimum image convention.
The reduced basis spans exactly the same periodic images as the full lattice, so any choice of transx and transy
is possible, e.g. a hexagonal unit cell with transx ≠ transy. The full lattice itself is not changed, as the
translations and the mapping of the grid points are defined with respect to it.
A basis that is already reduced is returned unchanged, so the results for such cells are not affected.

# Input
- `fulllattice` full lattice (lattice vectors as columns, including the translations)
- `invfulllattice` inverse of the full lattice

# Return values
- `reducedlattice` reduced basis (lattice vectors as columns)
- `invreducedlattice` inverse of the reduced basis
- `reductionmatrix` integer matrix converting fractional coordinates of the full lattice into fractional coordinates of the reduced basis
"""
function reduce_lattice_basis(fulllattice, invfulllattice)

    # Nothing to do for a reduced basis
    dimension = size(fulllattice, 1)
    if is_reduced_lattice_basis(fulllattice)
        return copy(fulllattice), copy(invfulllattice), Matrix{Int64}(I, dimension, dimension)
    end

    # Lagrange-Gauss reduction of the two in-plane vectors
    # The integer coefficients (columns of coefficients) express the reduced vectors in the original ones
    vector_1 = fulllattice[:,1]
    vector_2 = fulllattice[:,2]
    coefficients_1 = [1, 0]
    coefficients_2 = [0, 1]
    if dot(vector_1, vector_1) > dot(vector_2, vector_2)
        vector_1, vector_2 = vector_2, vector_1
        coefficients_1, coefficients_2 = coefficients_2, coefficients_1
    end
    for iteration in 1:1000
        factor = round(Int64, dot(vector_1, vector_2) / dot(vector_1, vector_1))
        vector_2 = vector_2 - factor * vector_1
        coefficients_2 = coefficients_2 - factor * coefficients_1
        if dot(vector_2, vector_2) < dot(vector_1, vector_1)
            vector_1, vector_2 = vector_2, vector_1
            coefficients_1, coefficients_2 = coefficients_2, coefficients_1
        else
            break
        end
        iteration == 1000 && error("Lattice basis reduction did not converge")
    end

    # Reduced basis (out-of-plane vectors are kept)
    reducedlattice = copy(fulllattice)
    reducedlattice[:,1] = vector_1
    reducedlattice[:,2] = vector_2
    if !is_reduced_lattice_basis(reducedlattice)
        error("Lattice basis reduction failed")
    end

    # Fractional coordinates: full lattice = reduced basis * inverse coefficient matrix
    # The coefficient matrix is unimodular, so its inverse is an integer matrix as well
    determinant = coefficients_1[1] * coefficients_2[2] - coefficients_2[1] * coefficients_1[2]
    reductionmatrix = Matrix{Int64}(I, dimension, dimension)
    reductionmatrix[1,1] = determinant * coefficients_2[2]
    reductionmatrix[1,2] = -determinant * coefficients_2[1]
    reductionmatrix[2,1] = -determinant * coefficients_1[2]
    reductionmatrix[2,2] = determinant * coefficients_1[1]

    # Return results
    return reducedlattice, inv(reducedlattice), reductionmatrix

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

    closest_image_within_plane(vector, fulllattice)

Function to find the closest periodic image of a vector which is already wrapped into the unit cell.

# Input
- `vector` in cartesian coordinates

# Return values
Returns the shortest vector in cartesian coordinates.
"""
function closest_image_within_plane(vector, fulllattice)

    # Check the neighbouring images within the surface plane
    # Every candidate is built in one reused buffer; only the shift of the closest image is stored
    candidate_vector = similar(vector, Float64)
    closest_distance = norm(vector)
    closest_shift_x = 0
    closest_shift_y = 0
    for shift_x in -1:1
        for shift_y in -1:1
            for i in eachindex(candidate_vector)
                candidate_vector[i] = vector[i] + shift_x * fulllattice[i,1] + shift_y * fulllattice[i,2]
            end
            candidate_distance = norm(candidate_vector)
            if candidate_distance < closest_distance
                closest_distance = candidate_distance
                closest_shift_x = shift_x
                closest_shift_y = shift_y
            end
        end
    end

    # The wrapped vector itself is the closest image
    if closest_shift_x == 0 && closest_shift_y == 0
        return vector
    end

    # Rebuild the closest image
    for i in eachindex(candidate_vector)
        candidate_vector[i] = vector[i] + closest_shift_x * fulllattice[i,1] + closest_shift_y * fulllattice[i,2]
    end

    # Return
    return candidate_vector

end

"""

    minimum_image_vector(vector, fulllattice, invfulllattice)

Function to reduce a single distance vector to its closest periodic image.
The range used here is [-0.5, 0.5] for the fractional coordinates.

# Input
- `vector` in cartesian coordinates

# Return values
Returns the shortest vector in cartesian coordinates.
"""
function minimum_image_vector(vector, fulllattice, invfulllattice)

    # Wrap the vector into the unit cell
    fractional_vector = invfulllattice * vector
    for i in eachindex(fractional_vector)
        fractional_vector[i] -= round(fractional_vector[i])
    end

    # Return the closest periodic image
    return closest_image_within_plane(fulllattice * fractional_vector, fulllattice)

end

"""

    minimum_image_from_fractional_coordinates(coordinates, reducedlattice, reductionmatrix)

Function to reduce a single distance vector to its closest periodic image.
The range used here is [-0.5, 0.5] for the fractional coordinates of the reduced basis.

# Input
- `coordinates` distance vector in fractional coordinates of the full lattice
- `reducedlattice` reduced basis of the full lattice, see [`reduce_lattice_basis`](@ref)
- `reductionmatrix` integer matrix converting fractional coordinates of the full lattice into the reduced basis

# Return values
Returns the shortest vector in cartesian coordinates.
"""
function minimum_image_from_fractional_coordinates(coordinates, reducedlattice, reductionmatrix)

    # Wrap the fractional coordinates (in the reduced basis) into the unit cell
    if isone(reductionmatrix)
        fractional_vector = coordinates .- round.(coordinates)
    else
        fractional_vector = reductionmatrix * coordinates
        for i in eachindex(fractional_vector)
            fractional_vector[i] -= round(fractional_vector[i])
        end
    end

    # Return the closest periodic image
    return closest_image_within_plane(reducedlattice * fractional_vector, reducedlattice)

end

"""

    apply_minimum_image_convention(coordinates, fulllattice, invfulllattice)

Function to reduce distance vectors to their closest periodic images.
Accepts a single distance vector or a matrix storing one distance vector per column. 
Importantly: The range used here is [-0.5, 0.5] for the fractional coordinates. Do not use this function for displacement.

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
        reduced_coordinates[:, column_id] = minimum_image_vector(@view(coordinates[:, column_id]), fulllattice, invfulllattice)
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