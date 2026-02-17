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
    fulllattice[1,:] *= transx + 1
    fulllattice[2,:] *= transy + 1
    
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
                tmp_matrix[:,point] += i * lattice[1,:] + j * lattice[2,:]
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

# Function to move a point to the unit cell or reduce a distance by applying pbc
function apply_pbc_to_coordinates(coordinates, fulllattice, invfulllattice)

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

# Function to move a point to the unit cell or reduce a distance by applying pbc
# The initial point is already provided in fractional coordinates
# Returns cartesian coordinates
function apply_pbc_to_fractional_coordinates(coordinates, fulllattice)

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

function apply_pbc_to_translation_value(value, boundary)
 
    # Apply PBC to the translation values
    # Boundary is only equal to the actual boundary in case the "0" is not a valid element. Otherwise boundary should specify the number of intersections.
    while value < 0
        value += boundary
    end
    while value > (boundary - 1)
        value -= boundary
    end

    return value

end