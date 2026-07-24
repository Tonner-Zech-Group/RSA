# Collection of all functions related to RSA initialization and preparation

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

# 
# Definition of structs
# Use @kwdef to define default values
#

@kwdef mutable struct rate_constants_struct
    ads::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)         # [Molecule, Gridtype]: rate constant
    rot::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)         # [Molecule, Gridtype]: rate constant
    dif::Array{Float64, 3} = Array{Float64}(undef, 0, 0, 0)        # [Molecule, Gridtype, Gridtype]: rate constant
    dif_rad::Array{Float64, 3} = Array{Float64}(undef, 0, 0, 0)    # [Molecule, Gridtype, Gridtype]: diffusion radius
    con::Array{Float64, 3} = Array{Float64}(undef, 0, 0, 0)        # [Molecule, Gridtype, Molecule]: rate constant
end

# This struct is generated to hold all information for a single gridpoint in a kMC run
# Contained information are:
# Number of possible events, Total rate constant of all events, possible rotations for an adsorbing adsorbate
mutable struct gridpoint_struct
    Nevents::Int64 # Total number of possible events
    Nads::Int64
    Ndif::Int64
    Nrot::Int64
    Ncon::Int64
    Trateconst::Float64 # Total rate constant
    Trateconst_ads::Float64 # Total rate constant with adsorption only
    Trate_molec::Vector{Float64} # Total rate constant per molecule
    Cumulative_rate_molecules::Vector{Float64}
    Cumulative_rate_molecules_ads::Vector{Float64} # Cumulative rate constant per molecule with ads only
    Trate_event::Vector{Vector{Float64}} # Total rate constant per event for every molecule; Obtained by Nfree_rotation and rate constants per event type
                                         # 1: Ads, 2: Rot, 3: Dif, 4: Con
    Cumulative_rate_events::Vector{Vector{Float64}}
    Trate_diffusion::Vector{Float64} # Total rate constant per diffusion per gridtype (target)
    Cumulative_rate_diffusions::Vector{Float64}
    Trate_conformer::Vector{Float64} # Total rate constant per conformer change per molecule (target)
    Cumulative_rate_conformers::Vector{Float64}

    target_diffusion::Vector{Vector{Vector{Int64}}} # Final sites, which can be reached by a diffusion event: [for Molecule][Gridtype][target points]
    bool_diffusion::Vector{Vector{Vector{Bool}}} # store whether diffusion target is free (true/false)
    reduced_diffusion::Vector{Vector{Vector{Int64}}}
    Ntarget_diffusion::Vector{Vector{Int64}} # number of free target sites per molecule and per gridtype (all information stored to recduce work in case conformer is changed)
    origin_diffusion::Vector{Vector{Vector{Int64}}} # Initial sites, which can be reach this point by a diffusion event: [for Molecule][Gridtype][initial points]

    info_event::Vector{Bool} # can a molecule in principle stick to this point       
    bool_occupied::Bool # Is a molecule adsorbed on this gridpoint
    adsorbed_molecule::Int64 # Which molecule is adsorbed on this point
    adsorbed_rotation::Int64 # In which rotation state is the molecule adsorbed 
    bool_rotation::Vector{Vector{Bool}} # Vector storing the possible rotation states per molecule
    reduced_rotation::Vector{Vector{Int64}} # Vector only storing the position of possible rotation states
    Nfree_rotation::Vector{Int64} # Number of free rotation states per molecule
    Nblocked_rotation::Vector{Vector{Int64}} # Vector counting how often a rotation was blocked by an adjacent adsorbate or by the adsorbate
end

#
# Function definitions
#

# A function generating a "gridpoint_struct" struct for every gridpoint and storing all values in a vector
function grid_initialization(Ngrids, grids, Nmolecules, molecules, lattice, rate_constants_info, neighbour_list)
    
    # Structure of the grid points vector is [gridtype][gridpoint]
    rsa_gridpoints = Vector{Vector{gridpoint_struct}}(undef, Ngrids)
    for grid_id in 1:Ngrids
        rsa_gridpoints[grid_id] = Vector{gridpoint_struct}(undef, grids[grid_id].Npoints)
    end

    # Loop over gridtypes
    for grid_id in 1:Ngrids

        # Generate tmp matrix
        pointsvec = Vector{gridpoint_struct}(undef, grids[grid_id].Npoints)

        # General information
        adsorbed_molecule = 0
        adsorbed_rotation = 0
        bool_occupied = false
        
        # General information per molecule
        info_event = Vector{Bool}(undef, Nmolecules)
        info_event .= false
        for molecule_id in 1:Nmolecules
            if ! isempty(findall(x -> x == grid_id, molecules[molecule_id].grids))
                info_event[molecule_id] = true
            end
        end

        # Generate the Nfree_rotation
        Nfree_rotation = Vector{Int64}(undef, Nmolecules)
        Nfree_rotation .= 0
        for molecule_id in 1:Nmolecules
            if info_event[molecule_id] == true
                Nfree_rotation[molecule_id] = molecules[molecule_id].Nrotations
            end
        end

        # Generate the bool_rotation
        bool_rotation = Vector{Vector{Bool}}(undef, Nmolecules)
        for molecule_id in 1:Nmolecules
            if info_event[molecule_id] == true
                bool_rotation[molecule_id] = [true for i in 1:molecules[molecule_id].Nrotations]
            else
                bool_rotation[molecule_id] = []
            end
        end

        # Generate reduced_rotation
        reduced_rotation = Vector{Vector{Int64}}(undef, Nmolecules)
        for molecule_id in 1:Nmolecules
            reduced_rotation[molecule_id] = findall(x->x==true, bool_rotation[molecule_id])
        end

        # Generate the Nblocked_rotation
        Nblocked_rotation = Vector{Vector{Int64}}(undef, Nmolecules)
        for molecule_id in 1:Nmolecules
            if info_event[molecule_id] == true
                Nblocked_rotation[molecule_id] = [0 for i in 1:molecules[molecule_id].Nrotations]
            else
                Nblocked_rotation[molecule_id] = []
            end
        end

        # Generate Trate_conformer and Cumulative_rate_conformers
        Trate_conformer = Vector{Float64}(undef, Nmolecules)
        Trate_conformer .= 0.0
        Cumulative_rate_conformers = Vector{Float64}(undef, Nmolecules)
        Cumulative_rate_conformers .= 0.0

        # Generate Trate_diffusion and Cumulative_rate_diffusions
        Trate_diffusion = Vector{Float64}(undef, Ngrids)
        Trate_diffusion .= 0.0
        Cumulative_rate_diffusions = Vector{Float64}(undef, Ngrids)
        Cumulative_rate_diffusions .= 0.0

        # Generate Trate_event
        # Only use positve rate constants
        Trate_event = Vector{Vector{Float64}}(undef, Nmolecules)
        for molecule_id in 1:Nmolecules
            Trate_event[molecule_id] = zeros(4)
            if rate_constants_info.ads[molecule_id, grid_id] < 0
                Trate_event[molecule_id][1] = 0.0
            else
                Trate_event[molecule_id][1] = Nfree_rotation[molecule_id] * rate_constants_info.ads[molecule_id, grid_id]
            end
        end

        # Generate Cumulative_rate_events
        Cumulative_rate_events = Vector{Vector{Float64}}(undef, Nmolecules)
        for molecule_id in 1:Nmolecules
            Cumulative_rate_events[molecule_id] = zeros(4)
            cum_sum = 0
            for event_id in 1:4
                cum_sum += Trate_event[molecule_id][event_id]
                Cumulative_rate_events[molecule_id][event_id] = cum_sum
            end
        end

        # Generate Trate_molec 
        Trate_molec = Vector{Float64}(undef, Nmolecules)
        Trate_molec .= 0.0
        for molecule_id in 1:Nmolecules
            if rate_constants_info.ads[molecule_id, grid_id] < 0
                Trate_molec[molecule_id] = 0.0
            else
                Trate_molec[molecule_id] = Trate_event[molecule_id][1]
            end
        end

        # Generate Cumulative_rate_molecules
        Cumulative_rate_molecules = Vector{Float64}(undef, Nmolecules)
        Cumulative_rate_molecules .= 0
        cum_sum = 0
        for molecule_id in 1:Nmolecules
            cum_sum += Trate_molec[molecule_id]
            Cumulative_rate_molecules[molecule_id] = cum_sum
        end

        # Generate Cumulative_rate_molecules_ads
        Cumulative_rate_molecules_ads = Vector{Float64}(undef, Nmolecules)
        Cumulative_rate_molecules_ads .= 0
        cum_sum = 0
        for molecule_id in 1:Nmolecules
            cum_sum += Trate_event[molecule_id][1]
            Cumulative_rate_molecules_ads[molecule_id] = cum_sum
        end

        # Generate Trateconst
        Trateconst = Cumulative_rate_molecules[Nmolecules]

        # Generate Trateconst_ads
        Trateconst_ads = Cumulative_rate_molecules_ads[Nmolecules]

        # Generate Nevents
        Nads = sum(Nfree_rotation)
        Nevents = Nads
        Ndif = 0
        Nrot = 0
        Ncon = 0
        
        # Initialize all points
        for points_id in 1:grids[grid_id].Npoints

            # Get the unique point and mapping of this grid point
            unique_point, mapping_x, mapping_y = grids[grid_id].mapping[:,points_id]

            # Generate the list of possible diffusion targets (these are relative to the unique point)
            # target_diffusion = deepcopy(neighbour_list[grid_id][unique_point]) that creates a lot of overhead especially for memory allocation
            target_diffusion = [[copy(end_grid_list) for end_grid_list in molecule_list] for molecule_list in neighbour_list[grid_id][unique_point]]

            # Generate Ntarget_diffusion
            Ntarget_diffusion = Vector{Vector{Int64}}(undef, Nmolecules)
            for molecule_id in 1:Nmolecules
                Ntarget_diffusion[molecule_id] = Vector{Int64}(undef, Ngrids)
                for grid_2_id in 1:Ngrids
                    Ntarget_diffusion[molecule_id][grid_2_id] = size(target_diffusion[molecule_id][grid_2_id], 1)
                end
            end

            # Generate bool_diffusion
            bool_diffusion = Vector{Vector{Vector{Bool}}}(undef, Nmolecules)
            for molecule_id in 1:Nmolecules
                bool_diffusion[molecule_id] = Vector{Vector{Bool}}(undef, Ngrids)
                for grid_2_id in 1:Ngrids
                    bool_diffusion[molecule_id][grid_2_id] = [false for i in 1:Ntarget_diffusion[molecule_id][grid_2_id]]
                end
            end

            # Generate reduced_diffusion
            reduced_diffusion = Vector{Vector{Vector{Int64}}}(undef, Nmolecules)
            for molecule_id in 1:Nmolecules
                reduced_diffusion[molecule_id] = Vector{Vector{Int64}}(undef, Ngrids)
                for grid_2_id in 1:Ngrids
                    reduced_diffusion[molecule_id][grid_2_id] = []
                end
            end

            # Change from unique points to final points
            for molecule_id in 1:Nmolecules
                for grid_2_id in 1:Ngrids
                    for points_2_id in 1:Ntarget_diffusion[molecule_id][grid_2_id]

                        # Get the information of the target point
                        target_unique_point, target_mapping_x, target_mapping_y = grids[grid_2_id].mapping[:,target_diffusion[molecule_id][grid_2_id][points_2_id]]

                        # Generate the information of the final point
                        final_unique_point = target_unique_point
                        final_mapping_x = mapping_x + target_mapping_x
                        final_mapping_y = mapping_y + target_mapping_y
                        final_point = map_translation_to_gridpoint(final_unique_point, grids[grid_2_id].Nuniquepoints, final_mapping_x, lattice.Ncellx, final_mapping_y, lattice.Ncelly)

                        # Replace the old point
                        target_diffusion[molecule_id][grid_2_id][points_2_id] = final_point

                    end
                end
            end
            
            # Generate the placeholder for origin_diffusion
            origin_diffusion = Vector{Vector{Vector{Int64}}}(undef, Nmolecules)
            for molecule_id in 1:Nmolecules
                origin_diffusion[molecule_id] = Vector{Vector{Int64}}(undef, Ngrids)
                for grid_2_id in 1:Ngrids
                    origin_diffusion[molecule_id][grid_2_id] = []
                end
            end
            
            # Get the final struct with all information previously collected
            # To replace the deepcopy (again a lot of overhead):
            # target_diffusion, bool_diffusion, reduced_diffusion, Ntarget_diffusion, origin_diffusion are created for every point and therfore do not need to be copied
            # All integer, float, and boolean values are copied by default: Nevents, Nads, Ndif, Nrot, Ncon, Trateconst, Trateconst_ads, bool_occupied, adsorbed_molecule, adsorbed_rotation
            # All vectors get a single "copy": Trate_molec, Cumulative_rate_molecules, Cumulative_rate_molecules_ads, info_event, Nfree_rotation
            # Vectors of vectors need a copy in a for loop: Trate_event, Cumulative_rate_events, bool_rotation, reduced_rotation, Nblocked_rotation 
            pointsvec[points_id] = gridpoint_struct(Nevents, Nads, Ndif, Nrot, Ncon, Trateconst, Trateconst_ads, 
                copy(Trate_molec), copy(Cumulative_rate_molecules), copy(Cumulative_rate_molecules_ads), 
                [copy(ele) for ele in Trate_event], [copy(ele) for ele in Cumulative_rate_events], 
                copy(Trate_diffusion), copy(Cumulative_rate_diffusions), copy(Trate_conformer), copy(Cumulative_rate_conformers), 
                target_diffusion, bool_diffusion, reduced_diffusion, Ntarget_diffusion, origin_diffusion, 
                copy(info_event), 
                bool_occupied, adsorbed_molecule, adsorbed_rotation, 
                [copy(ele) for ele in bool_rotation], [copy(ele) for ele in reduced_rotation], 
                copy(Nfree_rotation), 
                [copy(ele) for ele in Nblocked_rotation])

        end

        # Store in the final vector
        rsa_gridpoints[grid_id] = pointsvec

    end

    # As now all points are defined, we can use these information to create the origin_diffusion vectors
    for grid_id in 1:Ngrids
        for points_id in 1:grids[grid_id].Npoints
            
            # Use every target point in the target_diffusion vector and add the initial point to the corresponding origin_diffusion 
            for molecule_id in 1:Nmolecules
                for grid_2_id in 1:Ngrids
                    for ele_id in eachindex(rsa_gridpoints[grid_id][points_id].target_diffusion[molecule_id][grid_2_id])
                        target_point_id = rsa_gridpoints[grid_id][points_id].target_diffusion[molecule_id][grid_2_id][ele_id]
                        push!(rsa_gridpoints[grid_2_id][target_point_id].origin_diffusion[molecule_id][grid_id], points_id)
                    end
                end
            end

        end
    end

    # Return result
    return rsa_gridpoints

end

# Keep for debugging purposes, but not used in the current version of the code
# Function to calculate the overlap of two molecules in 2D
function check_for_overlap(molecule_coords_A, molecule_coords_B, molecule_elements_A, molecule_elements_B, fulllattice, invfulllattice)

    # Default
    overlap = false

    # Check
    for a in axes(molecule_coords_A,2)
        for b in axes(molecule_coords_B,2)
            
            # VdW sum
            sum_vdw = (atomic_information[molecule_elements_A[a], 3] + atomic_information[molecule_elements_B[b], 3])
            #println("VdW Radius: " * string(sum_vdw))

            # Calculate distance (Difference --> convert to fractional --> check pbc --> convert back --> distance)
            difference = molecule_coords_A[:,a] - molecule_coords_B[:,b]
            
            # Apply pbc
            difference = apply_pbc_to_coordinates(difference, fulllattice, invfulllattice)

            # Calculate distance
            distance = norm(difference[1:2])
            #println("Distance: " * string(distance))
            
            if distance < sum_vdw
                overlap = true
                break
            end

        end

        if overlap == true
            break
        end

    end

    # Return
    return overlap

end

# Current version
# Function to check whether atoms within two structures are overlapping
# This function is using distance and difference matrices for unitcell gridpoints, translation vectors and rotations
# To judge the overlap the distance is projected to the xy-plane
function check_for_atom_overlap(points_difference_vector, points_difference, molecules, molecule_id, rotation_id, Tmolecule_id, Trotation_id, lattice, rotation_difference_matrices, overlap_2d)
    
    #unitcell_point_id, unitcell_rotation_id, gridpoint_id, gridpoint_rotation_id, molecule_elements, 
    #GridpointsMapping, unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, fulllattice)

    # Input:
    # points_difference_vector  - the distance vector between the gridpoints
    # molecules                 - the coordinates
    # molecule_id, rotation_id  - IDs of the first molecule
    # Tmolecule_id, Trotation_id- IDs of the second (target) molecule
    # lattice                   - lattice vectors for transformation
    # 2d_overlap                - bool to request overlap only in 2 dimensions
    
    # Default
    overlap = false
    breakcondition_length = false

    # Get the rotation distance matrix once
    rotation_difference_matrix = rotation_difference_matrices[molecule_id,Tmolecule_id][rotation_id,Trotation_id]
    #length_sums_matrix = length_sums[molecule_id,Tmolecule_id] --> using this seems to slightly slow down computations
    #elements_1 = molecules[molecule_id].elements_sorted --> using this seems to slightly slow down computations
    #elements_2 = molecules[Tmolecule_id].elements_sorted --> using this seems to slightly slow down computations

    # Loop over all atoms in molecule 1
    for atom_1_id in axes(molecules[molecule_id].elements_sorted, 1)

        length_1 = molecules[molecule_id].radii_sorted[atom_1_id]
        vdw_1 = atomic_information[molecules[molecule_id].elements_sorted[atom_1_id], 3]
        #vdw_1 = atomic_information[elements_1[atom_1_id], 3]

        # Loop over all atoms in molecule 2
        for atom_2_id in axes(molecules[Tmolecule_id].elements_sorted, 1)

            length_2 = molecules[Tmolecule_id].radii_sorted[atom_2_id]
            vdw_2 = atomic_information[molecules[Tmolecule_id].elements_sorted[atom_2_id], 3]
            #vdw_2 = atomic_information[elements_2[atom_2_id], 3]

            # Prescreen using length vector
            if points_difference >  length_1 + length_2 # length_sums_matrix[atom_1_id,atom_2_id]
                
                # Set a flag in case even the molecule with the largest radius is not sufficient
                if atom_2_id == 1
                    breakcondition_length = true
                    break
                end
                
                # Get the next atom of molecule 1
                break
            end
 
            # Calculate the VdW sum
            sum_vdw = vdw_1 + vdw_2

            # Calculate the difference between the atoms (vector in fractional coordinates)
            fractional_distance_vector = points_difference_vector + rotation_difference_matrix[atom_1_id,atom_2_id]
            cartesian_distance_vector = apply_pbc_to_fractional_coordinates(fractional_distance_vector, lattice.transvectors)
            # Only consider the distance in plane for RSA
            if overlap_2d
                atoms_distance = norm(cartesian_distance_vector[1:2])
            else
                atoms_distance = norm(cartesian_distance_vector)
            end

            # Check for overlap
            if atoms_distance < sum_vdw
                overlap = true
                break
            end

        end

        # Check for overlap and length flag
        if overlap == true || breakcondition_length == true
            break
        end

    end

    # Return result
    return overlap

end

# CURRENTLY: Not used
# A function to calculate every sum of two atom radii
# Aiming at reducing the computation of the effected gridpoints for the eventlist
function calculate_length_sums(Nmolecules, molecules)
    
    # Create the matrix to store all elements
    length_sums = Matrix{Matrix{Float64}}(undef, Nmolecules, Nmolecules)

    for molecule_A_id in 1:Nmolecules
        for molecule_B_id in 1:Nmolecules

            # Generate tmp matrix
            tmp_length_sums = Matrix{Float64}(undef, molecules[molecule_A_id].Natoms, molecules[molecule_B_id].Natoms)

            # Calculate these values
            for atom_A_id in 1:molecules[molecule_A_id].Natoms
                for atom_B_id in 1:molecules[molecule_B_id].Natoms
                    length_1 = molecules[molecule_A_id].radii_sorted[atom_A_id]
                    length_2 = molecules[molecule_B_id].radii_sorted[atom_B_id]
                    tmp_length_sums[atom_A_id,atom_B_id] = length_1 + length_2
                end
            end

            # Store in the final matrix
            length_sums[molecule_A_id,molecule_B_id] = tmp_length_sums

        end
    end

    # Return the results
    return length_sums

end


# Function to calculate the difference vectors for all pairs of gridpoints within the unitcell
# PBC are not applied
# Final vectors are converted to fractional coordinates
function calculate_unit_cell_gridpoints_difference_vectors(Ngrids, grids, lattice)
    
    # Create the matrix to store all elements
    # Level 1: Grid A vs Grid B
    # Level 2: Gridpoint (A) vs Gridpoint (B)
    unit_cell_gridpoints_difference = Matrix{Matrix{Vector{Float64}}}(undef, Ngrids, Ngrids)

    # Loop over gridpoints
    for grid_B_id in 1:Ngrids
        for grid_A_id in grid_B_id:Ngrids
            # Generate tmp matrix
            gridpoints_difference_matrix = Matrix{Vector{Float64}}(undef, grids[grid_A_id].Nuniquepoints, grids[grid_B_id].Nuniquepoints)

            # Calculate the values
            if grid_A_id == grid_B_id
                # The final matrix is anti-symmetric - so only calculate half its values
                Threads.@threads for column_id in 1:grids[grid_A_id].Nuniquepoints
                    for row_id in column_id + 1:grids[grid_A_id].Nuniquepoints
                        gridpoints_difference_matrix[row_id, column_id] = lattice.inversevectors * (grids[grid_A_id].uniquepoints[:, column_id] - grids[grid_A_id].uniquepoints[:, row_id])
                        gridpoints_difference_matrix[column_id, row_id] = -gridpoints_difference_matrix[row_id, column_id]
                    end
                end

                # Set the diagonal to zero-vectors
                Threads.@threads for column_id in 1:grids[grid_A_id].Nuniquepoints
                    gridpoints_difference_matrix[column_id, column_id] = [0.0, 0.0, 0.0]
                end
            
            elseif grid_A_id != grid_B_id
                # The final matrix is not symmetric - so calculate all its values
                Threads.@threads for column_id in 1:grids[grid_B_id].Nuniquepoints
                    for row_id in 1:grids[grid_A_id].Nuniquepoints
                        gridpoints_difference_matrix[row_id, column_id] = lattice.inversevectors * (grids[grid_B_id].uniquepoints[:, column_id] - grids[grid_A_id].uniquepoints[:, row_id])
                    end
                end

            end

            # Store the result in the final matrix
            if grid_B_id == grid_A_id
                unit_cell_gridpoints_difference[grid_A_id, grid_A_id] = gridpoints_difference_matrix
            elseif grid_B_id != grid_A_id
                unit_cell_gridpoints_difference[grid_A_id, grid_B_id] = gridpoints_difference_matrix
                unit_cell_gridpoints_difference[grid_B_id, grid_A_id] = -permutedims(gridpoints_difference_matrix,[2,1])
            end

        end
    end

    # Return the result
    return unit_cell_gridpoints_difference

end

# Function to calculate the distance vectors for all pairs of translation vectors
# To obtain a "difference" between two gridpoints two values of this distance matrix must be substracted
# PBC are not applied
# Final vectors are converted to fractional coordinates
function calculate_translation_distance_vectors(lattice)
    
    # Default
    translation_distance_vectors = Matrix{Vector{Float64}}(undef, lattice.transx + 1, lattice.transy + 1)
    lattice_x = lattice.vectors[1,1:3]
    lattice_y = lattice.vectors[2,1:3]

    # The matrix is not symmetric or anti-symmetric
    Threads.@threads for column_id in 1:lattice.transy + 1
        for row_id in 1:lattice.transx + 1
            translation_distance_vectors[row_id, column_id] = lattice.inversevectors * ((row_id - 1) * lattice_x + (column_id - 1) * lattice_y)
        end
    end

    # Return results
    return translation_distance_vectors

end

# Function to calculate the difference vectors for all pairs of atoms within any pair of rotation values
# PBC are not applied
# Final vectors are converted to fractional coordinates
function calculate_rotation_difference_matrices(Nmolecules, molecules, lattice)
    
    # Create the initial matrix to store all results
    rotation_difference_matrices = Matrix{Matrix{Matrix{Vector{Float64}}}}(undef, Nmolecules, Nmolecules)

    # First loop over the molecules
    for molecule_A_id in 1:Nmolecules
        #for molecule_B_id in molecule_A_id:Nmolecules
        for molecule_B_id in molecule_A_id:Nmolecules

            # Generate tmp matrix
            tmp_molecule_difference_matrix = Matrix{Matrix{Vector{Float64}}}(undef, molecules[molecule_B_id].Nrotations, molecules[molecule_A_id].Nrotations)

            # Based on whether the molecule_ids are identical adapt the calculation
            if molecule_A_id == molecule_B_id

                # Matrix is anti-symmetric
                Threads.@threads for rotation_A_id in 1:molecules[molecule_A_id].Nrotations
                    # Get one coordinate matrix
                    structure_2 = molecules[molecule_A_id].coordinates_rotated[rotation_A_id]
                
                    for rotation_B_id in rotation_A_id:molecules[molecule_B_id].Nrotations
                        # Get second coordinate matrix
                        structure_1 = molecules[molecule_B_id].coordinates_rotated[rotation_B_id]

                        # Define blank matrix
                        tmp_matrix = Array{Vector{Float64}}(undef, molecules[molecule_B_id].Natoms, molecules[molecule_A_id].Natoms)
            
                        # Calculate its values
                        # This matrix is anti-symmetric in case the rotation_id and molecule_id match (two times the same structure)
                        if rotation_A_id == rotation_B_id
                            for atom_A_id in 1:molecules[molecule_A_id].Natoms
                                for atom_B_id in atom_A_id:molecules[molecule_B_id].Natoms
                    
                                    if atom_A_id == atom_B_id
                                        tmp_matrix[atom_B_id, atom_A_id] = lattice.inversevectors * (structure_2[:,atom_A_id] - structure_1[:,atom_B_id])
                                    else
                                        tmp_matrix[atom_B_id, atom_A_id] = lattice.inversevectors * (structure_2[:,atom_A_id] - structure_1[:,atom_B_id])
                                        tmp_matrix[atom_A_id,atom_B_id] = -tmp_matrix[atom_B_id, atom_A_id]
                                    end

                                end
                            end
                        else
                            for atom_A_id in 1:molecules[molecule_A_id].Natoms
                                for atom_B_id in 1:molecules[molecule_B_id].Natoms
                                    tmp_matrix[atom_B_id, atom_A_id] = lattice.inversevectors * (structure_2[:,atom_A_id] - structure_1[:,atom_B_id])
                                end
                            end
                        end

                        # Store it in the final matrix
                        if rotation_A_id == rotation_B_id
                            tmp_molecule_difference_matrix[rotation_B_id, rotation_A_id] = tmp_matrix
                        else
                            tmp_molecule_difference_matrix[rotation_B_id, rotation_A_id] = tmp_matrix
                            tmp_molecule_difference_matrix[rotation_A_id,rotation_B_id] = -permutedims(tmp_matrix,[2,1])
                        end

                    end
                end

            else
                
                # Matrix is not symmetric
                Threads.@threads for rotation_A_id in 1:molecules[molecule_A_id].Nrotations
                    # Get one coordinate matrix
                    structure_2 = molecules[molecule_A_id].coordinates_rotated[rotation_A_id]
                
                    for rotation_B_id in 1:molecules[molecule_B_id].Nrotations
                        # Get second coordinate matrix
                        structure_1 = molecules[molecule_B_id].coordinates_rotated[rotation_B_id]

                        # Define blank matrix
                        tmp_matrix = Array{Vector{Float64}}(undef, molecules[molecule_B_id].Natoms, molecules[molecule_A_id].Natoms)
            
                        # Calculate its values
                        # This matrix is non-symmetric
                        for atom_A_id in 1:molecules[molecule_A_id].Natoms
                            for atom_B_id in 1:molecules[molecule_B_id].Natoms
                                tmp_matrix[atom_B_id, atom_A_id] = lattice.inversevectors * (structure_2[:,atom_A_id] - structure_1[:,atom_B_id])
                            end
                        end

                        # Store it in the final matrix
                        tmp_molecule_difference_matrix[rotation_B_id, rotation_A_id] = tmp_matrix

                    end
                end

            end

            # Store the final results
            if molecule_A_id == molecule_B_id
                rotation_difference_matrices[molecule_B_id, molecule_B_id] = tmp_molecule_difference_matrix
            else
                rotation_difference_matrices[molecule_B_id, molecule_A_id] = tmp_molecule_difference_matrix
                rotation_difference_matrices[molecule_A_id, molecule_B_id] = -permutedims(tmp_molecule_difference_matrix,[2,1])
                # Row and column are swaped here as we transposed the original matrix!
                for row_id in 1:molecules[molecule_B_id].Nrotations
                    for column_id in 1:molecules[molecule_A_id].Nrotations
                        rotation_difference_matrices[molecule_A_id, molecule_B_id][column_id,row_id] = permutedims(rotation_difference_matrices[molecule_A_id, molecule_B_id][column_id,row_id],[2,1])
                    end
                end
                
            end

        end
    end

    # Return the result
    return rotation_difference_matrices

end

# This function generates the event list for every gridpoint within the unitcell
# The eventlist is structured by the label of the gridpoint and possible rotations of the molecule
# Each eventlist entry is a vector: 1st value is the affected gridpoint; ever other value is the affected orientation
function generate_event_list_per_unique_gridpoint(Nmolecules, molecules, Ngrids, grids, lattice, events, 
                unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices)
    

    # As reminder: For any distance between two atoms: Take the diffrence of the unitcell gridpoins + two translation vectors + rotations difference. At the end convert to cartesian coordinates.
    # u = unit_cell_gridpoints_difference[gridA,gridB][pointA, GridpointsMapping[1,pointB]]
    # a = translation_distance_vectors[GridpointsMapping[2,483]+1, GridpointsMapping[3,483]+1]
    # b = translation_distance_vectors[GridpointsMapping[2,121]+1, GridpointsMapping[3,121]+1]
    # c = rotation_difference_matrices[moleculeA,molecueB][1,100][2,1]
    # diff = u + b - a + c
    # fulllattice (transvectors) * diff

    # General
    overlap_2d = events.overlap2d


    # Generate the array to store all information
    # All levels are allocated just to make sure that entry 1, 2, etc has a well defined meaning (first, second, third, gridpoint) even if the point is actually never affected
    # Levels: Type of Molecule (Vector) -> Type of Grid (Vector) -> Unique gridpoints (Vector) -> Rotation (Vector) -> Affected points (Matrix: Molek, Grid, Point, Rot)
    Affected_Points_Rotations = Vector{Vector{Vector{Vector{Matrix{Int64}}}}}(undef,Nmolecules)
    for molecule_id in 1:Nmolecules
        Affected_Points_Rotations[molecule_id] = Vector{Vector{Vector{Matrix{Int64}}}}(undef,Ngrids)
        for grid_id in 1:Ngrids
            Affected_Points_Rotations[molecule_id][grid_id] = Vector{Vector{Matrix{Int64}}}(undef, grids[grid_id].Nuniquepoints)
            for unique_grid_id in 1:grids[grid_id].Nuniquepoints
                Affected_Points_Rotations[molecule_id][grid_id][unique_grid_id] = Vector{Matrix{Int64}}(undef, molecules[molecule_id].Nrotations)
                for rotation_id in 1:molecules[molecule_id].Nrotations
                    Affected_Points_Rotations[molecule_id][grid_id][unique_grid_id][rotation_id] = Matrix{Int64}(undef, 4 ,0)
                end
            end
        end
    end
    

    # Prepare the prescreening based on the molecular radius
    molecule_radius = Matrix{Float64}(undef, Nmolecules, Nmolecules)
    for molecule_1_id in 1:Nmolecules
        for molecule_2_id in 1:Nmolecules
            molecule_radius[molecule_1_id, molecule_2_id] = molecules[molecule_1_id].maxradius + molecules[molecule_2_id].maxradius
        end
    end


    # Now calculate all entries of the affected points array
    # Use the same for loop structure
    for molecule_id in 1:Nmolecules
        for grid_id in 1:Ngrids
            
            # Check whether the molecule (molecule_id) can be present on this grid (grid_id)
            if ! any(value -> value == grid_id, molecules[molecule_id].grids)
                continue
            end
            
            for point_id in 1:grids[grid_id].Nuniquepoints
                for rotation_id in 1:molecules[molecule_id].Nrotations

                    # These loops define the starting molecule
                    # Molecule type: molecule_id
                    # on Grid type: grid_id
                    # on unique gridpoint: point_id
                    # with rotation: rotation_id
                    # Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id]

                    
                    # More loops to define the Target molecule
                    for Tmolecule_id in 1:Nmolecules
                        for Tgrid_id in 1:Ngrids

                            # Check whether the molecule (Tmolecule_id) can be present on this grid (Tgrid_id)
                            if ! any(value -> value == Tgrid_id, molecules[Tmolecule_id].grids)
                                continue
                            end

                            for Tpoint_id in 1:grids[Tgrid_id].Npoints

                                # Skip
                                # if both grid points are equal
                                if grid_id == Tgrid_id && point_id == Tpoint_id
                                    continue
                                end

                                # Prescreening: Distance between gridpoints vs. molecular radius
                                # the translation_distance_vectors for the "point_id" is not needed as this is always a unique point and the distance vector therefor euqal to zero
                                mapping = grids[Tgrid_id].mapping[:,Tpoint_id]
                                fractional_points_difference_vector = unit_cell_gridpoints_difference[grid_id,Tgrid_id][point_id, mapping[1]] + 
                                                                translation_distance_vectors[mapping[2]+1, mapping[3]+1]
                                cartesian_points_difference_vector = apply_pbc_to_fractional_coordinates(fractional_points_difference_vector, lattice.transvectors)

                                if overlap_2d
                                    points_difference = norm(cartesian_points_difference_vector[1:2])
                                else
                                    points_difference = norm(cartesian_points_difference_vector)
                                end
                                if points_difference > molecule_radius[molecule_id, Tmolecule_id]
                                    continue
                                end
                                
                                # Generate an tmp matrix to store the next results
                                tmp_overlap_results = Matrix{Int64}(undef, 4 ,0)
                                Noverlaps = 0
                                
                                for Trotation_id in 1:molecules[Tmolecule_id].Nrotations

                                    # These loops define the target molecule
                                    # Molecule type: Tmolecule_id
                                    # on Grid type: Tgrid_id
                                    # on gridpoint: Tpoint_id
                                    # with rotation: Trotation_id

                                    # Check for overlap
                                    overlap = check_for_atom_overlap(fractional_points_difference_vector, points_difference, 
                                                                molecules, molecule_id, rotation_id, Tmolecule_id, Trotation_id,
                                                                lattice, rotation_difference_matrices, overlap_2d)


                                    # Are these two molecules overlapping?
                                    if overlap
                                        Noverlaps += 1
                                        tmp_overlap_results = hcat(tmp_overlap_results, [Tmolecule_id, Tgrid_id, Tpoint_id, Trotation_id])
                                    end

                                end

                                # Check whether all rotations are affected
                                #if Noverlaps == molecules[Tmolecule_id].Nrotations
                                #    Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id] = hcat(Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id], [Tmolecule_id, Tgrid_id, Tpoint_id, -1])
                                #else
                                    Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id] = hcat(Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id], tmp_overlap_results)
                                #end

                            end
                        end
                    end


                end
            end
        end
    end


    # Return the result
    return Affected_Points_Rotations

end

# A function to rearrange the rate constants of all events to simplify later usage
function sort_rate_constants_per_event(Nmolecules, Ngrids, events)

    # Generate the default matrices with negative rate constant and radii
    rate_constants_info = rate_constants_struct()
    rate_constants_info.ads = Matrix{Float64}(undef, Nmolecules, Ngrids) 
    rate_constants_info.rot = Matrix{Float64}(undef, Nmolecules, Ngrids) 
    rate_constants_info.dif = Array{Float64}(undef, Nmolecules, Ngrids, Ngrids)
    rate_constants_info.dif_rad = Array{Float64}(undef, Nmolecules, Ngrids, Ngrids)
    rate_constants_info.con = Array{Float64}(undef, Nmolecules, Ngrids, Nmolecules)

    # Set all values to "-1"
    rate_constants_info.ads .= -1.0
    rate_constants_info.rot .= -1.0
    rate_constants_info.dif .= -1.0
    rate_constants_info.dif_rad .= -1.0
    rate_constants_info.con .= -1.0

    # Read the events and update these values
    for event_id in 1:events.Nadsorptions
        rate_constants_info.ads[events.adsorptions[event_id].molecule, events.adsorptions[event_id].grid] = events.adsorptions[event_id].weigth
    end
    for event_id in 1:events.Nrotations
        rate_constants_info.rot[events.rotations[event_id].molecule, events.rotations[event_id].grid] = events.rotations[event_id].weigth
    end
    for event_id in 1:events.Ndiffusions
        rate_constants_info.dif[events.diffusions[event_id].molecule, events.diffusions[event_id].grid_start, events.diffusions[event_id].grid_end] = events.diffusions[event_id].weigth 
        rate_constants_info.dif_rad[events.diffusions[event_id].molecule, events.diffusions[event_id].grid_start, events.diffusions[event_id].grid_end] = events.diffusions[event_id].radius 
    end
    for event_id in 1:events.Nconformers
        rate_constants_info.con[events.conformers[event_id].molecule_start, events.conformers[event_id].grid, events.conformers[event_id].molecule_end] = events.conformers[event_id].weigth 
    end

    # Return the result
    return rate_constants_info

end

# A function to derive a neighbour list for every unique grid point
function create_neighbour_list(Nmolecules, Ngrids, rate_constants_info, lattice, grids, events, unit_cell_gridpoints_difference, translation_distance_vectors)

    # General
    overlap_2d = events.overlap2d

    # Create neighbour list vector
    neighbour_list = Vector{Vector{Vector{Vector{Vector{Int64}}}}}(undef,0)
    
    # fyi: We don't make use of the symmetry here... 
    # Loop over start grid type
    for start_grid_id in 1:Ngrids

        # Create tmp vector
        start_point_vector = Vector{Vector{Vector{Vector{Int64}}}}(undef,0)
        
        # Loop over unique points (start grid point)
        for start_point_id in 1:grids[start_grid_id].Nuniquepoints

            # Create tmp vector
            molecule_vector = Vector{Vector{Vector{Int64}}}(undef,0)

            # Loop over molecule
            for molecule_id in 1:Nmolecules

                # Create tmp vector
                grid_type_vector = Vector{Vector{Int64}}(undef,0)

                # Loop over end grid types
                for end_grid_id in 1:Ngrids

                    # Check whether the rate constants is negative
                    if rate_constants_info.dif[molecule_id, start_grid_id, end_grid_id] <= 0.0
                        end_point_vector = Vector{Int64}(undef,0)
                        push!(grid_type_vector, end_point_vector)
                        continue
                    end

                    # Get the diffusion radius
                    diff_radius = rate_constants_info.dif_rad[molecule_id, start_grid_id, end_grid_id]
                    
                    # Create tmp vector
                    end_point_vector = Vector{Int64}(undef,0)
                    
                    # Loop over points (end grid points)
                    for end_point_id in 1:grids[end_grid_id].Npoints

                        # Skip the identical point
                        if start_point_id == end_point_id
                            continue
                        end

                        # Get the distance between the grid points
                        mapping = grids[end_grid_id].mapping[:,end_point_id]
                        fractional_points_difference_vector = unit_cell_gridpoints_difference[start_grid_id,end_grid_id][start_point_id, mapping[1]] + 
                                                                translation_distance_vectors[mapping[2]+1, mapping[3]+1]
                        cartesian_points_difference_vector = apply_pbc_to_fractional_coordinates(fractional_points_difference_vector, lattice.transvectors)

                        if overlap_2d
                            distance = norm(cartesian_points_difference_vector[1:2])
                        else
                            distance = norm(cartesian_points_difference_vector)
                        end

                        # Add the point to the matrix in case it is within the diffusion radius
                        if distance < diff_radius
                            # add point to list [grid type][uniquepoint]: [for Molecule][Gridtype][target points]
                            push!(end_point_vector, end_point_id)
                        end

                    end

                    # Add the end point vector to the gridtype vector
                    push!(grid_type_vector, end_point_vector)

                end

                # Add the grid type vector to the molecule vector
                push!(molecule_vector, grid_type_vector)

            end

            # Add the molecule vector to the start point vector
            push!(start_point_vector, molecule_vector)

        end

        # Add the start point vector to the final neighbour list
        push!(neighbour_list, start_point_vector)

    end

    # Return results
    return neighbour_list

end

# This function calculates the inital matrices needed in the kMC run
# Finally, the event list for every unitcell gridpoint is generated
# It first calculates all rotation structes and all distance & difference matrices
function rsa_initialization(Nmolecules, molecules, Ngrids, grids, lattice, events)

    # Create all arrays and perform the initial array calculations

    # Calculate the difference vectors for all points within the unitcell
    unit_cell_gridpoints_difference = calculate_unit_cell_gridpoints_difference_vectors(Ngrids, grids, lattice)

    # Calculate the distance vectors for every combination of translation vectors
    # Hint: Two values of the resulting matrix have to be substracted to get the actual difference vector as the "distance vector" is relative to the origin
    translation_distance_vectors = calculate_translation_distance_vectors(lattice)

    # Calculate the difference vectors between atoms for every possible combination of rotating two structures
    rotation_difference_matrices = calculate_rotation_difference_matrices(Nmolecules, molecules, lattice)

    # Calculate the sum of atomic radii --> currently not used as it is not speeding up the simulation
    #length_sums = calculate_length_sums(Nmolecules, molecules)

    # Generate a list for every gridpoint of the unitcell & rotation of the molecule
    # This list contains every affected gridpoint & rotation to update all affected gridpoints as the consequence of a kMC event
    Affected_Points_Rotations = generate_event_list_per_unique_gridpoint(Nmolecules, molecules, Ngrids, grids, lattice, events, 
                            unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices)

    # Use events to rearrange rate constants
    rate_constants_info = sort_rate_constants_per_event(Nmolecules, Ngrids, events)

    # Create a neighbour list for every unique gridpoint
    neighbour_list = create_neighbour_list(Nmolecules, Ngrids, rate_constants_info, lattice, grids, events, unit_cell_gridpoints_difference, translation_distance_vectors)

    # Return the obtained matrices
    return unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list

end
