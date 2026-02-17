"""
# RSA.jl

A *Random Sequential Adsorption* (RSA) package for the modelling of molecular adsorption events within area-selective deposition. 

## Features
* Random adsorption of adsorbates on a surface grid
* Support of multiple adsorbates and surface grids
* Support of diffusion and rotation events
* Support of adsorbate conversion events
* Analysis of surface coverage and effective gap size
* Multithreading support for running multiple simulations simultaneously

## Documentation
The complete documentation is available online: https://github.com/Tonner-Zech-Group/RSA
"""
module RSA
#
# Using statements
#
using DataFrames
using CSV
using LinearAlgebra
using Plots
using Printf
using Dates
using HDF5

#
# General head section of any file within this module 
#

#
# Include statements
#
include("constants.jl")
include("math.jl")
include("statebased.jl")
include("io.jl")
include("evaluation.jl")
include("pbc.jl")
include("analysis_plotting.jl")

#
# Export statements
#
export perform_multiple_rsa_runs

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

#= This is no longer used as it heavily affects the runtime of RSA steps 
@kwdef mutable struct rsa_timings_struct
    RSA_update_rateconstants::Float64 = 0.0
    RSA_select_events::Float64 = 0.0
    RSA_update_event_list::Float64 = 0.0
    RSA_compile_information::Float64 = 0.0
end
=#

"""
Mutable struct to store all timings of multiple RSA runs.

# Fields
- `Processing_Input`: Time for reading and processing input files.
- `Calculation_gridpoint_distance`: Time for calculating distances between all grid points.
- `Calculation_translation_distances`: Time for calculating all translation distances.
- `Calculation_rotation_distances`: Time for calculating distances between all rotations.
- `Calculation_affected_points`: Time used to find all affected points for every event.
- `Calculation_overlap`: Time used to judge the overlap of two molecules (also included in previous timing).
- `Calculation_rateconstants`: Time used to rearrange weigths of the events.
- `Calculation_neighbour_list`: Time to calculate a neighbour list for every grid point.
- `RSA_runs`: Time spend in RSA runs.
"""
@kwdef mutable struct timings_struct
    Processing_Input::Float64 = 0.0
    Calculation_gridpoint_distance::Float64 = 0.0
    Calculation_translation_distances::Float64 = 0.0
    Calculation_rotation_distances::Float64 = 0.0
    #Calculation_length_sums::Float64 = 0.0 currently not used
    Calculation_affected_points::Float64 = 0.0
    Calculation_overlap::Float64 = 0.0
    Calculation_rateconstants::Float64 = 0.0
    Calculation_neighbour_list::Float64 = 0.0
    RSA_runs::Float64 = 0.0
    #RSA::Vector{rsa_timings_struct} = Vector{rsa_timings_struct}(undef,0)
end

"""
Mutable struct to store all results of multiple RSA runs.

# Fields
- `randomseed`: Vector storing the random number of every step.
- `status`: Matrix storing the currently adsorbed molecules in columns (molecule type, grid type, gridpoint, rotation).
- `Nsteps`: Number of performed RSA steps.
- `Nevents`: Matrix (molecule type, grid type) counting the number of performed events. Every element is a vector counting the event types (ads, rot, dif, con).
- `stepinfo`: Matrix storing information for every RSA step in columns (Total number of possible events, Number of possible ads/rot/dif/con events, selected grid type, selected gridpoint, selected molecule type, selected event type, selected subevent, selected event, selected subevent).

# Technical fields
- `stepsize`: Factor to increase internal matrices. 
- `size`: Current size of internal matrices.
"""
@kwdef mutable struct rsa_run_results_struct
    randomseed::Vector{Float64} = Vector{Float64}(undef, 0)
    status::Matrix{Int64} = Matrix{Int64}(undef, 4, 0) # molecule type, on grid type, on gridpoint, in rotation
    Nsteps::Int64 = 0
    Nevents::Matrix{Vector{Int64}} = Matrix{Int64}(undef, 0, 0) # molecule type, grid type; Vector: 1: Ads, 2: Rot, 3: Dif, 4: Con
    stepsize::Int64 = 100_000
    size::Int64 = 0
    stepinfo::Matrix{Int64} = Matrix{Int64}(undef, 12, 0) # Nevents possible, Nads & Nrot & Ndif & Ncon possible, selected grid type, selected gridpoint, selected molecule, selected event type, selected sub event, selected event, selected event 2 
end

@kwdef mutable struct rate_constants_struct
    ads::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)         # [Molecule, Gridtype]: rate constant
    rot::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)         # [Molecule, Gridtype]: rate constant
    dif::Array{Float64} = Array{Float64}(undef, 0, 0, 0)        # [Molecule, Gridtype, Gridtype]: rate constant
    dif_rad::Array{Float64} = Array{Float64}(undef, 0, 0, 0)    # [Molecule, Gridtype, Gridtype]: diffusion radius
    con::Array{Float64} = Array{Float64}(undef, 0, 0, 0)        # [Molecule, Gridtype, Molecule]: rate constant
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
            target_diffusion = deepcopy(neighbour_list[grid_id][unique_point])

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
            pointsvec[points_id] = deepcopy(gridpoint_struct(Nevents, Nads, Ndif, Nrot, Ncon, Trateconst, Trateconst_ads, Trate_molec, Cumulative_rate_molecules, Cumulative_rate_molecules_ads, Trate_event, Cumulative_rate_events, Trate_diffusion, Cumulative_rate_diffusions, Trate_conformer, Cumulative_rate_conformers, target_diffusion, bool_diffusion, reduced_diffusion, Ntarget_diffusion, origin_diffusion, info_event, bool_occupied, adsorbed_molecule, adsorbed_rotation, bool_rotation, reduced_rotation, Nfree_rotation, Nblocked_rotation))

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

# Old version; still usefull to manually check whether two molecules overlap
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
                unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, timings)
    

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
                                    start_time = time()
                                    overlap = check_for_atom_overlap(fractional_points_difference_vector, points_difference, 
                                                                molecules, molecule_id, rotation_id, Tmolecule_id, Trotation_id,
                                                                lattice, rotation_difference_matrices, overlap_2d)

                                    end_time = time()
                                    timings.Calculation_overlap += end_time - start_time
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
    return Affected_Points_Rotations, timings

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
function rsa_initialization(Nmolecules, molecules, Ngrids, grids, lattice, events, timings)

    # Create all arrays and perform the initial array calculations

    # Calculate the difference vectors for all points within the unitcell
    start_time = time()
    unit_cell_gridpoints_difference = calculate_unit_cell_gridpoints_difference_vectors(Ngrids, grids, lattice)
    end_time = time()
    timings.Calculation_gridpoint_distance = end_time-start_time
    @printf "Timings - Calculation of gridpoint distances [in s]: %.2f \n" timings.Calculation_gridpoint_distance

    # Calculate the distance vectors for every combination of translation vectors
    # Hint: Two values of the resulting matrix have to be substracted to get the actual difference vector as the "distance vector" is relative to the origin
    start_time = time()
    translation_distance_vectors = calculate_translation_distance_vectors(lattice)
    end_time = time()
    timings.Calculation_translation_distances = end_time-start_time
    @printf "Timings - Calculation of translation vectors distances [in s]: %.2f \n" timings.Calculation_translation_distances

    # Calculate the difference vectors between atoms for every possible combination of rotating two structures
    start_time = time()
    rotation_difference_matrices = calculate_rotation_difference_matrices(Nmolecules, molecules, lattice)
    end_time = time()
    timings.Calculation_rotation_distances = end_time-start_time
    @printf "Timings - Calculation of rotation differences [in s]: %.2f \n" timings.Calculation_rotation_distances

    # Calculate the sum of atomic radii --> currently not used as it is not speeding up the simulation
    #start_time = time()
    #length_sums = calculate_length_sums(Nmolecules, molecules)
    #end_time = time()
    #timings.Calculation_length_sums = end_time-start_time

    # Generate a list for every gridpoint of the unitcell & rotation of the molecule
    # This list contains every affected gridpoint & rotation to update all affected gridpoints as the consequence of a kMC event
    start_time = time()
    Affected_Points_Rotations, timings = generate_event_list_per_unique_gridpoint(Nmolecules, molecules, Ngrids, grids, lattice, events, 
                            unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, timings)
    end_time = time()
    timings.Calculation_affected_points = end_time-start_time
    @printf "Timings - Calculation of affected gridpoints per event [in s]: %.2f \n" timings.Calculation_affected_points

    # Use events to rearrange rate constants
    start_time = time()
    rate_constants_info = sort_rate_constants_per_event(Nmolecules, Ngrids, events)
    end_time = time()
    timings.Calculation_rateconstants = end_time-start_time
    @printf "Timings - Calculation of rate constants [in s]: %.2f \n" timings.Calculation_rateconstants

    # Create a neighbour list for every unique gridpoint
    start_time = time()
    neighbour_list = create_neighbour_list(Nmolecules, Ngrids, rate_constants_info, lattice, grids, events, unit_cell_gridpoints_difference, translation_distance_vectors)
    end_time = time()
    timings.Calculation_neighbour_list = end_time-start_time
    @printf "Timings - Calculation of neighbour lists [in s]: %.2f \n" timings.Calculation_neighbour_list

    # Return the obtained matrices
    return unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list, timings

end

# Update the event list based on the selected event
# Update: Number of possible events, 
function rsa_update_event_list!(rsa_gridpoints, selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2, Nmolecules, molecules, Ngrids, grids, lattice, Affected_Points_Rotations, rate_constants_info)

     # Strategy:
     # + Update the bool_rotation vector (storing the number of accessible rotations) all the time
     # + The bool_diffusion vector (storing accessible diffusion sites) is only updated for occupied grid points
     #      This demands that in case a molecule occupies a new gridpoint all diffusion targets are checked for their accessibility
     # + For rate constant calculations:
     #    - take care that only events fitting the present molecule are calculated --> set everything else to zero
     #    - do not simply "count" true events in the bool vectors (they are updated even for molecules which are currently not present)
    
    
    #
    # Step 0 - Collecting information
    # Get a list of all points to update and check
    # 

    #start_time = time() 

    # Information on the selected point and molecule
    self_grid = selected_grid_type
    self_point = selected_grid_point
    self_molecule = selected_molecule
    self_rotation = rsa_gridpoints[self_grid][self_point].adsorbed_rotation
    self_unique_point, self_transx, self_transy = grids[self_grid].mapping[1:3, self_point]
    
    # Information on the target point based on the selected event
    if selected_event_type == 1
        target_rotation = selected_event
    elseif selected_event_type == 2
        target_rotation = selected_event
    elseif selected_event_type == 3
        target_rotation = selected_event_2
        target_grid = selected_subevent
        target_point = selected_event
        target_unique_point, target_transx, target_transy = grids[target_grid].mapping[1:3, target_point]
    elseif selected_event_type == 4
        target_rotation = selected_event
        target_molecule = selected_subevent
    end

    # Self and target point are added to the affected points list
    list_affected_points = Set{Vector{Int64}}()
    push!(list_affected_points, [self_grid, self_point])
    if selected_event_type == 3
        push!(list_affected_points, [target_grid, target_point])
    end

    #end_time = time()
    #println("Time of initalization: " * string(end_time-start_time))
    #start_time = time()  

    # Get the affected points for self and targets
    # These points are directly mapped to the actual points
    if selected_event_type == 1
        self_affected_points = map_affected_grid_points_to_actual_points(Affected_Points_Rotations, self_molecule, self_grid, self_unique_point, target_rotation, grids, lattice, self_transx, self_transy)
        target_affected_points = Matrix{Int64}(undef,4,0)
    elseif selected_event_type == 2
        self_affected_points = map_affected_grid_points_to_actual_points(Affected_Points_Rotations, self_molecule, self_grid, self_unique_point, self_rotation, grids, lattice, self_transx, self_transy)
        target_affected_points = map_affected_grid_points_to_actual_points(Affected_Points_Rotations, self_molecule, self_grid, self_unique_point, target_rotation, grids, lattice, self_transx, self_transy)
    elseif selected_event_type == 3
        self_affected_points = map_affected_grid_points_to_actual_points(Affected_Points_Rotations, self_molecule, self_grid, self_unique_point, self_rotation, grids, lattice, self_transx, self_transy)
        target_affected_points = map_affected_grid_points_to_actual_points(Affected_Points_Rotations, self_molecule, target_grid, target_unique_point, target_rotation, grids, lattice, target_transx, target_transy)
    elseif selected_event_type == 4
        self_affected_points = map_affected_grid_points_to_actual_points(Affected_Points_Rotations, self_molecule, self_grid, self_unique_point, self_rotation, grids, lattice, self_transx, self_transy)
        target_affected_points = map_affected_grid_points_to_actual_points(Affected_Points_Rotations, target_molecule, self_grid, self_unique_point, target_rotation, grids, lattice, self_transx, self_transy)
    end

    # Transform the affected points from unique points to actual points
    #=
    for point_id in axes(self_affected_points,2)
        mapping_unique_point, mapping_transx, mapping_transy = grids[self_affected_points[2,point_id]].mapping[1:3, self_affected_points[3,point_id]]
        new_transx = mapping_transx + self_transx
        new_transy = mapping_transy + self_transy
        affected_point = map_translation_to_gridpoint(mapping_unique_point, grids[self_affected_points[2,point_id]].Nuniquepoints, new_transx, lattice.Ncellx, new_transy, lattice.Ncelly)
        self_affected_points[3,point_id] = affected_point
    end
    for point_id in axes(target_affected_points,2)
        if selected_event_type == 3
            mapping_unique_point, mapping_transx, mapping_transy = grids[target_affected_points[2,point_id]].mapping[1:3, target_affected_points[3,point_id]]
            new_transx = mapping_transx + target_transx
            new_transy = mapping_transy + target_transy
            affected_point = map_translation_to_gridpoint(mapping_unique_point, grids[target_affected_points[2,point_id]].Nuniquepoints, new_transx, lattice.Ncellx, new_transy, lattice.Ncelly)
            target_affected_points[3,point_id] = affected_point
        else
            mapping_unique_point, mapping_transx, mapping_transy = grids[target_affected_points[2,point_id]].mapping[1:3, target_affected_points[3,point_id]]
            new_transx = mapping_transx + self_transx
            new_transy = mapping_transy + self_transy
            affected_point = map_translation_to_gridpoint(mapping_unique_point, grids[target_affected_points[2,point_id]].Nuniquepoints, new_transx, lattice.Ncellx, new_transy, lattice.Ncelly)
            target_affected_points[3,point_id] = affected_point
        end
    end
    =#

    #return self_affected_points, target_affected_points

    # Create a delta list
    if selected_event_type != 1
        release_matrix, block_matrix = delta_matrix(self_affected_points, target_affected_points) 
    else
        release_matrix = Matrix{Int64}(undef, 4, 0)
        block_matrix = self_affected_points
    end
    
    #return release_matrix, block_matrix

    #end_time = time()
    #println("Time of delta matrix: " * string(end_time-start_time))
    #start_time = time()  

    #
    # Step 1 - Update local information
    # 

    # Update whether a molecule is adsorbed, which molecule is adsorbed, and in which rotation
    if selected_event_type == 1
        rsa_gridpoints[self_grid][self_point].bool_occupied = true
        rsa_gridpoints[self_grid][self_point].adsorbed_molecule = self_molecule
        rsa_gridpoints[self_grid][self_point].adsorbed_rotation = target_rotation
    elseif selected_event_type == 2
        rsa_gridpoints[self_grid][self_point].adsorbed_rotation = target_rotation
    elseif selected_event_type == 3
        # Self
        rsa_gridpoints[self_grid][self_point].bool_occupied = false
        rsa_gridpoints[self_grid][self_point].adsorbed_molecule = 0
        rsa_gridpoints[self_grid][self_point].adsorbed_rotation = 0
        # Target
        rsa_gridpoints[target_grid][target_point].bool_occupied = true
        rsa_gridpoints[target_grid][target_point].adsorbed_molecule = self_molecule
        rsa_gridpoints[target_grid][target_point].adsorbed_rotation = target_rotation
    elseif selected_event_type == 4
        rsa_gridpoints[self_grid][self_point].adsorbed_molecule = target_molecule
        rsa_gridpoints[self_grid][self_point].adsorbed_rotation = target_rotation
    end
    
    #end_time = time()
    #println("Update time of local information: " * string(end_time-start_time))
    #start_time = time()  

    #
    # Step 2 - Update possible rotations for all grid points
    # 
  
    # Update self and target point
    # Nblocked_rotation is not updated here. This is only necessary for diffusion events and done later on.
    if selected_event_type == 1
        rsa_gridpoints[self_grid][self_point].bool_rotation[self_molecule][target_rotation] = false
    elseif selected_event_type == 2
        rsa_gridpoints[self_grid][self_point].bool_rotation[self_molecule][target_rotation] = false
        rsa_gridpoints[self_grid][self_point].bool_rotation[self_molecule][self_rotation] = true
    elseif selected_event_type == 3
        rsa_gridpoints[target_grid][target_point].bool_rotation[self_molecule][target_rotation] = false
        rsa_gridpoints[self_grid][self_point].bool_rotation[self_molecule][self_rotation] = true
    elseif selected_event_type == 4
        rsa_gridpoints[self_grid][self_point].bool_rotation[target_molecule][target_rotation] = false
        rsa_gridpoints[self_grid][self_point].bool_rotation[self_molecule][self_rotation] = true
    end

    # Update all affected points
    for release_id in axes(release_matrix,2)
        
        # Update the point
        #=
        if release_matrix[4, release_id] == -1
            for rotation_id in 1:molecules[release_matrix[1, release_id]].Nrotations
                rsa_gridpoints[release_matrix[2, release_id]][release_matrix[3, release_id]].Nblocked_rotation[release_matrix[1, release_id]][rotation_id] -= 1
                if rsa_gridpoints[release_matrix[2, release_id]][release_matrix[3, release_id]].Nblocked_rotation[release_matrix[1, release_id]][rotation_id] == 0
                    rsa_gridpoints[release_matrix[2, release_id]][release_matrix[3, release_id]].bool_rotation[release_matrix[1, release_id]][rotation_id] = true
                end
            end
        else
        =#
            rsa_gridpoints[release_matrix[2, release_id]][release_matrix[3, release_id]].Nblocked_rotation[release_matrix[1, release_id]][release_matrix[4, release_id]] -= 1
            if rsa_gridpoints[release_matrix[2, release_id]][release_matrix[3, release_id]].Nblocked_rotation[release_matrix[1, release_id]][release_matrix[4, release_id]] == 0
                rsa_gridpoints[release_matrix[2, release_id]][release_matrix[3, release_id]].bool_rotation[release_matrix[1, release_id]][release_matrix[4, release_id]] = true
            end
        #end
        
        # Add the point to the list of affected points
        new_element = release_matrix[2:3, release_id]
        if new_element ∉ list_affected_points
            push!(list_affected_points, new_element)
        end

    end
    for block_id in axes(block_matrix,2)

        # Update the points
        #=
        if block_matrix[4, block_id] == -1
            for rotation_id in 1:molecules[block_matrix[1, block_id]].Nrotations
                rsa_gridpoints[block_matrix[2, block_id]][block_matrix[3, block_id]].Nblocked_rotation[block_matrix[1, block_id]][rotation_id] += 1
                rsa_gridpoints[block_matrix[2, block_id]][block_matrix[3, block_id]].bool_rotation[block_matrix[1, block_id]][rotation_id] = false
            end
        else
        =#
            rsa_gridpoints[block_matrix[2, block_id]][block_matrix[3, block_id]].Nblocked_rotation[block_matrix[1, block_id]][block_matrix[4, block_id]] += 1
            rsa_gridpoints[block_matrix[2, block_id]][block_matrix[3, block_id]].bool_rotation[block_matrix[1, block_id]][block_matrix[4, block_id]] = false
        #end

        # Add the point to the list of affected points
        new_element = block_matrix[2:3, block_id]
        if new_element ∉ list_affected_points
            push!(list_affected_points, new_element)
        end

    end

    #end_time = time()
    #println("Update time of possible rotations: " * string(end_time-start_time))
    #start_time = time()  

    #return list_affected_points

    #println("Affected points list")
    #println(list_affected_points)

    # Update the reduced rotation vector and the number of free rotations
    for affected_point in list_affected_points
        for molecule_id in 1:Nmolecules
            rsa_gridpoints[affected_point[1]][affected_point[2]].reduced_rotation[molecule_id] = findall(x->x==true, rsa_gridpoints[affected_point[1]][affected_point[2]].bool_rotation[molecule_id])
            rsa_gridpoints[affected_point[1]][affected_point[2]].Nfree_rotation[molecule_id] = length(rsa_gridpoints[affected_point[1]][affected_point[2]].reduced_rotation[molecule_id])
        end
    end

    #end_time = time()
    #println("Update time of reduced rotations: " * string(end_time-start_time))
    #start_time = time()  


    #
    # Step 3 - Update the bool_diffusion vector for all grid points
    # Here, we have to consider all diffusion targets of the self and target point, and all affected points 
    # So simply loop over list_affected_points as all these points are included there
    # 

    # Create a list of all points to update
    list_diffusion_points = Set{Vector{Int64}}()
    for affected_point in list_affected_points
        #affected_point = list_affected_points[:,affected_point_id]
        
        # Add the affected point to the list
        if rsa_gridpoints[affected_point[1]][affected_point[2]].bool_occupied == true
            if affected_point ∉ list_diffusion_points
                push!(list_diffusion_points, affected_point)
            end
        end

        # Add diffusion targets to the list
        # Here I work with the origin_diffusion instead of the target_diffusion as diffusion events might be possible only in one direction
        for molecule_id in 1:Nmolecules
            for diffusion_grid_id in 1:Ngrids
                for diffusion_point_id in eachindex(rsa_gridpoints[affected_point[1]][affected_point[2]].origin_diffusion[molecule_id][diffusion_grid_id])
                    
                    diffusion_point = [diffusion_grid_id, rsa_gridpoints[affected_point[1]][affected_point[2]].origin_diffusion[molecule_id][diffusion_grid_id][diffusion_point_id]]
                    if rsa_gridpoints[diffusion_point[1]][diffusion_point[2]].bool_occupied == true
                        if rsa_gridpoints[diffusion_point[1]][diffusion_point[2]].adsorbed_molecule == molecule_id
                            if diffusion_point ∉ list_diffusion_points
                                push!(list_diffusion_points, diffusion_point)
                            end
                        end
                    end

                end
            end
        end
        
    end

    #end_time = time()
    #println("Creation time of diffusion points list: " * string(end_time-start_time))
    #start_time = time()  

    #return list_diffusion_points

    # For every diffusion target update the bool_diffusion vector
    # Only updated for the molecule type currently adsorbed
    # Only checked for points in the list_affected_points (thats not working!)
    for diffusion_point in list_diffusion_points
        dif_grid, dif_point = diffusion_point
        dif_ads_molecule = rsa_gridpoints[dif_grid][dif_point].adsorbed_molecule
        dif_rotation = rsa_gridpoints[dif_grid][dif_point].adsorbed_rotation
        dif_target_points = rsa_gridpoints[dif_grid][dif_point].target_diffusion[dif_ads_molecule]
        dif_unique_point, dif_transx, dif_transy = grids[dif_grid].mapping[1:3, dif_point]
        
        # Loop over grids
        for grid_id in 1:Ngrids

            # Loop over points
            for point_id in eachindex(dif_target_points[grid_id])

                # Get the target grid point
                target_point = dif_target_points[grid_id][point_id]

                # Check whether the point is in the list of affected points --> as the target points are also changed in other stepts (and I don't update the diffusion for unoccupied steps) I have to check all targets here.
                #if [grid_id, target_point] in eachcol(list_affected_points)

                    # Set the default
                    rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = true

                    # Is the target point occupied?
                    if rsa_gridpoints[grid_id][target_point].bool_occupied == true
                        rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = false
                        continue
                    end

                    # New version
                    # Only consider the current rotation of the adsorbate
                    # Is the current rotation blocked
                    if rsa_gridpoints[grid_id][target_point].bool_rotation[dif_ads_molecule][dif_rotation] == false

                        # Is the rotation blocked only once
                        if rsa_gridpoints[grid_id][target_point].Nblocked_rotation[dif_ads_molecule][dif_rotation] == 1

                            # Is the rotation blocked by the current adsorbate
                            blocked_unique_point, blocked_transx, blocked_transy = grids[grid_id].mapping[1:3, target_point]
                            corrected_transx = blocked_transx - dif_transx
                            corrected_transy = blocked_transy - dif_transy
                            corrected_point = map_translation_to_gridpoint(blocked_unique_point, grids[grid_id].Nuniquepoints, corrected_transx, lattice.Ncellx, corrected_transy, lattice.Ncelly)

                            if [dif_ads_molecule, grid_id, corrected_point, dif_rotation] in eachcol(Affected_Points_Rotations[dif_ads_molecule][dif_grid][dif_unique_point][dif_rotation])
                                rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = true
                                continue
                            #elseif [dif_ads_molecule, grid_id, corrected_point, -1] in eachcol(Affected_Points_Rotations[dif_ads_molecule][dif_grid][dif_unique_point][dif_rotation])
                            #    rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = true
                            #    continue
                            else
                                rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = false
                                continue
                            end

                        else
                            
                            # This rotation is blocked multiple times; a diffusion is not possible
                            rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = false
                            continue

                        end

                    end

                    # Old version
                    # Check for a single rotation blocked by the current adsorbate and enable a combination of diffusion and rotation
                    # Problem: In the selection of the RSA any information regarding which rotation should be selected for a diffusion is missing!
                    # Is the number of free rotations zero
                    #= 
                    if rsa_gridpoints[grid_id][target_point].Nfree_rotation[dif_ads_molecule] == 0

                        # Check for every rotation how often it is blocked
                        # As Nfree rotation is zero, all elements are larger than 0
                        # Find all elements equal 1. If there is none, all elements are at least equal or larger than 2
                        elements_equal_1 = findall(x -> x == 1, rsa_gridpoints[grid_id][target_point].Nblocked_rotation[dif_ads_molecule])

                        if isempty(elements_equal_1)
                            rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = false
                            continue
                        else
                            # Default
                            rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = false

                            # Transform the blocked point to its unique point
                            # Take care that the translation to the dif_point remains correct
                            blocked_unique_point, blocked_transx, blocked_transy = grids[grid_id].mapping[1:3, target_point]
                            corrected_transx = blocked_transx - dif_transx
                            corrected_transy = blocked_transy - dif_transy
                            corrected_point = map_translation_to_gridpoint(blocked_unique_point, grids[grid_id].Nuniquepoints, corrected_transx, lattice.Ncellx, corrected_transy, lattice.Ncelly)

                            # Check for every point only blocked by one adsorbate, whether it is blocked by this adsorbate
                            for blocked_rotation_id in eachindex(elements_equal_1)

                                # Is this point in the list of affected points
                                if [dif_ads_molecule, grid_id, corrected_point, elements_equal_1[blocked_rotation_id]] in eachcol(Affected_Points_Rotations[dif_ads_molecule][dif_grid][dif_unique_point][dif_rotation])
                                    rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = true
                                    continue
                                elseif [dif_ads_molecule, grid_id, corrected_point, -1] in eachcol(Affected_Points_Rotations[dif_ads_molecule][dif_grid][dif_unique_point][dif_rotation])
                                    rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = true
                                    continue
                                end

                            end

                        end                       

                    end
                    =#

                #end

            end

        end

    end

    #end_time = time()
    #println("Update time of diffusion vector: " * string(end_time-start_time))
    #start_time = time()  

    # Update the reduced diffusion vector and the number of free diffusion targets
    for diffusion_point in list_diffusion_points
        dif_grid, dif_point = diffusion_point
        dif_ads_molecule = rsa_gridpoints[dif_grid][dif_point].adsorbed_molecule
        for grid_id in 1:Ngrids
            rsa_gridpoints[dif_grid][dif_point].reduced_diffusion[dif_ads_molecule][grid_id] = findall(x->x==true, rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id])
            rsa_gridpoints[dif_grid][dif_point].Ntarget_diffusion[dif_ads_molecule][grid_id] = length(rsa_gridpoints[dif_grid][dif_point].reduced_diffusion[dif_ads_molecule][grid_id])
        end
    end

    #end_time = time()
    #println("Update time of reduced diffusion: " * string(end_time-start_time))
    #start_time = time() 


    #
    # Step 4 - Update rate constants and number of possible events
    #
    #println("Diffusion targets points list")
    #println(list_diffusion_points)
    list_union = union(list_affected_points, list_diffusion_points) 
    for element in list_union
        rsa_gridpoints[element[1]][element[2]] = update_rate_constants_and_events(rsa_gridpoints, element[1], element[2], Nmolecules, Ngrids, rate_constants_info)
    end

    #end_time = time()
    #println("Update time of rate constants: " * string(end_time-start_time)) 

end

# A function to update the rate constants and number of events for a given grid point
function update_rate_constants_and_events(rsa_gridpoints, gridtype, gridpoint, Nmolecules, Ngrids, rate_constants_info)

    # Get the grid point to update
    point = rsa_gridpoints[gridtype][gridpoint]

    # First check whether a molecule is adsorbed on this point
    if point.bool_occupied == true

        # Update conformer changes
        point.Ncon = 0
        # Loop over target molecules
        for molecule_id in 1:Nmolecules
            if molecule_id == point.adsorbed_molecule
                point.Trate_conformer[molecule_id] = 0.0
            else
                if rate_constants_info.con[point.adsorbed_molecule, gridtype, molecule_id] < 0
                    point.Trate_conformer[molecule_id] = 0.0
                else
                    point.Trate_conformer[molecule_id] = point.Nfree_rotation[molecule_id] * rate_constants_info.con[point.adsorbed_molecule, gridtype, molecule_id]
                    point.Ncon += point.Nfree_rotation[molecule_id]
                end
            end
        end
        
        # Update cummulative rate constant
        cum_sum = 0.0
        for molecule_id in 1:Nmolecules
            cum_sum += point.Trate_conformer[molecule_id]
            point.Cumulative_rate_conformers[molecule_id] = cum_sum
        end

        # Update diffusions
        point.Ndif = 0
        # Loop over target gridtypes
        for grid_id in 1:Ngrids
            if rate_constants_info.dif[point.adsorbed_molecule, gridtype, grid_id] < 0
                point.Trate_diffusion[grid_id] = 0.0
            else
                point.Trate_diffusion[grid_id] = point.Ntarget_diffusion[point.adsorbed_molecule][grid_id] * rate_constants_info.dif[point.adsorbed_molecule, gridtype, grid_id]
                point.Ndif += point.Ntarget_diffusion[point.adsorbed_molecule][grid_id]
            end
        end

        # Update cummulative rate constants
        cum_sum = 0.0
        for grid_id in 1:Ngrids
            cum_sum += point.Trate_diffusion[grid_id]
            point.Cumulative_rate_diffusions[grid_id] = cum_sum
        end

        # Update Trate_event vector
        point.Nrot = 0
        for molecule_id in 1:Nmolecules
            if molecule_id == point.adsorbed_molecule
                point.Trate_event[molecule_id][1] = 0.0
                if rate_constants_info.rot[molecule_id, gridtype] < 0
                    point.Trate_event[molecule_id][2] = 0.0
                else
                    point.Trate_event[molecule_id][2] = point.Nfree_rotation[molecule_id] * rate_constants_info.rot[molecule_id, gridtype]
                    point.Nrot += point.Nfree_rotation[molecule_id]
                end
                point.Trate_event[molecule_id][3] = point.Cumulative_rate_diffusions[end]
                point.Trate_event[molecule_id][4] = point.Cumulative_rate_conformers[end]
            else
                point.Trate_event[molecule_id][:] .= 0.0
            end
        end

        # Update cummulative rate constant
        for molecule_id in 1:Nmolecules
            cum_sum = 0
            for event_id in 1:4
                cum_sum += point.Trate_event[molecule_id][event_id]
                point.Cumulative_rate_events[molecule_id][event_id] = cum_sum
            end
        end

        # Update number of possible adsorption events
        point.Nads = 0

    else
        
        # Update conformer changes
        point.Trate_conformer .= 0.0
        point.Cumulative_rate_conformers .= 0.0

        # Update diffusions
        point.Trate_diffusion .= 0.0
        point.Cumulative_rate_diffusions .= 0.0

        # Update Trate_event vector
        point.Nads = 0
        for molecule_id in 1:Nmolecules
            if rate_constants_info.ads[molecule_id, gridtype] < 0
                point.Trate_event[molecule_id][:] .= 0.0
            else
                point.Trate_event[molecule_id][1] = point.Nfree_rotation[molecule_id] * rate_constants_info.ads[molecule_id, gridtype]
                point.Trate_event[molecule_id][2:4] .= 0.0
                point.Nads += point.Nfree_rotation[molecule_id]
            end
        end

        # Update cummulative rate constant
        for molecule_id in 1:Nmolecules
            cum_sum = 0
            for event_id in 1:4
                cum_sum += point.Trate_event[molecule_id][event_id]
                point.Cumulative_rate_events[molecule_id][event_id] = cum_sum
            end
        end

        # Update number of possible events
        point.Ndif = 0
        point.Nrot = 0
        point.Ncon = 0

    end

    # Update total rate constants per molecule
    for molecule_id in 1:Nmolecules
        point.Trate_molec[molecule_id] = point.Cumulative_rate_events[molecule_id][end]
    end

    # Update cummulative rate per molecule
    cum_sum = 0
    for molecule_id in 1:Nmolecules
        cum_sum += point.Trate_molec[molecule_id]
        point.Cumulative_rate_molecules[molecule_id] = cum_sum
    end

    # Update cummulate rate constants (adsorption only)
    cum_sum = 0
    for molecule_id in 1:Nmolecules
        cum_sum += point.Trate_event[molecule_id][1]
        point.Cumulative_rate_molecules_ads[molecule_id] = cum_sum
    end

    # Update total rate constants
    point.Trateconst = point.Cumulative_rate_molecules[end]
    point.Trateconst_ads = point.Cumulative_rate_molecules_ads[end]

    # Update number of total events
    point.Nevents = point.Nads + point.Nrot + point.Ndif + point.Ncon 

    # Return the updated gridpoint
    return point

end

# A function to calculate the rateconstant as sum over all gridpoints
# In addition a vector storing the cumulative rate constants is generated
function calculate_total_rateconstant(rsa_gridpoints, Ngrids, grids, force_adsorption)

    # Default rate constants
    total_rate_constant = 0.0
    total_grid_rate_constant = zeros(Ngrids)
    cumulative_grid_rate_constants = zeros(Ngrids)
    cumulative_points_rate_constants = Vector{Vector{Float64}}(undef, Ngrids)

    # Default events possible
    total_events_possible = 0
    ads_events_possible = 0
    dif_events_possible = 0
    rot_events_possible = 0
    con_events_possible = 0

    # First loop over grids
    for grid_id in 1:Ngrids

        # Get the submatrix of the gridpoints
        subgrid_matrix = rsa_gridpoints[grid_id]

        # Defaults
        rate_constant = 0
        cumulative_points_rates = Vector{Float64}(undef, grids[grid_id].Npoints)

        # Second loop over gridpoints of this grid
        for point_id in 1:grids[grid_id].Npoints

            # Update rate constants
            if force_adsorption == true
                rate_constant += subgrid_matrix[point_id].Trateconst_ads
            else
                rate_constant += subgrid_matrix[point_id].Trateconst
            end
            cumulative_points_rates[point_id] = rate_constant
            
            # Update events possible
            if force_adsorption == true
                total_events_possible += subgrid_matrix[point_id].Nads
            else
                total_events_possible += subgrid_matrix[point_id].Nevents
            end
            ads_events_possible += subgrid_matrix[point_id].Nads
            dif_events_possible += subgrid_matrix[point_id].Ndif
            rot_events_possible += subgrid_matrix[point_id].Nrot
            con_events_possible += subgrid_matrix[point_id].Ncon

        end

        # Update grid rate constants
        total_rate_constant += rate_constant
        total_grid_rate_constant[grid_id] = rate_constant
        cumulative_points_rate_constants[grid_id] = cumulative_points_rates
        cumulative_grid_rate_constants[grid_id] = total_rate_constant
        
    end

    # Return result
    return total_rate_constant, cumulative_grid_rate_constants, total_grid_rate_constant, cumulative_points_rate_constants, total_events_possible, ads_events_possible, rot_events_possible, dif_events_possible, con_events_possible

end

# A function to select the gridpoint and element based on a random number
function select_rsa_event(total_rate_constant, cumulative_grid_rate_constants, total_grid_rate_constant, cumulative_points_rate_constants, rsa_gridpoints, Ngrids, grids, Nmolecules, rate_constants_info, force_adsorption)

    # Get a random number in the range (0,1]
    random_number = 1-rand()

    # The event is selected based on the hierachy: grid type -> grid point -> molecule -> event type -> event
    
    # Select the grid type
    grid_type_value = random_number * total_rate_constant
    selected_grid_type = 0
    for grid_type_id in 1:Ngrids
        if grid_type_value <= cumulative_grid_rate_constants[grid_type_id]
            selected_grid_type = grid_type_id
            break
        end
    end

    # Select the grid point
    if selected_grid_type > 1
        grid_point_value = grid_type_value - cumulative_grid_rate_constants[selected_grid_type - 1]
    else
        grid_point_value = grid_type_value
    end
    selected_grid_point = 0
    for point_id in 1:grids[selected_grid_type].Npoints
        if grid_point_value <= cumulative_points_rate_constants[selected_grid_type][point_id]
            selected_grid_point = point_id
            break
        end
    end

    # Select the molecule
    if selected_grid_point > 1
        molecule_value = grid_point_value - cumulative_points_rate_constants[selected_grid_type][selected_grid_point - 1]
    else
        molecule_value = grid_point_value 
    end
    if force_adsorption == false
        selected_molecule = 0
        for molecule_id in 1:Nmolecules
            if molecule_value <= rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_molecules[molecule_id]
                selected_molecule = molecule_id
                break
            end
        end
    elseif force_adsorption == true
        selected_molecule = 0
        for molecule_id in 1:Nmolecules
            if molecule_value <= rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_molecules_ads[molecule_id]
                selected_molecule = molecule_id
                break
            end
        end
    end

    # Select the event type
    if force_adsorption == false
        if selected_molecule > 1
            event_value = molecule_value - rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_molecules[selected_molecule - 1]
        else
            event_value = molecule_value
        end
        selected_event_type = 0
        for event_id in 1:4
            if event_value <= rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_events[selected_molecule][event_id]
                selected_event_type = event_id
                break
            end
        end
    elseif force_adsorption == true
        if selected_molecule > 1
            event_value = molecule_value - rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_molecules_ads[selected_molecule - 1]
        else
            event_value = molecule_value
        end
        selected_event_type = 1
    end
    
    # Select subevent (necessary of diffusion - target grid type - and conformer changes - target molecule.)
    # 3 -> dif
    # 4 -> con
    selected_subevent = 0
    subevent_value = 0.0
    if selected_event_type == 3
        subevent_value = event_value - rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_events[selected_molecule][2]
        for subevent_id in 1:Ngrids
            if subevent_value <= rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_diffusions[subevent_id]
                selected_subevent = subevent_id
                break
            end
        end 
    elseif selected_event_type == 4
        subevent_value = event_value - rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_events[selected_molecule][3]
        for subevent_id in 1:Nmolecules
            if subevent_value <= rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_conformers[subevent_id]
                selected_subevent = subevent_id
                break
            end
        end
    end

    # Select the event
    # 1 -> ads
    # 2 -> rot
    # 3 -> dif
    # 4 -> con
    selected_event = 0
    selected_event_2 = 0
    if selected_event_type == 1
        # To get the target rotation
        rate_constant_per_event = rate_constants_info.ads[selected_molecule, selected_grid_type]
        #event_value = random_number * rsa_gridpoints[selected_grid_type][selected_grid_point].Trate_event[selected_molecule][1]
        rotation_value = event_value
        selected_event_relative = ceil(Int64, rotation_value / rate_constant_per_event)
        selected_event = rsa_gridpoints[selected_grid_type][selected_grid_point].reduced_rotation[selected_molecule][selected_event_relative]
    elseif selected_event_type == 2
        # To get the target rotation
        rate_constant_per_event = rate_constants_info.rot[selected_molecule, selected_grid_type]
        #event_value = random_number * rsa_gridpoints[selected_grid_type][selected_grid_point].Trate_event[selected_molecule][2]
        rotation_value = event_value - rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_events[selected_molecule][1]
        selected_event_relative = ceil(Int64, rotation_value / rate_constant_per_event)
        selected_event = rsa_gridpoints[selected_grid_type][selected_grid_point].reduced_rotation[selected_molecule][selected_event_relative]
    elseif selected_event_type == 3
        # To get the target point
        rate_constant_per_event = rate_constants_info.dif[selected_molecule, selected_grid_type, selected_subevent]
        #event_value = random_number * rsa_gridpoints[selected_grid_type][selected_grid_point].Trate_diffusion[selected_subevent]
        if selected_subevent > 1
            target_point_value = subevent_value - rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_diffusions[selected_subevent - 1]
        else
            target_point_value = subevent_value
        end
        selected_event_relative = ceil(Int64, target_point_value / rate_constant_per_event)
        relative_point = rsa_gridpoints[selected_grid_type][selected_grid_point].reduced_diffusion[selected_molecule][selected_subevent][selected_event_relative]
        selected_event = rsa_gridpoints[selected_grid_type][selected_grid_point].target_diffusion[selected_molecule][selected_subevent][relative_point]
        # To get the target rotation
        # In the current version the rotation must stay the same
        #rate_constant_per_event = rate_constants_info.ads[selected_molecule, selected_subevent]
        #event_value = random_number * rsa_gridpoints[selected_subevent][selected_event].Trate_event[selected_molecule][1]
        #selected_event_relative = ceil(Int64, event_value / rate_constant_per_event)
        #selected_event_2 = rsa_gridpoints[selected_grid_type][selected_grid_point].reduced_rotation[selected_molecule][selected_event_relative]
        selected_event_2 = rsa_gridpoints[selected_grid_type][selected_grid_point].adsorbed_rotation
    elseif selected_event_type == 4
        # To get the target rotation
        rate_constant_per_event = rate_constants_info.con[selected_molecule, selected_grid_type, selected_subevent]
        if selected_subevent > 1
            rotation_value = subevent_value - rsa_gridpoints[selected_grid_type][selected_grid_point].Cumulative_rate_conformers[selected_subevent - 1]
        else
            rotation_value = subevent_value
        end
        #event_value = random_number * rsa_gridpoints[selected_grid_type][selected_grid_point].Trate_conformer[selected_subevent]
        selected_event_relative = ceil(Int64, rotation_value / rate_constant_per_event)
        selected_event = rsa_gridpoints[selected_grid_type][selected_grid_point].reduced_rotation[selected_subevent][selected_event_relative]
    end
    

    # IO
    #=
    println("Random number: "*string(random_number))
    println("Grid value: "*string(grid_type_value))
    println("Selected grid type: "*string(selected_grid_type))
    println("Point value: "*string(grid_point_value))
    println("Selected grid point: "*string(selected_grid_point))
    println("Molecule value: "*string(molecule_value))
    println("Selected molecule: "*string(selected_molecule))
    println("Event value: "*string(event_value))
    println("Selected event type: "*string(selected_event_type))
    println("Subevent value: "*string(subevent_value))
    println("Selected subevent: "*string(selected_subevent))
    println("Selected event: "*string(selected_event))
    println("Selected event 2: "*string(selected_event_2))
    =#

    # Return results
    return random_number, selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2

end

# Invoke all routines for a single rsa step
function perform_rsa_step!(rsa_gridpoints, Ngrids, grids, Nmolecules, molecules, lattice, Affected_Points_Rotations, rate_constants_info, rsa_run_results, force_adsorption)

    # Increase the step counter
    rsa_run_results.Nsteps += 1

    # Create total rate constants
    #start_time = time()
    total_rate_constant, cumulative_grid_rate_constants, total_grid_rate_constant, cumulative_points_rate_constants, total_events_possible, ads_events_possible, rot_events_possible, dif_events_possible, con_events_possible = calculate_total_rateconstant(rsa_gridpoints, Ngrids, grids, force_adsorption)
    #=
    println("Total rate constant: "*string(total_rate_constant))
    println("Cummulative grid rate constant: "*string(cumulative_grid_rate_constants))
    println("Total grid rate constant: "*string(total_grid_rate_constant))
    println("Total events possible: "*string(total_events_possible))
    =#

    # Check that in case of an forced adsorption an event is actually possible
    if force_adsorption == true && total_events_possible == 0
        # As an adsorption event can not be forced the rate constants are recalculated without forced adsorption
        force_adsorption = false
        total_rate_constant, cumulative_grid_rate_constants, total_grid_rate_constant, cumulative_points_rate_constants, total_events_possible, ads_events_possible, rot_events_possible, dif_events_possible, con_events_possible = calculate_total_rateconstant(rsa_gridpoints, Ngrids, grids, force_adsorption)
    end
    #end_time = time()
    #rsa_timings[1] += end_time-start_time


    # Stop the RSA run in case no reaction event is possible
    if total_events_possible == 0
        push!(rsa_run_results.randomseed, -1)
        rsa_run_results.stepinfo, rsa_run_results.size = fill_preallocated_status_matrix!(rsa_run_results.size, rsa_run_results.stepsize, rsa_run_results.stepinfo, rsa_run_results.Nsteps, total_events_possible, ads_events_possible, rot_events_possible, dif_events_possible, con_events_possible, 0, 0, 0, 0, 0, 0, 0)
        return
    end

    # Select event
    #start_time = time()
    random_number, selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2 = select_rsa_event(total_rate_constant, cumulative_grid_rate_constants, total_grid_rate_constant, cumulative_points_rate_constants, rsa_gridpoints, Ngrids, grids, Nmolecules, rate_constants_info, force_adsorption)
    #end_time = time()
    #rsa_timings[2] += end_time-start_time

    # Store information
    #start_time = time()
    # Collect the random numbers
    push!(rsa_run_results.randomseed, random_number)
    # Count the selected events
    rsa_run_results.Nevents[selected_molecule, selected_grid_type][selected_event_type] += 1
    # Store all step information
    rsa_run_results.stepinfo, rsa_run_results.size = fill_preallocated_status_matrix!(rsa_run_results.size, rsa_run_results.stepsize, rsa_run_results.stepinfo, rsa_run_results.Nsteps, total_events_possible, ads_events_possible, rot_events_possible, dif_events_possible, con_events_possible, selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2)

    # Update the status matrix
    if selected_event_type == 1
        rsa_run_results.status = hcat(rsa_run_results.status, [selected_molecule, selected_grid_type, selected_grid_point, selected_event])
    elseif selected_event_type == 2
        change_column = findfirst_column(rsa_run_results.status, [selected_molecule, selected_grid_type, selected_grid_point], 3)
        rsa_run_results.status[4,change_column] = selected_event
    elseif selected_event_type == 3
        change_column = findfirst_column(rsa_run_results.status, [selected_molecule, selected_grid_type, selected_grid_point], 3)
        rsa_run_results.status[2:4,change_column] = [selected_subevent, selected_event, selected_event_2]
    elseif selected_event_type == 4
        change_column = findfirst_column(rsa_run_results.status, [selected_molecule, selected_grid_type, selected_grid_point], 3)
        rsa_run_results.status[1:4,change_column] = [selected_subevent, selected_grid_type, selected_grid_point, selected_event]
    end
    #end_time = time()
    #rsa_timings[4] += end_time-start_time
    

    # Update event list
    #start_time = time()
    rsa_update_event_list!(rsa_gridpoints, selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2, Nmolecules, molecules, Ngrids, grids, lattice, Affected_Points_Rotations, rate_constants_info)
    #end_time = time()
    #rsa_timings[3] += end_time-start_time

    # Return additional information for IO etc
    #return 
    
end

# Invoke all routines for a complete rsa run (initialization is not included)
function perform_rsa_run(Ngrids, grids, Nmolecules, molecules, lattice, events, Affected_Points_Rotations, rate_constants_info, neighbour_list)

    # Initialization of the grid
    rsa_gridpoints = grid_initialization(Ngrids, grids, Nmolecules, molecules, lattice, rate_constants_info, neighbour_list)
    
    # Initialization of event counter
    non_adsorption_events = 0

    # Initialize forced adsorption
    force_adsorption = false

    # Struct to store all run information
    rsa_run_results = rsa_run_results_struct()
    rsa_run_results.Nevents = Matrix{Vector{Int64}}(undef, Nmolecules, Ngrids)
    for molecule_id in 1:Nmolecules
        for grid_id in 1:Ngrids
            rsa_run_results.Nevents[molecule_id, grid_id] = [0 for i in 1:4]
        end
    end
 
    # Perform rsa steps until no adsorption is possible any more
    while true

        # Perform RSA step
        perform_rsa_step!(rsa_gridpoints, Ngrids, grids, Nmolecules, molecules, lattice, Affected_Points_Rotations, rate_constants_info, rsa_run_results, force_adsorption)

        # Break conditions:
        # No further events possible
        if rsa_run_results.stepinfo[1, rsa_run_results.Nsteps] == 0
            break
        end

        # Break conditions:
        # Reached maximum number of RSA steps
        if events.break_steps == true
            if rsa_run_results.Nsteps >= events.steps
                break
            end
        end

        # Count number of consecutive non-adsorption events
        if rsa_run_results.stepinfo[9, rsa_run_results.Nsteps] == 1
            non_adsorption_events = 0
        else
            non_adsorption_events += 1
        end

        # Break conditions:
        # non-adsorption events reached convergence value
        if events.break_convergence == true
            if non_adsorption_events >= events.coverage_convergence
                break
            end
        end

        # Force adsorption in the next step
        force_adsorption = false
        if events.force_adsorption == true
            if non_adsorption_events >= events.Nforce_adsorption && rsa_run_results.stepinfo[2, rsa_run_results.Nsteps] > 0
                force_adsorption = true
            end
        end

    end

    # Reduce the stepinfo to the number of steps
    rsa_run_results.stepinfo = rsa_run_results.stepinfo[:,1:rsa_run_results.Nsteps]

    # Return all valuable information for IO and further analysis
    return rsa_run_results
    
end

"""

    perform_multiple_rsa_runs(NRuns::Integer, inputfile_path::String)
    perform_multiple_rsa_runs(NRuns::Integer, inputfile_path::String; hdf5::Bool = false)

Performs `NRuns` RSA simulations controlled by the input file found under the given path `inputfile_path`.

# Optional input
- `hdf5`: Flag to request the generation of a HDF5 file storing all information.

# Return values
Returns input information and results as structs in the follwing order:
- `rsa_results`: Vector storing all results of the RSA simulations. Every vector entry is a rsa\\_run\\_results_struct struct containing all information for a **single** RSA run. 
- `Nmolecules`: Total number of molecule types used.
- `molecules`: Vector storing all information concerning the molecule types. Every vector entry is a molecule_struct struct containing all information for a **single** molecule type.
- `Ngrids`: Total number of grid types used.
- `grids`: Vector storing all information concerning the grid types. Every vector entry is a grid_struct struct containing all information for a **single** grid type.
- `lattice`: A lattice_struct containing all information of the used lattice.
- `events`: A events_struct containing all information of the events and general RSA settings.
- `timings`: A timings_struct containing some timing information.


"""
function perform_multiple_rsa_runs(NRuns, inputfile_path; hdf5=false)

    #
    # RSA initialization
    #

    # Currently there are no restart features
    Nrun = 1

    # Generate structs
    timings = timings_struct()
    
    #= removed for better performance
    timings.RSA = Vector{rsa_timings_struct}(undef,NRuns)
    for rsa_rund_id in 1:NRuns
        timings.RSA[rsa_rund_id] = rsa_timings_struct()
    end
    =#

    # Input read input files
    start_time = time()
    Nmolecules, molecules, Ngrids, grids, lattice, events = read_input(inputfile_path)
    end_time = time()
    timings.Processing_Input = end_time-start_time
    @printf "Timings - Input file processing [in s]: %.2f \n" timings.Processing_Input

    # Create HDF5 output file
    if hdf5 == true
        hdf5_file = create_hdf5_output_file(inputfile_path)
        write_hdf5_input_information(hdf5_file, Nrun, Nmolecules, molecules, Ngrids, grids, lattice, events)
    end

    # Generate all relevant matrices
    unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list, timings = rsa_initialization(Nmolecules, molecules, Ngrids, grids, lattice, events, timings)
    
    # Add to the HDF5 file
    if hdf5 == true
        write_preparation_information(hdf5_file, Nrun, Ngrids, grids, lattice, Nmolecules, molecules, unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list)
    end

    # For debugging:
    #return unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list, timings

    # Allocate all matrices to store information for every run
    rsa_results = Vector{rsa_run_results_struct}(undef, NRuns)
    
    #
    # Perform rsa runs
    #
    println("Starting the RSA simulations ...")
    start_time = time()
    # The execution of kMC runs is parallized while each run is serial
    Threads.@threads for run_id in 1:NRuns

        # Perform the RSA run
        rsa_results[run_id] = perform_rsa_run(Ngrids, grids, Nmolecules, molecules, lattice, events, Affected_Points_Rotations, rate_constants_info, neighbour_list)

    end
    println("... done!")
    end_time = time()
    timings.RSA_runs = end_time-start_time
    @printf "Timings - RSA simulations [in s]: %.2f \n" timings.RSA_runs

    # Add results to HDF5 file
    if hdf5 == true
        write_hdf5_timings(hdf5_file, Nrun, timings)
        write_rsa_results(hdf5_file, Nrun, NRuns, rsa_results)
        println("All information are stored in the following HDF5 file:")
        println(hdf5_file)
    end


    #
    # Return all results
    #
    return rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events, timings

    # For debugging
    #return rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events, unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list, timings

end


end # module kMC