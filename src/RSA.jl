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
using LinearAlgebra
using TimerOutputs
using TimerOutputs: NoTimerOutput

#
# General head section of any file within this module 
#

#
# Include statements
#
include("constants.jl")
include("math.jl")
include("io.jl")
include("hdf5.jl")
include("pbc.jl")
include("analysis_plotting.jl")
include("init.jl")

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

"""
Mutable struct to store all results of multiple RSA runs.

# Fields
- `randomseed`: Vector storing the random number of every step.
- `status`: Matrix storing the currently adsorbed molecules in columns (molecule type, grid type, gridpoint, rotation).
- `Nsteps`: Number of performed RSA steps.
- `Nevents`: Matrix (molecule type, grid type) counting the number of performed events. Every element is a vector counting the event types (ads, rot, dif, con).
- `stepinfo`: Matrix storing information for every RSA step in columns (Total number of possible events, Number of possible ads/rot/dif/con events, selected grid type, selected gridpoint, selected molecule type, selected event type, selected subevent, selected event, selected subevent).
- `seed_offset`: Number stating the entries in the status & stepinfo matrix obtained by initializing the surface within a restart run

# Technical fields
- `size`: Current size of internal matrices.
"""
@kwdef mutable struct rsa_run_results_struct
    randomseed::Vector{Float64} = Vector{Float64}(undef, 0)
    status::Matrix{Int64} = Matrix{Int64}(undef, 4, 0) # molecule type, on grid type, on gridpoint, in rotation
    Nsteps::Int64 = 0
    Nevents::Matrix{Vector{Int64}} = Matrix{Int64}(undef, 0, 0) # molecule type, grid type; Vector: 1: Ads, 2: Rot, 3: Dif, 4: Con
    size::Int64 = 0
    stepinfo::Matrix{Int64} = Matrix{Int64}(undef, 12, 0) # Nevents possible, Nads & Nrot & Ndif & Ncon possible, selected grid type, selected gridpoint, selected molecule, selected event type, selected sub event, selected event, selected event 2 
    seed_offset::Int64 = 0 # Number of entries in the status & stepinfo for the initialization of a restart
end

#
# Function definitions
#

# Release/block delta for one state transition, computed once and cached.
#
# The delta does not actually depend on where in the lattice the event happens. For a
# rotation or a conformer change the two orientations sit on the same grid point; for a
# diffusion they sit a fixed number of cells apart. So the delta is a function of the
# transition alone -- (event type, both orientations, relative cell offset) -- and can be
# computed once in coordinates relative to the origin cell and reused everywhere. Each
# use then only translates the delta columns.
#
# Returns the cached (release, block) pair in relative coordinates. Both are owned by
# the cache and must not be mutated by the caller.
function cached_affected_delta(event_buffer, delta_key, Affected_Points_Rotations,
                               self_molecule, self_grid, self_unique_point, self_rotation,
                               target_molecule, target_grid, target_unique_point, target_rotation,
                               offset_transx, offset_transy, grids, lattice)

    entry = get(event_buffer.delta_cache, delta_key, nothing)
    if entry !== nothing
        return entry
    end

    # Self is taken straight from the event list: its stored point indices are already
    # relative to the origin cell. The target is brought into the same frame by
    # translating it by the relative cell offset.
    self_relative = Affected_Points_Rotations[self_molecule][self_grid][self_unique_point][self_rotation]
    Ntarget = map_affected_grid_points_to_actual_points!(event_buffer.target_relative, Affected_Points_Rotations, target_molecule, target_grid, target_unique_point, target_rotation, grids, lattice, offset_transx, offset_transy)
    target_relative = @view event_buffer.target_relative[:, 1:Ntarget]

    # This runs once per distinct transition, so the copies below are cheap; they have to
    # be copies because the buffers are reused on the next miss.
    Nrelease, Nblock = delta_matrix!(event_buffer.release, event_buffer.block, event_buffer.keys_self, event_buffer.keys_target, self_relative, target_relative)
    entry = (event_buffer.release[:, 1:Nrelease], event_buffer.block[:, 1:Nblock])

    # Bounded cache: past the cap the delta is recomputed per use instead of stored, so
    # a system with an unexpectedly large transition space cannot exhaust memory.
    if length(event_buffer.delta_cache) < rsa_delta_cache_max_entries
        event_buffer.delta_cache[delta_key] = entry
    end

    return entry

end

# Update the event list based on the selected event
# Update: Number of possible events, 
function rsa_update_event_list!(rsa_gridpoints, selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2, Nmolecules, molecules, Ngrids, grids, lattice, Affected_Points_Rotations, rate_constants_info, event_buffer)

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

    # Information on the selected point and molecule
    self_grid = selected_grid_type
    self_point = selected_grid_point
    self_molecule = selected_molecule
    self_rotation = rsa_gridpoints[self_grid][self_point].adsorbed_rotation
    self_unique_point, self_transx, self_transy = grids[self_grid].mapping[1:3, self_point]
    
    # Information on the target based on the selected event
    # By default the target is the adsorbate itself; every event changes only part of it
    target_molecule = self_molecule
    target_grid = self_grid
    target_point = self_point
    target_unique_point = self_unique_point
    target_transx = self_transx
    target_transy = self_transy
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

    # Cell offset between adsorbate and target (only non-zero for a diffusion)
    # Reduced modulo the cell counts as in map_translation_to_gridpoint, which keeps the delta cache keys bounded
    offset_transx = mod(target_transx - self_transx, lattice.Ncellx)
    offset_transy = mod(target_transy - self_transy, lattice.Ncelly)
    
    # The stamp counter is used to identify points which need to be updated
    event_buffer.counter += 1
    visit_stamp = event_buffer.counter
    affected_stamp = event_buffer.affected_stamp
    diffusion_stamp = event_buffer.diffusion_stamp
    list_affected_points = event_buffer.affected_points
    list_diffusion_points = event_buffer.diffusion_points
    empty!(list_affected_points)
    empty!(list_diffusion_points)

    # Self and target point are added to the affected points list
    mark_point!(list_affected_points, affected_stamp, visit_stamp, self_grid, self_point)
    if selected_event_type == 3
        mark_point!(list_affected_points, affected_stamp, visit_stamp, target_grid, target_point)
    end

    # Get the release and block lists
    #
    # An adsorption has nothing to release and blocks everything its own orientation
    # covers, so it needs the mapped self list directly and no delta at all.
    if selected_event_type == 1

        Nblock = map_affected_grid_points_to_actual_points!(event_buffer.block, Affected_Points_Rotations, self_molecule, self_grid, self_unique_point, target_rotation, grids, lattice, self_transx, self_transy)
        release_matrix = @view event_buffer.release[:, 1:0]
        block_matrix = @view event_buffer.block[:, 1:Nblock]

    else

        if event_buffer.memoize_delta

            delta_key = (selected_event_type, self_molecule, self_grid, self_unique_point, self_rotation,
                         target_molecule, target_grid, target_unique_point, target_rotation,
                         offset_transx, offset_transy)
            release_relative, block_relative = cached_affected_delta(event_buffer, delta_key, Affected_Points_Rotations,
                                                                    self_molecule, self_grid, self_unique_point, self_rotation,
                                                                    target_molecule, target_grid, target_unique_point, target_rotation,
                                                                    offset_transx, offset_transy, grids, lattice)

            # The cached delta is in coordinates relative to the origin cell; only its
            # columns have to be translated, not both full affected-point lists.
            Nrelease = translate_affected_columns!(event_buffer.release, release_relative, grids, lattice, self_transx, self_transy)
            Nblock = translate_affected_columns!(event_buffer.block, block_relative, grids, lattice, self_transx, self_transy)
            release_matrix = @view event_buffer.release[:, 1:Nrelease]
            block_matrix = @view event_buffer.block[:, 1:Nblock]

        else

            # Legacy path, kept reachable via memoize_delta = false
            Nself = map_affected_grid_points_to_actual_points!(event_buffer.target_relative, Affected_Points_Rotations, self_molecule, self_grid, self_unique_point, self_rotation, grids, lattice, self_transx, self_transy)
            self_affected_points = copy(@view event_buffer.target_relative[:, 1:Nself])
            Ntarget = map_affected_grid_points_to_actual_points!(event_buffer.target_relative, Affected_Points_Rotations, target_molecule, target_grid, target_unique_point, target_rotation, grids, lattice, target_transx, target_transy)
            target_affected_points = @view event_buffer.target_relative[:, 1:Ntarget]
            Nrelease, Nblock = delta_matrix!(event_buffer.release, event_buffer.block, event_buffer.keys_self, event_buffer.keys_target, self_affected_points, target_affected_points)
            release_matrix = @view event_buffer.release[:, 1:Nrelease]
            block_matrix = @view event_buffer.block[:, 1:Nblock]

        end

    end

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
        rsa_gridpoints[release_matrix[2, release_id]][release_matrix[3, release_id]].Nblocked_rotation[release_matrix[1, release_id]][release_matrix[4, release_id]] -= 1
        if rsa_gridpoints[release_matrix[2, release_id]][release_matrix[3, release_id]].Nblocked_rotation[release_matrix[1, release_id]][release_matrix[4, release_id]] == 0
            rsa_gridpoints[release_matrix[2, release_id]][release_matrix[3, release_id]].bool_rotation[release_matrix[1, release_id]][release_matrix[4, release_id]] = true
        end
        
        # Add the point to the list of affected points
        mark_point!(list_affected_points, affected_stamp, visit_stamp, release_matrix[2, release_id], release_matrix[3, release_id])

    end
    for block_id in axes(block_matrix,2)

        # Update the points
        rsa_gridpoints[block_matrix[2, block_id]][block_matrix[3, block_id]].Nblocked_rotation[block_matrix[1, block_id]][block_matrix[4, block_id]] += 1
        rsa_gridpoints[block_matrix[2, block_id]][block_matrix[3, block_id]].bool_rotation[block_matrix[1, block_id]][block_matrix[4, block_id]] = false

        # Add the point to the list of affected points
        mark_point!(list_affected_points, affected_stamp, visit_stamp, block_matrix[2, block_id], block_matrix[3, block_id])

    end

    # Update the reduced rotation vector and the number of free rotations
    for affected_point in list_affected_points
        for molecule_id in 1:Nmolecules
            rsa_gridpoints[affected_point[1]][affected_point[2]].reduced_rotation[molecule_id] = findall(x->x==true, rsa_gridpoints[affected_point[1]][affected_point[2]].bool_rotation[molecule_id])
            rsa_gridpoints[affected_point[1]][affected_point[2]].Nfree_rotation[molecule_id] = length(rsa_gridpoints[affected_point[1]][affected_point[2]].reduced_rotation[molecule_id])
        end
    end


    #
    # Step 3 - Update the bool_diffusion vector for all grid points
    # Here, we have to consider all diffusion targets of the self and target point, and all affected points 
    # So simply loop over list_affected_points as all these points are included there
    # 

    # Create a list of all points to update
    for affected_point in list_affected_points
        
        affected_grid = affected_point[1]
        affected_point_id = affected_point[2]
        affected_struct = rsa_gridpoints[affected_grid][affected_point_id]

        # Add the affected point to the list
        if affected_struct.bool_occupied == true
            mark_point!(list_diffusion_points, diffusion_stamp, visit_stamp, affected_grid, affected_point_id)
        end

        # Add diffusion targets to the list
        # Here I work with the origin_diffusion instead of the target_diffusion as diffusion events might be possible only in one direction
        for molecule_id in 1:Nmolecules
            origin_diffusion = affected_struct.origin_diffusion[molecule_id]
            for diffusion_grid_id in 1:Ngrids
                origin_points = origin_diffusion[diffusion_grid_id]
                for diffusion_point_id in eachindex(origin_points)

                    origin_point = origin_points[diffusion_point_id]
                    origin_struct = rsa_gridpoints[diffusion_grid_id][origin_point]
                    if origin_struct.bool_occupied == true
                        if origin_struct.adsorbed_molecule == molecule_id
                            mark_point!(list_diffusion_points, diffusion_stamp, visit_stamp, diffusion_grid_id, origin_point)
                        end
                    end

                end
            end
        end
        
    end

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

        # The affected-points matrix of this adsorbate is fixed for the loops below
        dif_affected_points = Affected_Points_Rotations[dif_ads_molecule][dif_grid][dif_unique_point][dif_rotation]
        
        # Loop over grids
        for grid_id in 1:Ngrids

            # Loop over points
            for point_id in eachindex(dif_target_points[grid_id])

                # Get the target grid point
                target_point = dif_target_points[grid_id][point_id]

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

                            if contains_column(dif_affected_points, dif_ads_molecule, grid_id, corrected_point, dif_rotation)
                                rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id][point_id] = true
                                continue
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

            end

        end

    end


    # Update the reduced diffusion vector and the number of free diffusion targets
    for diffusion_point in list_diffusion_points
        dif_grid, dif_point = diffusion_point
        dif_ads_molecule = rsa_gridpoints[dif_grid][dif_point].adsorbed_molecule
        for grid_id in 1:Ngrids
            rsa_gridpoints[dif_grid][dif_point].reduced_diffusion[dif_ads_molecule][grid_id] = findall(x->x==true, rsa_gridpoints[dif_grid][dif_point].bool_diffusion[dif_ads_molecule][grid_id])
            rsa_gridpoints[dif_grid][dif_point].Ntarget_diffusion[dif_ads_molecule][grid_id] = length(rsa_gridpoints[dif_grid][dif_point].reduced_diffusion[dif_ads_molecule][grid_id])
        end
    end


    #
    # Step 4 - Update rate constants and number of possible events
    #
    #println("Diffusion targets points list")
    #println(list_diffusion_points)
    for element in list_affected_points
        update_rate_constants_and_events!(rsa_gridpoints, element[1], element[2], Nmolecules, Ngrids, rate_constants_info, event_buffer)
    end
    for element in list_diffusion_points
        if affected_stamp[element[1]][element[2]] != visit_stamp
            update_rate_constants_and_events!(rsa_gridpoints, element[1], element[2], Nmolecules, Ngrids, rate_constants_info, event_buffer)
        end
    end


end

# A function to update the rate constants and number of events for a given grid point
function update_rate_constants_and_events!(rsa_gridpoints, gridtype, gridpoint, Nmolecules, Ngrids, rate_constants_info, event_buffer)

    # Get the grid point to update
    point = rsa_gridpoints[gridtype][gridpoint]

    # Counter values before the update, so the running totals can be adjusted by the
    # difference instead of being re-summed over the whole grid every step
    previous_nads = point.Nads
    previous_ndif = point.Ndif
    previous_nrot = point.Nrot
    previous_ncon = point.Ncon
    previous_nevents = point.Nevents

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
            cum_sum = 0.0
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
            cum_sum = 0.0
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
    cum_sum = 0.0
    for molecule_id in 1:Nmolecules
        cum_sum += point.Trate_molec[molecule_id]
        point.Cumulative_rate_molecules[molecule_id] = cum_sum
    end

    # Update cummulate rate constants (adsorption only)
    cum_sum = 0.0
    for molecule_id in 1:Nmolecules
        cum_sum += point.Trate_event[molecule_id][1]
        point.Cumulative_rate_molecules_ads[molecule_id] = cum_sum
    end

    # Update total rate constants
    point.Trateconst = point.Cumulative_rate_molecules[end]
    point.Trateconst_ads = point.Cumulative_rate_molecules_ads[end]

    # Update number of total events
    point.Nevents = point.Nads + point.Nrot + point.Ndif + point.Ncon 

    # Mirror the two rate constants into the dense arrays and fold the counter changes
    # into the running totals. Integer arithmetic, so the totals stay exact.
    event_buffer.trate[gridtype][gridpoint] = point.Trateconst
    event_buffer.trate_ads[gridtype][gridpoint] = point.Trateconst_ads
    event_buffer.total_nads += point.Nads - previous_nads
    event_buffer.total_ndif += point.Ndif - previous_ndif
    event_buffer.total_nrot += point.Nrot - previous_nrot
    event_buffer.total_ncon += point.Ncon - previous_ncon
    event_buffer.total_nevents += point.Nevents - previous_nevents

    # Return nothing - gridpoints are updated in place
    return nothing

end

# A function to calculate the rateconstant as sum over all gridpoints
# In addition a vector storing the cumulative rate constants is generated
function calculate_total_rateconstant!(rate_constant_buffers, rsa_gridpoints, Ngrids, grids, force_adsorption, event_buffer)

    # Default rate constants
    total_grid_rate_constant, cumulative_grid_rate_constants, cumulative_points_rate_constants = rate_constant_buffers
    total_rate_constant = 0.0

    # First loop over grids
    for grid_id in 1:Ngrids

        # The dense mirror of the per-point rate constant. Reading it instead of
        # dereferencing every gridpoint_struct turns this pass into a contiguous scan;
        # the summation order is unchanged (ascending point, then ascending grid), so the
        # accumulated values are bit-identical to summing the struct fields directly.
        point_rates = force_adsorption == true ? event_buffer.trate_ads[grid_id] : event_buffer.trate[grid_id]

        # Defaults
        rate_constant = 0.0
        cumulative_points_rates = cumulative_points_rate_constants[grid_id]

        # Second loop over gridpoints of this grid
        @inbounds for point_id in 1:grids[grid_id].Npoints
            rate_constant += point_rates[point_id]
            cumulative_points_rates[point_id] = rate_constant
        end

        # Update grid rate constants
        total_rate_constant += rate_constant
        total_grid_rate_constant[grid_id] = rate_constant
        cumulative_grid_rate_constants[grid_id] = total_rate_constant
        
    end

    # Events possible are maintained incrementally by update_rate_constants_and_events!
    # and refreshed once per run, so no per-step pass over the grid is needed
    total_events_possible = force_adsorption == true ? event_buffer.total_nads : event_buffer.total_nevents

    # Return result
    return total_rate_constant, cumulative_grid_rate_constants, total_grid_rate_constant, cumulative_points_rate_constants, total_events_possible, event_buffer.total_nads, event_buffer.total_nrot, event_buffer.total_ndif, event_buffer.total_ncon

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
function perform_rsa_step!(rsa_gridpoints, Ngrids, grids, Nmolecules, molecules, lattice, Affected_Points_Rotations, rate_constants_info, rsa_run_results, force_adsorption, rate_constant_buffers, event_buffer, timer)

    # Increase the step counter
    rsa_run_results.Nsteps += 1

    # Create total rate constants
    total_rate_constant, cumulative_grid_rate_constants, total_grid_rate_constant, cumulative_points_rate_constants, total_events_possible, ads_events_possible, rot_events_possible, dif_events_possible, con_events_possible = @timeit timer "Rate constant calculation" calculate_total_rateconstant!(rate_constant_buffers, rsa_gridpoints, Ngrids, grids, force_adsorption, event_buffer)
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
        total_rate_constant, cumulative_grid_rate_constants, total_grid_rate_constant, cumulative_points_rate_constants, total_events_possible, ads_events_possible, rot_events_possible, dif_events_possible, con_events_possible = @timeit timer "Rate constant calculation 2" calculate_total_rateconstant!(rate_constant_buffers, rsa_gridpoints, Ngrids, grids, force_adsorption, event_buffer)
    end


    # Stop the RSA run in case no reaction event is possible
    if total_events_possible == 0
        push!(rsa_run_results.randomseed, -1)
        rsa_run_results.stepinfo, rsa_run_results.size = @timeit timer "Step info bookkeeping" fill_preallocated_status_matrix(rsa_run_results.size, rsa_run_results.stepinfo, rsa_run_results.Nsteps, total_events_possible, ads_events_possible, rot_events_possible, dif_events_possible, con_events_possible, 0, 0, 0, 0, 0, 0, 0)
        return
    end

    # Select event
    random_number, selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2 = @timeit timer "Event selection" select_rsa_event(total_rate_constant, cumulative_grid_rate_constants, total_grid_rate_constant, cumulative_points_rate_constants, rsa_gridpoints, Ngrids, grids, Nmolecules, rate_constants_info, force_adsorption)

    # Store information
    # Collect the random numbers
    @timeit timer "Random number bookkeeping" push!(rsa_run_results.randomseed, random_number)
    # Count the selected events
    rsa_run_results.Nevents[selected_molecule, selected_grid_type][selected_event_type] += 1
    # Store all step information
    rsa_run_results.stepinfo, rsa_run_results.size = @timeit timer "Step info bookkeeping" fill_preallocated_status_matrix(rsa_run_results.size, rsa_run_results.stepinfo, rsa_run_results.Nsteps, total_events_possible, ads_events_possible, rot_events_possible, dif_events_possible, con_events_possible, selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2)

    # Update the status matrix
    @timeit timer "Adsorbate status bookkeeping" begin
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
    # timeit end
    end

    # Update event list
    @timeit timer "Event list update" rsa_update_event_list!(rsa_gridpoints, selected_grid_type, selected_grid_point, selected_molecule, selected_event_type, selected_subevent, selected_event, selected_event_2, Nmolecules, molecules, Ngrids, grids, lattice, Affected_Points_Rotations, rate_constants_info, event_buffer)

    # Return additional information for IO etc
    #return 
    
end

# Invoke all routines for a complete rsa run
function perform_rsa_run!(rsa_gridpoints, Ngrids, grids, Nmolecules, molecules, lattice, events, Affected_Points_Rotations, rate_constants_info, event_buffer, timer; seed::Matrix{Int64}=Matrix{Int64}(undef, 4, 0))

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

    # Allocate a buffer to calculate the total rate constant (total_grid_rate_constant), total cummulative rate constant (cumulative_grid_rate_constants),
    # and the cummilative rate constant for ever grid points (cumulative_points_rate_constants)
    rate_constant_buffers = (zeros(Ngrids), zeros(Ngrids), [Vector{Float64}(undef, grids[grid_id].Npoints) for grid_id in 1:Ngrids])
 
    # Rebuild the dense rate mirrors and running counter totals for this run. Covers
    # both a freshly initialised grid and one restored by reset_gridpoints!.
    refresh_rate_mirror!(event_buffer, rsa_gridpoints, Ngrids, grids)

    # Seed the calculation based on the current origin seed
    # This functions does nothing in case the seed is empty (therefore only used for restart runs)
    seed_grid!(rsa_gridpoints, seed, Ngrids, grids, Nmolecules, molecules, lattice, Affected_Points_Rotations, rate_constants_info, rsa_run_results, event_buffer, timer)
    rsa_run_results.seed_offset = size(seed, 2)

    # Perform rsa steps until no adsorption is possible any more
    while true

        # Perform RSA step
        @timeit timer "Steps" perform_rsa_step!(rsa_gridpoints, Ngrids, grids, Nmolecules, molecules, lattice, Affected_Points_Rotations, rate_constants_info, rsa_run_results, force_adsorption, rate_constant_buffers, event_buffer, timer)

        # Break conditions:
        # No further events possible
        if rsa_run_results.stepinfo[1, rsa_run_results.Nsteps] == 0
            break
        end

        # Break conditions:
        # Reached maximum number of RSA steps
        # Steps of the seeding are not included in the break condition
        if events.break_steps == true
            if rsa_run_results.Nsteps - rsa_run_results.seed_offset >= events.steps
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
- `memoize_delta`: Cache the release/block delta per state transition instead of recomputing it every step (default `true`). 

# Return values
Returns input information and results as structs in the follwing order:
- `rsa_results`: Vector storing all results of the RSA simulations. Every vector entry is a rsa\\_run\\_results_struct struct containing all information for a **single** RSA run. 
- `Nmolecules`: Total number of molecule types used.
- `molecules`: Vector storing all information concerning the molecule types. Every vector entry is a molecule_struct struct containing all information for a **single** molecule type.
- `Ngrids`: Total number of grid types used.
- `grids`: Vector storing all information concerning the grid types. Every vector entry is a grid_struct struct containing all information for a **single** grid type.
- `lattice`: A lattice_struct containing all information of the used lattice.
- `events`: A events_struct containing all information of the events and general RSA settings.


"""
function perform_multiple_rsa_runs(NRuns, inputfile_path; timer::Union{TimerOutput,NoTimerOutput}=NoTimerOutput(), hdf5=false, memoize_delta::Bool=true)

    #
    # RSA initialization
    #

    # TimerOutputs
    @timeit timer "Initialization" begin

    # Input read input files
    Nmolecules, molecules, Ngrids, grids, lattice, events = read_input(inputfile_path)

    # TimerOutputs
    end

    # Is this a restart run or the first run
    if events.restart_flag == true
        # Continue as restart run
        return perform_multiple_rsa_runs_restart(NRuns, timer, Nmolecules, molecules, Ngrids, grids, lattice, events, memoize_delta)
    else
        # Continue as first run
        return perform_multiple_rsa_runs_first(NRuns, inputfile_path, timer, hdf5, Nmolecules, molecules, Ngrids, grids, lattice, events, memoize_delta)
    end

end

# The main function controlling multiple rsa runs starting from scratch
function perform_multiple_rsa_runs_first(NRuns, inputfile_path, timer, hdf5, Nmolecules, molecules, Ngrids, grids, lattice, events, memoize_delta::Bool=true)

    # TimerOutputs
    @timeit timer "Initialization" begin

    # Set the current generation
    generation = 1

    # Check that the eventlist is reasonable
    check_eventlist!(molecules, grids, events)

    # Create HDF5 output file
    if hdf5 == true
        hdf5_file = create_hdf5_output_file(inputfile_path)
        write_hdf5_input_information(hdf5_file, generation, Nmolecules, molecules, Ngrids, grids, lattice, events)
    end

    # Generate all relevant matrices
    unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list = @timeit timer "Event initialization" rsa_initialization(Nmolecules, molecules, Ngrids, grids, lattice, events)
    
    # Add to the HDF5 file
    if hdf5 == true
        write_preparation_information(hdf5_file, generation, Ngrids, grids, lattice, Nmolecules, molecules, unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list)
    end

    # Allocate all matrices to store information for every run
    rsa_results = Vector{rsa_run_results_struct}(undef, NRuns)

    # Build the initial grid once and reuse it across all runs
    rsa_gridpoints = @timeit timer "Grid initialization" grid_initialization(Ngrids, grids, Nmolecules, molecules, lattice, rate_constants_info, neighbour_list)

    # Allocate the event-list event buffers
    event_buffer = rsa_event_buffer_initialization(Ngrids, grids, Affected_Points_Rotations, memoize_delta)

    # TimerOutputs
    end
    
    #
    # Perform rsa runs
    #

    # TimerOutput
    @timeit timer "Simulations" begin
    
    # The execution is currently not parallelized as this version of the grid reset is not thread safe. The grid reset needs to be reimplemented in a thread safe manner to allow parallel execution of the runs.
    for run_id in 1:NRuns

        # Restore the clean initial grid state before every run except the first
        if run_id != 1
            @timeit timer "Grid reset" reset_gridpoints!(rsa_gridpoints, Ngrids, Nmolecules, molecules, rate_constants_info)
        end

        # Perform the RSA run
        rsa_results[run_id] = @timeit timer "Runs" perform_rsa_run!(rsa_gridpoints, Ngrids, grids, Nmolecules, molecules, lattice, events, Affected_Points_Rotations, rate_constants_info, event_buffer, timer)

    end

    # TimerOutputs
    end

    # TimerOutput
    @timeit timer "Output storage" begin

    # Add results to HDF5 file
    if hdf5 == true
        write_rsa_results(hdf5_file, generation, NRuns, rsa_results)
        #println("All information are stored in the following HDF5 file:")
        #println(hdf5_file)
    end

    # TimerOutputs
    end

    #
    # Return all results
    #
    return rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events

end

# The main function controlling multiple rsa runs based on a restart calculation
function perform_multiple_rsa_runs_restart(NRuns, timer, Nmolecules, molecules, Ngrids, grids, lattice, events, memoize_delta::Bool=true)

    # TimerOutputs
    @timeit timer "Initialization" begin

    # Read the hdf5 file
    origin_rsa_results, origin_Nmolecules, origin_molecules, origin_Ngrids, origin_grids, origin_lattice, origin_events = read_hdf5_output_file(events.restart_file, events.restart_generation)

    # Resolve the requested origin runs
    resolve_restart_runs!(events, origin_rsa_results)

    # Check the molecule, grid, and lattice information
    validate_restart_compatibility(Nmolecules, molecules, Ngrids, grids, lattice, origin_Nmolecules, origin_molecules, origin_Ngrids, origin_grids, origin_lattice)
    
    # Set the current generation
    hdf5_file = events.restart_file
    generation = read_hdf5_number_of_generations(hdf5_file) + 1

    # Update HDF5 output file
    create_hdf5_restart_subgroups(hdf5_file, generation)
    reset_hdf5_number_of_generations(hdf5_file, generation)
    
    # Update the grids an adsorbate can be located on based on the seed
    seed_adsorbate_grids!(molecules, events.restart_runs, origin_rsa_results)

    # Check that the eventlist is reasonable
    check_eventlist!(molecules, grids, events)
    
    # Generate all relevant matrices
    unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list = rsa_initialization(Nmolecules, molecules, Ngrids, grids, lattice, events)
    
    # Add to the HDF5 file
    write_hdf5_input_information(hdf5_file, generation, Nmolecules, molecules, Ngrids, grids, lattice, events)
    write_preparation_information(hdf5_file, generation, Ngrids, grids, lattice, Nmolecules, molecules, unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list)

    # Get the total number of runs
    Ntotal = length(events.restart_runs) * NRuns

    # Allocate all matrices to store information for every run
    rsa_results = Vector{rsa_run_results_struct}(undef, Ntotal)

    # Build the initial grid once and reuse it across all runs
    rsa_gridpoints = @timeit timer "Grid initialization" grid_initialization(Ngrids, grids, Nmolecules, molecules, lattice, rate_constants_info, neighbour_list)

    # Allocate the event-list event buffers
    event_buffer = rsa_event_buffer_initialization(Ngrids, grids, Affected_Points_Rotations, memoize_delta)

    # TimerOutputs
    end

    #
    # Perform rsa runs
    #

    # TimerOutput
    @timeit timer "Simulations" begin

    # The execution is currently not parallelized as this version of the grid reset is not thread safe. The grid reset needs to be reimplemented in a thread safe manner to allow parallel execution of the runs.
    run_id = 0

    # First loop over the selected origins
    for origin_id in events.restart_runs

        # Get the seed of the origin run
        seed_status = origin_rsa_results[origin_id].status
        
        # Second loop over NRuns
        for _ in 1:NRuns

            # Increase run_id
            run_id += 1

            # Restore the clean initial grid state before every run except the first
            if run_id != 1
                @timeit timer "Grid reset" reset_gridpoints!(rsa_gridpoints, Ngrids, Nmolecules, molecules, rate_constants_info)
            end

            # Perform the RSA run
            rsa_results[run_id] = @timeit timer "Runs" perform_rsa_run!(rsa_gridpoints, Ngrids, grids, Nmolecules, molecules, lattice, events, Affected_Points_Rotations, rate_constants_info, event_buffer, timer; seed = seed_status)

        end
    
    end

    # TimerOutputs
    end

    # TimerOutput
    @timeit timer "Output" begin
    
    # Add results to HDF5 file
    write_rsa_results(hdf5_file, generation, Ntotal, rsa_results)

    # TimerOutputs
    end

    #
    # Return all results
    #
    return rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events

end

end # module RSA