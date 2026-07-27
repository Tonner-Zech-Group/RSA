#
# General head section of any file within this module 
#

#
# Using statements
#
using Dates

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
"""
Mutable struct to store all information concerning **one** molecule.

# General fields
- `label`: Label of the molecule.
- `grids`: Vector storing all possible grids the molecule can reach.

# Fixpoint fields
- `fixpointtype`: Keyword defining the type of the used fixpoint.
- `fixpointatoms`: Vector storing the atoms used for the fixpoint.
- `fixpoint`: Vector storing the fixpoint.

# Rotation fields
- `rotationmodus`: Keyword selecting how rotations are generated.
- `rotationangle`: Value of rotation angle to generate all rotations.
- `rotationvalues`: Vector storing all explicit values to generate rotations.
- `Nrotations`: Total number of rotations.
- `rotations`: Vector storing the finally used rotation values.

# Structural fields
- `coordinate_path`: Path to the xyz coordinate file.
- `Natoms`: Total number of atoms.
- `elements`: Vector containing the elements of each atom.
- `elements_sorted`: Element vector sorted based on the distance to the fixpoint.
- `dimension`: Dimension of the given coordinates.
- `coordinates`: Matrix storing the initial xyz coordinates.
- `coordinates_sorted`: Coordinates sorted based on the distance to the fixpoint.
- `coordinates_rotated`: Vector of matrices storing all rotated structures.
- `maxradius`: Largest distance to the fixpoint including vdW radius of the atom.
- `radii`: Distance to the fixpoint of every atom.
- `radii_sorted`: Sorted distance to the fixpoint of every atom.
"""
@kwdef mutable struct molecule_struct
    # General infos
    label::Int64 = 0
    grids::Vector{Int64} = Vector{Int64}(undef, 0)
    
    # Define fixpoint type, used atoms, and position
    fixpointtype::String = "centroid"
    fixpointatoms::Vector{Int64} = Vector{Int64}(undef, 0)
    fixpoint::Matrix{Float64} = zeros(3,1)

    # Define rotation modus, rotation angle, and final rotations
    rotationmodus::String = "values"
    rotationangle::Float64 = 360.0
    rotationvalues::Vector{Float64} = Vector{Float64}(undef, 0)
    Nrotations::Int64 = 0
    rotations::Vector{Float64} = Vector{Float64}(undef, 0)

    # Structural information
    coordinate_path::String = ""
    Natoms::Int64 = 0
    elements::Vector{Int64} = Vector{Int64}(undef, 0)
    elements_sorted::Vector{Int64} = Vector{Int64}(undef, 0)
    dimension::Int64 = 3
    coordinates::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    coordinates_sorted::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    coordinates_rotated::Vector{Matrix{Float64}} = Vector{Matrix{Float64}}(undef, 0)
    maxradius::Float64 = 0.0
    radii::Vector{Float64} = Vector{Float64}(undef, 0)
    radii_sorted::Vector{Float64} = Vector{Float64}(undef, 0)

end

"""
Mutable struct to store all information concerning **one** grid.

# General fields
- `label`: Label of the grid.

# Structural fields
- `Nuniquepoints`: Number of uniques points within the unit cell.
- `uniquepoints`: Matrix storing the coordinates of the unique points.
- `dimension`: Dimension of the coordinates.
- `Npoints`: Total number of points in the supercell.
- `points`: Matrix storing the coordinates of all grid points.
- `mapping`: Matrix storing the mapping (unique point, x translation, y translation) of every grid point.
"""
@kwdef mutable struct grid_struct
    # General infos
    label::Int64 = 0

    # Structural information
    Nuniquepoints::Int64 = 0
    uniquepoints::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    dimension::Int64 = 3
    Npoints::Int64 = 0
    points::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    mapping::Matrix{Int64} = Matrix{Int64}(undef, 0, 0)
end

"""
Mutable struct to store all information concerning the lattice.

# Structural fields
- `vectors`: Vector spanning the unit cell.
- `dimension`: Dimension of the vectors.

# Translation fields
- `transx`: Number of translations along the first vector.
- `transy`: Number of translations along the second vector.
- `Ncellx`: Total number of cells along the first vector (is equal to transx + 1).
- `Ncelly`: Total number of cells along the second vector (is equal to transy + 1).
- `transvectors`: Lattice vectors of the supercell.
- `inversevectors`: Inverse lattice vectors of the supercell.
"""
@kwdef mutable struct lattice_struct
    # General information
    vectors::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    dimension::Int64 = 3

    # Translations
    transx::Int64 = 0
    transy::Int64 = 0
    Ncellx::Int64 = 1
    Ncelly::Int64 = 1
    transvectors::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    inversevectors::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
end
@kwdef mutable struct event_diffusion_struct
    molecule::Int64 = 1
    grid_start::Int64 = 1
    grid_end::Int64 = 1
    weigth::Float64 = 1.0
    radius::Float64 = 1.0
end
@kwdef mutable struct event_adsorption_struct
    molecule::Int64 = 1
    grid::Int64 = 1
    weigth::Float64 = 1.0
end
@kwdef mutable struct event_rotation_struct
    molecule::Int64 = 1
    grid::Int64 = 1
    weigth::Float64 = 1.0
end
@kwdef mutable struct event_conformer_change_struct
    molecule_start::Int64 = 1
    molecule_end::Int64 = 1
    grid::Int64 = 1
    weigth::Float64 = 1.0
end

"""
Mutable struct to store all information concerning the events and general RSA settings.

# General fields
- `steps`: Maximum number of steps to be performed in a RSA simulation.
- `break_steps`: Boolean to flag whether the max steps keyword is used.
- `coverage_convergence`: Maximum number of steps without adsorption for coverage convergence. 
- `break_convergence`: Boolean to flag whether the coverage convergence keyword is used.
- `Nforce_adsorption`: Number of non-adsorption steps after which an adsorption event is forced.
- `force_adsorption`: Boolean to flag whether the forced adsorption keyword is used.
- `overlap2d`: Boolean to flag whether 2D overlap is used.
- `overlap3d`: Boolean to flag whether 3D overlap is used.

# Restart fields
- `restart_flag`: Boolean to flag whether this is a restart run.
- `restart_generation`: Generation number of RSA runs the restart is based on.
- `restart_runs`: Vector storing the run numbers of RSA runs used in the restart.
- `restart_file`: Path to the restart hdf5 file.

# Event fields
- `Nevents`: Total number of defined events.
- `Nadsorptions`: Number of defined adsorption events.
- `adsorptions`: Vector storing the adsorption events.
- `Ndiffusions`: Number of defined diffusion events.
- `diffusions`: Vector storing the diffusion events.
- `Nrotations`: Number of defined rotation events.
- `rotations`: Vector storing the rotation events.
- `Nconformers`: Number of defined conformer or conversion events.
- `conformers`: Vector storing the conversion events.
"""
@kwdef mutable struct events_struct
    # General settings
    steps::Int64 = 0
    break_steps::Bool = false
    coverage_convergence::Int64 = 0
    break_convergence::Bool = false
    Nforce_adsorption::Int64 = 0
    force_adsorption::Bool = false
    overlap2d::Bool = true
    overlap3d::Bool = false

    # Restart settings
    restart_flag::Bool = false
    restart_generation::Int64 = 0
    restart_runs::Vector{Int64} = Vector{Int64}(undef, 0)
    restart_file::String = ""
    
    # Currently not used
    weigth_scale::Float64 = 0.0 # not used
    force_weigth_scale::Bool = false # not used

    # Event settings
    Nevents::Int64 = 0
    Nadsorptions::Int64 = 0
    adsorptions::Vector{event_adsorption_struct} = Vector{event_adsorption_struct}(undef, 0)
    Ndiffusions::Int64 = 0
    diffusions::Vector{event_diffusion_struct} = Vector{event_diffusion_struct}(undef, 0)
    Nrotations::Int64 = 0
    rotations::Vector{event_rotation_struct} = Vector{event_rotation_struct}(undef, 0)
    Nconformers::Int64 = 0
    conformers::Vector{event_conformer_change_struct} = Vector{event_conformer_change_struct}(undef, 0)
end

# General read function
# Call the function to read the input file
# Call the function to read xyz file
# Rearrange and complete input
function read_input(path::String)

    # Use the path to the general input file
    Nmolecules, molecules, Ngrids, grids, lattice, events = read_input_file(path)

    # Control that the dimensions of the lattice, grid, molecule, and overlap match
    if events.overlap3d
        if !isempty(findall(x->x!=3, lattice.dimension))
            println("The dimension of the lattice should be 3. (3 lattice vectors with 3 coordinates)")
            error("Input File Error")
        end
        if !isempty(findall(x->x!=3, getproperty.(grids, :dimension)))
            println("The dimension of the grids should be 3. (3 coordinates per gridpoint)")
            error("Input File Error")
        end
        if !isempty(findall(x->x!=3, getproperty.(molecules, :dimension)))
            println("The dimension of the molecules should be 3. (3 coordinates per atom)")
            error("Input File Error")
        end
    elseif events.overlap2d
        if !isempty(findall(x->x!=lattice.dimension, getproperty.(grids, :dimension)))
            println("The dimension of the lattice and grids do not match. Both need the same number of coordinates per lattice vector and gridpoint.")
            error("Input File Error")
        end
        if !isempty(findall(x->x!=lattice.dimension, getproperty.(molecules, :dimension)))
            println("The dimension of the lattice and molecules do not match. Both need the same number of coordinates per lattice vector and atom.")
            error("Input File Error")
        end
    end

    # Use the translations to generate all gridpoints for every grid
    for grid_id in 1:Ngrids
        # Get all points
        grids[grid_id].points, grids[grid_id].Npoints, grids[grid_id].mapping = 
                    replicate_gridpoints_with_translation(lattice.vectors, grids[grid_id].uniquepoints, lattice.transx, lattice.transy)
    end

    # For every molecule
    if events.overlap2d
        overlap_case = "2D"
    else
        overlap_case = "3D"
    end

    for molecule_id in 1:Nmolecules
        # Get the distance of every atom to the fixpoint
        molecules[molecule_id].maxradius, molecules[molecule_id].radii = get_largest_vdW_distance_to_point(molecules[molecule_id].Natoms, molecules[molecule_id].elements, molecules[molecule_id].coordinates, molecules[molecule_id].fixpoint, overlap_case)

        # Sort this distances from largest to smallest
        permutation = sortperm(molecules[molecule_id].radii, rev=true)
        molecules[molecule_id].radii_sorted = molecules[molecule_id].radii[permutation]

        # Sort coodinates and element vector based on this permutation
        molecules[molecule_id].elements_sorted = molecules[molecule_id].elements[permutation]
        molecules[molecule_id].coordinates_sorted = molecules[molecule_id].coordinates[:,permutation]

        # Create all rotated structures
        molecules[molecule_id].coordinates_rotated = create_all_rotations_centered(molecules[molecule_id].coordinates_sorted , molecules[molecule_id].Nrotations, molecules[molecule_id].rotations, molecules[molecule_id].fixpoint, lattice.dimension)

    end

    # Use the molecule and grid labels to correct the eventlist
    # The eventlist is using the molecules and grids in the same order as the corresponding vector
    for event_id in 1:events.Nadsorptions
        # Get the label of the current molecule
        tmp_molec_label = events.adsorptions[event_id].molecule
        tmp_grid_label = events.adsorptions[event_id].grid
        # Search for the molecule with this label
        new_molec_label = find_label(Nmolecules, molecules, tmp_molec_label)
        new_grid_label = find_label(Ngrids, grids, tmp_grid_label)
        # Replace old label
        events.adsorptions[event_id].molecule = new_molec_label
        events.adsorptions[event_id].grid = new_grid_label
    end
    for event_id in 1:events.Ndiffusions

        tmp_molec_label = events.diffusions[event_id].molecule
        new_molec_label = find_label(Nmolecules, molecules, tmp_molec_label)
        events.diffusions[event_id].molecule = new_molec_label
        
        tmp_grid_label = events.diffusions[event_id].grid_start
        new_grid_label = find_label(Ngrids, grids, tmp_grid_label)
        events.diffusions[event_id].grid_start = new_grid_label

        tmp_grid_label = events.diffusions[event_id].grid_end
        new_grid_label = find_label(Ngrids, grids, tmp_grid_label)
        events.diffusions[event_id].grid_end = new_grid_label
    end
    for event_id in 1:events.Nrotations
        tmp_molec_label = events.rotations[event_id].molecule
        new_molec_label = find_label(Nmolecules, molecules, tmp_molec_label)
        events.rotations[event_id].molecule = new_molec_label

        tmp_grid_label = events.rotations[event_id].grid
        new_grid_label = find_label(Ngrids, grids, tmp_grid_label)
        events.rotations[event_id].grid = new_grid_label
    end
    for event_id in 1:events.Nconformers
        tmp_molec_label = events.conformers[event_id].molecule_start
        new_molec_label = find_label(Nmolecules, molecules, tmp_molec_label)
        events.conformers[event_id].molecule_start = new_molec_label

        tmp_molec_label = events.conformers[event_id].molecule_end
        new_molec_label = find_label(Nmolecules, molecules, tmp_molec_label)
        events.conformers[event_id].molecule_end = new_molec_label

        tmp_grid_label = events.conformers[event_id].grid
        new_grid_label = find_label(Ngrids, grids, tmp_grid_label)
        events.conformers[event_id].grid = new_grid_label
    end

    # Check that one convergence criterium is selected
    if events.break_steps == false && events.break_convergence == false
        println("No convergence criterium selected!")
        println("Use either 'Steps' or 'Coverageconvergence' to define an endpoint for this simulation.")
        error("Input Conversion Error")
    end

    # Check that the event list is meaningful
    # Check 1: On which gridtype can a molecule adsorb
    for adsorption_id in 1:events.Nadsorptions
        molecule_id = events.adsorptions[adsorption_id].molecule
        grid_id = events.adsorptions[adsorption_id].grid
        if ! any(value -> value == grid_id, molecules[molecule_id].grids)
            push!(molecules[molecule_id].grids, grid_id) 
        end
    end

    # Check 2: Whether conformer change is possible on this grid
    for conformer_id in 1:events.Nconformers
        molecule_start_id = events.conformers[conformer_id].molecule_start
        molecule_end_id =events.conformers[conformer_id].molecule_end
        grid_id = events.conformers[conformer_id].grid

        # check whether this event is meaningful
        if any(value -> value == grid_id, molecules[molecule_start_id].grids)
            if ! any(value -> value == grid_id, molecules[molecule_end_id].grids)
                push!(molecules[molecule_end_id].grids, grid_id)
            end 
        else
            # Extension: Search conformer change events recursively
            # For now: Throw an error and ask the user to rearrange the event list
            println("I have problems reading events defined in the eventlist.")
            println("Problematic case: " * string(molecules[molecule_start_id].label) * " con " * string(molecules[molecule_end_id].label) * " " * string(grids[grid_id].label))
            println("Please define the event for the molecule reaching the initial grid type first.")
            error("Event List Error")
        end
    end

    # Check 3: To which gridtype a molecule can diffuse
    for diffusion_id in 1:events.Ndiffusions
        molecule_id = events.diffusions[diffusion_id].molecule
        grid_start_id = events.diffusions[diffusion_id].grid_start
        grid_end_id = events.diffusions[diffusion_id].grid_end

        # check whether this event is meaningful
        if any(value -> value == grid_start_id, molecules[molecule_id].grids)
            if ! any(value -> value == grid_end_id, molecules[molecule_id].grids)
                push!(molecules[molecule_id].grids, grid_end_id)
            end 
        else
            # Extension: Search diffusion events recursively
            # For now: Throw an error and ask the user to rearrange the event list
            println("I have problems reading events defined in the eventlist.")
            println("Problematic case: " * string(molecules[molecule_id].label) * " dif " * string(grids[grid_start_id].label) * " " * string(grids[grid_end_id].label))
            println("Please define the event for the molecule reaching the initial grid type first.")
            error("Event List Error")
        end
    end
    
    # Check 4: Rotations only on grid types the molecule can reach
    for rotation_id in 1:events.Nrotations
        molecule_id = events.rotations[rotation_id].molecule
        grid_id = events.rotations[rotation_id].grid
        if ! any(value -> value == grid_id, molecules[molecule_id].grids)
            println("I have problems reading events defined in the eventlist.")
            println("A molecule should not be able to rotate on a grid type it can never reach.")
            println("Problematic case: " * string(molecules[molecule_id].label) * " rot " * string(grids[grid_id].label))
            error("Event List Error")
        end
    end

    # Return everything
    return Nmolecules, molecules, Ngrids, grids, lattice, events
end

# A function to find a label in the known molecules or grids
function find_label(Nitems, structure, label)

    label_position = 0
    for item_id in 1:Nitems
        if structure[item_id].label == label
            label_position = item_id
            break
        end
    end

    if label_position == 0
        if typeof(structure) == Vector{molecule_struct}
            println("The following label was not found in the molecules struct: " * string(label))
            error("Input Conversion Error")
        elseif typeof(structure) == Vector{grid_struct}
            println("The following label was not found in the grids struct: " * string(label))
            error("Input Conversion Error")
        end
    end

    # Return the position
    return label_position
    
end

function read_xyz_file(path::String)
    
    # Open the file
    io_id = open(path)
    
    # Read the number of atoms
    line = strip(readline(io_id))
    number_of_atoms = parse(Int64, line)
    
    # Ignore the second line
    line = readline(io_id)

    # Define matrices
    atoms_elements = Vector{Int64}(undef, number_of_atoms)
    #atoms_coordinates = Matrix{Float64}(undef, 3, number_of_atoms)
    atoms_coordinates = Matrix{Float64}(undef, 0, 0)

    # Read every atom
    for i in 1:number_of_atoms

        # Read the next line
        line = readline(io_id)
        stripped = strip(line)

        # Get the information
        #symbol, coordx, coordy, coordz  = split(stripped)
        input_line = split(stripped)
        input_dimension = size(input_line, 1)
        symbol = input_line[1]
        coordinates = reshape(parse.(Float64, input_line[2:input_dimension]), input_dimension-1, 1)

        # Store coordinates
        # First line
        if i == 1
            atoms_coordinates = coordinates
        else
            atoms_coordinates = hcat(atoms_coordinates, coordinates)
        end

        # Store element
        position = findall(x -> x == symbol, atomic_information)
        if isempty(position)
            
            println("Element in input file is unknown!")
            println("Element: " * symbol * "; Line: " * string(i+1))
            error("xyz File Error")

        end
        atoms_elements[i] = position[1][1]

    end

    # Close the file
    close(io_id)

    # Return results
    dimension = size(atoms_coordinates,1)
    return number_of_atoms, atoms_elements, atoms_coordinates, dimension

end

# Function to check for any comments, empty lines, or "end" statements
function input_check(textline::AbstractString)
    
    # Skip empty lines
    if strip(textline) == ""
        return 1
    end
    
    # Skip comments indicated by # and !
    if occursin("#", textline)
        return 1
    end
    if occursin("!", textline)
        return 1
    end

    # "end" statement found
    if occursin("End", textline)
        return 2
    end
    if occursin("end", textline)
        return 2
    end

    # Standard return
    return 0

end

# Function to read the general input file
# Read also the xyz file of every molecule
function read_input_file(path::String)
    
    # Allowed keywords
    keyword_blocks = ["Molecule", "molecule", "Lattice", "lattice", "Grid", "grid", "Events", "events"]
    
    molecule_keywords = ["Label", "label", "Rotationmodus", "rotationmodus", "Rotationangle", "rotationangle", "Structure", "structure", "Fixpointtype", "fixpointtype", "Fixpointatoms", "fixpointatoms"]
    molecule_blockkeywords = ["Rotationvalues", "rotationvalues"]
    
    lattice_keywords = ["Transx", "transx", "Transy", "transy"]
    lattice_blockkeywords = ["Vectors", "vectors"]
    
    grid_keywords = ["Label", "label"]
    grid_blockkeywords = ["Points", "points"]    
 
    event_keywords = ["Steps", "steps", "Coverageconvergence", "coverageconvergence", "Forceadsorption", "forceadsorption", "Weigthscale", "weigthscale", "Overlap", "overlap", "Restart", "restart", "Restartruns", "restartruns", "Restartfile", "restartfile"]
    event_blockkeywords = ["Eventlist", "eventlist"]
    event_eventkeywords = ["Ads", "ads", "Dif", "dif", "Rot", "rot", "Con", "con"]

    # Define default vectors to store all structs
    Nmolecules = 0
    molecules = Vector{molecule_struct}(undef, 0)
    Ngrids = 0
    grids = Vector{grid_struct}(undef, 0)
    lattice = lattice_struct()
    events = events_struct()

    # Open input File in "read only"
    io_id = open(path)
    
    # Check every line for keywords
    while ! eof(io_id) 
        
        # Read the next line
        line = readline(io_id)

        # Check input
        checkvalue = input_check(line)
        if checkvalue == 1 || checkvalue == 2 
            continue
        end

        # Check for the main blockkeywords
        stripped = strip(line)
        if stripped in keyword_blocks

            if stripped == "Molecule" || stripped == "molecule"

                Nmolecules, molecules = read_molecule_block!(io_id, Nmolecules, molecules, molecule_blockkeywords, molecule_keywords)

            elseif stripped == "Lattice" || stripped == "lattice"
                
                lattice = read_lattice_block(io_id, lattice_blockkeywords, lattice_keywords)

            elseif stripped == "Grid" || stripped == "grid"

                Ngrids, grids = read_grid_block!(io_id, Ngrids, grids, grid_blockkeywords, grid_keywords)

            elseif stripped == "Events" || stripped == "events"

                events = read_event_block(io_id, event_keywords, event_blockkeywords, event_eventkeywords)

            end

        else
            println("The following block keyword is not known: " * stripped)
            error("Input File Error")
        end
    
    end

    # Close the file again
    close(io_id)

    # Return keywords
    return Nmolecules, molecules, Ngrids, grids, lattice, events

end

# Function to read molecule block
function read_molecule_block!(io_id, Nmolecules, molecules, molecule_blockkeywords, molecule_keywords)
    
    # Get a new molecule struct
    molecule = molecule_struct()
    Nmolecules += 1
    
    # Set "End" counter to 1
    endcounter = 1

    # Read more lines until the closing end is reached
    while endcounter != 0

        # Read the next line
        line = readline(io_id)
        stripped = strip(line)

        # Check input
        checkvalue = input_check(stripped)
        if checkvalue == 1
            continue
        elseif checkvalue == 2
            endcounter -= 1
            if endcounter == 0
                
                # Check whether the molecule has a unique label
                if molecule.label == 0
                    molecule.label = Nmolecules
                end
                if Nmolecules > 1 
                    if any(molecule.label == ele.label for ele in molecules)
                        println("Molecules with identical label found. Label: " * string(molecule.label) * ".")
                        error("Input File Error")
                    end
                end
                
                # Create list of angles for this molecule
                if molecule.rotationmodus == "values"
                    #println("Rotation values explicitly stated for molecule " * string(molecule.label) * ".")
                    molecule.rotations = molecule.rotationvalues
                elseif molecule.rotationmodus == "angle"
                    #println("Rotation angle explicitly stated.")
                    rotationstep = 0.0
                    while rotationstep < 360.0
                        push!(molecule.rotations, rotationstep)
                        rotationstep += molecule.rotationangle
                    end
                else
                    println("No rotation modus found. Specify the 'rotationmodus' keyword.")
                    println("Use either 'values' to specify individual values or 'angle' to set the rotation angle.")
                    println("You wrongly specified: " * molecule.rotationmodus)
                    error("Input File Error")
                end

                # Store number of rotations
                molecule.Nrotations = size(molecule.rotations, 1)

                # Use the obtained path to read the molecule coordinates
                molecule.Natoms, molecule.elements, molecule.coordinates, molecule.dimension = read_xyz_file(molecule.coordinate_path)

                # Define fixpoint of the adsorbate based on the coordinates
                if molecule.fixpointtype == "centroid"
                    molecule.fixpoint = calculate_centroid(molecule.coordinates)
                elseif molecule.fixpointtype == "atoms"
                    molecule.fixpoint = caclulate_partial_centroid(molecule.coordinates, molecule.fixpointatoms)
                else
                    println("The fixpointtype you specified is not known: " * molecule.fixpointtype )
                    println("Use either \"centroid\" or \"atoms\".")
                    error("Input File Error")
                end

                # Add the molecule to the list of molecules
                push!(molecules, molecule)

                # Leave the while loop
                break
            end

            # In case this was not the final end
            continue
        end

        # Check for a blockkeyword
        if stripped in molecule_blockkeywords
            if stripped == "rotationvalues" || stripped == "Rotationvalues"
                # Increase end counter --> Not necessary as the end is catched in the while loop
                # endcounter += 1

                # Every line is a rotation value
                # First line
                line = strip(readline(io_id))
                # Every other line
                while line != "end" && line != "End" && line != "END"
                    push!(molecule.rotationvalues, parse(Float64, line))
                    line = strip(readline(io_id))
                end

                # Reset reading
                continue
            end
        elseif ! occursin("=", stripped) 
            println("The following molecule block keyword is not known: " * stripped)
            error("Input File Error")
        end

        # Check for a keyword
        keyword, value = split(stripped, "=")

        # Check whether keyword is known
        stripped = strip(keyword)
        if stripped in molecule_keywords
            #println("Found keyword: " * stripped)
            if stripped == "Label" || stripped == "label"
                molecule.label = parse(Int64, value)
            elseif stripped == "Rotationmodus" || stripped == "rotationmodus"
                molecule.rotationmodus = strip(value)
            elseif stripped == "Rotationangle" || stripped == "rotationangle"
                molecule.rotationangle = parse(Float64, value)
            elseif stripped == "Structure" || stripped == "structure"
                molecule.coordinate_path = string(strip(value))
            elseif stripped == "Fixpointtype" || stripped == "fixpointtype"
                molecule.fixpointtype = strip(value)
            elseif stripped == "Fixpointatoms" || stripped == "fixpointatoms"
                molecule.fixpointatoms = parse.(Int64, split(value))
            end
        else
            println("The following molecule keyword is not known: " * stripped)
            error("Input File Error")
        end

    end

    # Return results
    return Nmolecules, molecules

end

# Function to read grid block
function read_grid_block!(io_id, Ngrids, grids, grid_blockkeywords, grid_keywords)
    # Get a new grid struct
    grid = grid_struct()
    Ngrids += 1
    
    # Set "End" counter to 1
    endcounter = 1

    # Read more lines until the closing end is reached
    while endcounter != 0

        # Read the next line
        line = readline(io_id)
        stripped = strip(line)

        # Check input
        checkvalue = input_check(stripped)
        if checkvalue == 1
            continue
        elseif checkvalue == 2
            endcounter -= 1
            if endcounter == 0

                # Check whether the grid has a unique label
                if grid.label == 0
                    grid.label = Ngrids
                end
                if Ngrids > 1
                    if grid.label in grids[].label
                        println("Grid with identical label found. Label: " * string(grid.label) * ".")
                        error("Input File Error")
                    end
                end

                # Add to the list of grids
                push!(grids, grid)

                # Leave the while loop
                break
            end

            # In case this was not the final end
            continue
        end

        # Check for a blockkeyword
        if stripped in grid_blockkeywords
            if stripped == "Points" || stripped == "points"
                # Increase end counter --> Not necessary as the end is catched in the while loop
                # endcounter += 1

                # Every line is a grid point
                # First line
                line = readline(io_id)
                string_inputvector = split(line)
                dimension = size(string_inputvector,1)
                grid.uniquepoints = reshape(parse.(Float64, string_inputvector), dimension, 1)
                
                # Every other line
                line = strip(readline(io_id))
                while line != "end" && line != "End" && line != "END"
                    string_inputvector = split(line)
                    grid.uniquepoints = hcat(grid.uniquepoints, parse.(Float64, string_inputvector))
                    line = strip(readline(io_id))
                end

                # Store number of unique grid points
                grid.Nuniquepoints = size(grid.uniquepoints, 2)

                # Store dimension of grid
                grid.dimension = size(grid.uniquepoints, 1)
    
                # Reset reading
                continue
            end
        elseif ! occursin("=", stripped) 
            println("The following grid block keyword is not known: " * stripped)
            error("Input File Error")
        end

        # Check for a keyword
        keyword, value = split(stripped, "=")

        # Check whether keyword is known
        stripped = strip(keyword)
        if stripped in grid_keywords
            #println("Found keyword: " * stripped)
            if stripped == "Label" || stripped == "label"
                grid.label = parse(Int64, value)
            end
        else
            println("The following grid keyword is not known: " * stripped)
            error("Input File Error")
        end

    end

    # Return results
    return Ngrids, grids
end

# Function to read lattice block
function read_lattice_block(io_id, lattice_blockkeywords, lattice_keywords)
    # Get a new lattice struct
    lattice = lattice_struct()
                
    # Set "End" counter to 1
    endcounter = 1

    # Read more lines until the closing end is reached
    while endcounter != 0

        # Read the next line
        line = readline(io_id)
        stripped = strip(line)

        # Check input
        checkvalue = input_check(stripped)
        if checkvalue == 1
            continue
        elseif checkvalue == 2
            endcounter -= 1
            if endcounter == 0
                # Generate the full lattice based on the translations
                lattice.transvectors, lattice.inversevectors = replicate_lattice_with_translation(lattice.vectors, lattice.transx, lattice.transy)

                # Leave the while loop
                break
            end

            # In case this was not the final end
            continue
        end

        # Check for a blockkeyword
        if stripped in lattice_blockkeywords
            if stripped == "Vectors" || stripped == "vectors"
                # Increase end counter
                #endcounter += 1

                # Every line is a lattice vector
                # First line
                line = strip(readline(io_id))
                string_inputvector = split(line)
                lattice.dimension = size(string_inputvector,1)
                lattice.vectors = reshape(parse.(Float64, string_inputvector), 1, lattice.dimension)
                
                # Every other line
                line = strip(readline(io_id))
                while line != "end" && line != "End" && line != "END"
                    string_inputvector = split(line)
                    lattice.vectors = vcat(lattice.vectors, reshape(parse.(Float64, string_inputvector), 1, lattice.dimension))
                    line = strip(readline(io_id))
                end

                # Check that the dimension of the lattice
                if lattice.dimension != size(lattice.vectors, 1)
                    println("The dimension of the lattice and the number of specified coordinates per lattice vector do not match.")
                    println("The defined matrix is not a square matrix.")
                    println("Dimension: " * string(lattice.dimension))
                    println("Coordinates: " * string(size(lattice.vectors, 1)))
                    error("Input File Error")
                end

                # Reset reading
                continue
    
            end
        elseif ! occursin("=", stripped) 
            println("The following lattice block keyword is not known: " * stripped)
            error("Input File Error")
        end

        # Check for a keyword
        keyword, value = split(stripped, "=")

        # Check whether keyword is known
        stripped = strip(keyword)
        if stripped in lattice_keywords
            #println("Found keyword: " * stripped)
            if stripped == "Transx" || stripped == "transx"
                lattice.transx = parse(Int64, value)
                lattice.Ncellx = lattice.transx + 1
            elseif stripped == "Transy" || stripped == "transy"
                lattice.transy = parse(Int64, value)
                lattice.Ncelly = lattice.transy + 1
            end
        else
            println("The following lattice keyword is not known: " * stripped)
            error("Input File Error")
        end

    end

    # Return the lattice input
    return lattice
end

# Function to read events block
function read_event_block(io_id, event_keywords, event_blockkeywords, event_eventkeywords)
    # Get a new events struct
    events = events_struct()

    # Set "End" counter to 1
    endcounter = 1

    # Read more lines until the closing end is reached
    while endcounter != 0

        # Read the next line
        line = readline(io_id)
        stripped = strip(line)

        # Check input
        checkvalue = input_check(stripped)
        if checkvalue == 1
            continue
        elseif checkvalue == 2
            endcounter -= 1
            if endcounter == 0
                # Leave the while loop
                break
            end

            # In case this was not the final end
            continue
        end

        # Check for a blockkeyword
        if stripped in event_blockkeywords
            if stripped == "Eventlist" || stripped == "eventlist"
                # Every line is an event
                line = strip(readline(io_id))
                while line != "end" && line != "End" && line != "END"
                    
                    # Split the input line
                    eventinput = split(line)

                    # Check whether keyword is known
                    if eventinput[2] in event_eventkeywords
                        if eventinput[2] == "Ads" || eventinput[2] == "ads"
                            molecule_label = parse(Int64, eventinput[1])
                            grid_label = parse(Int64, eventinput[3])
                            weigth = parse(Float64, eventinput[4])
                
                            # Create a new event
                            ads_event = event_adsorption_struct(molecule_label, grid_label, weigth)
                
                            # Add the event to the list of events
                            push!(events.adsorptions, ads_event)
                            events.Nevents += 1
                            events.Nadsorptions +=1
                        elseif eventinput[2] == "Dif" || eventinput[2] == "dif"
                            molecule_label = parse(Int64, eventinput[1])
                            grid_1_label = parse(Int64, eventinput[3])
                            grid_2_label = parse(Int64, eventinput[4])
                            weigth = parse(Float64, eventinput[5])
                            diff_radius = parse(Float64, eventinput[6])
                
                            # Create a new event
                            diff_event = event_diffusion_struct(molecule_label, grid_1_label, grid_2_label, weigth, diff_radius)
                
                            # Add the event to the list of events
                            push!(events.diffusions, diff_event)
                            events.Nevents += 1
                            events.Ndiffusions += 1
                        elseif eventinput[2] == "Rot" || eventinput[2] == "rot"
                            molecule_label = parse(Int64, eventinput[1])
                            grid_label = parse(Int64, eventinput[3])
                            weigth = parse(Float64, eventinput[4])
                
                            # Create a new event
                            rot_event = event_rotation_struct(molecule_label, grid_label, weigth)
                
                            # Add the event to the list of events
                            push!(events.rotations, rot_event)
                            events.Nevents += 1
                            events.Nrotations += 1
                        elseif eventinput[2] == "Con" || eventinput[2] == "con"
                            molecule_label_1 = parse(Int64, eventinput[1])
                            molecule_label_2 = parse(Int64, eventinput[3])
                            grid_label = parse(Int64, eventinput[4])
                            weigth = parse(Float64, eventinput[5])
                
                            # Create a new event
                            kon_event = event_conformer_change_struct(molecule_label_1, molecule_label_2, grid_label, weigth)
                
                            # Add the event to the list of events
                            push!(events.conformers, kon_event)
                            events.Nevents += 1
                            events.Nconformers += 1
                        end
                    else
                        println("The following event keyword is not known: " * eventinput[2])
                        error("Input File Error")
                    end
                    
                    # Read the next line
                    line = strip(readline(io_id))
                end              
    
                # Reset reading
                continue
    
            end
        elseif ! occursin("=", stripped) 
            println("The following event block keyword is not known: " * stripped)
            error("Input File Error")
        end

        # Check for a keyword
        keyword, value = split(stripped, "=")

        # Check whether keyword is known
        stripped = strip(keyword)
        if stripped in event_keywords
            #println("Found keyword: " * stripped)
            if stripped == "Steps" || stripped == "steps"
                events.steps = parse(Int64, value)
                if events.steps > 0
                    events.break_steps = true
                end
            elseif stripped == "Coverageconvergence" || stripped == "coverageconvergence"
                events.coverage_convergence = parse(Int64, value)
                if events.coverage_convergence > 0
                    events.break_convergence = true
                end
            elseif stripped == "Forceadsorption" || stripped == "forceadsorption"
                events.Nforce_adsorption = parse(Int64, value)
                if events.Nforce_adsorption > 0
                    events.force_adsorption = true
                end
            elseif stripped == "Weigthscale" || stripped == "weigthscale"
                events.weigth_scale = parse(Float64, value)
                if events.weigth_scale > 0.0
                    events.force_weigth_scale = true
                end
            elseif stripped == "Overlap" || stripped == "overlap"
                if strip(value) == "3D" || strip(value) == "3d"
                    events.overlap2d = false
                    events.overlap3d = true
                elseif strip(value) == "2D" || strip(value) == "2d"
                    events.overlap2d = true
                    events.overlap3d = false
                else
                    println("The following overlap keyword value is not known: " * value)
                    error("Input File Error")
                end
            elseif stripped == "Restart" || stripped == "restart"
                events.restart_generation = parse(Int64, value)
                if events.restart_generation > 0
                    events.restart_flag = true
                end
            elseif stripped == "Restartruns" || stripped == "restartruns"
                events.restart_runs = parse.(Int64, split(value))
            elseif stripped == "Restartfile" || stripped == "restartfile"
                events.restart_file = string(strip(value))
            end
        else
            println("The following event keyword is not known: " * stripped)
            error("Input File Error")
        end

    end

    # Return the event input
    return events
end
