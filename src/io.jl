#
# General head section of any file within this module 
#

#
# Include statements
#

#
# Export statements
#
export read_hdf5_output_file

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
 
    event_keywords = ["Steps", "steps", "Coverageconvergence", "coverageconvergence", "Forceadsorption", "forceadsorption", "Weigthscale", "weigthscale", "Overlap", "overlap"]
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
                    if molecule.label in molecules[].label
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
            end
        else
            println("The following event keyword is not known: " * stripped)
            error("Input File Error")
        end

    end

    # Return the event input
    return events
end

# A function to create the hdf5 file with all groups and datasets
function create_hdf5_output_file(inputfile_path)

    # Use the input file path to generate the output file name and path
    file_name = basename(inputfile_path)
    dir_name = dirname(inputfile_path)

    file_prefix = ""
    outputfile_path = ""
    if occursin(".", file_name)
        file_prefix = split(file_name,".")[1]
    else
        file_prefix = file_name
    end
    outputfile_path = joinpath(dir_name, file_prefix * ".h5")

    # Check whether the file already exists
    if isfile(outputfile_path)
        println("The following file already exists: " * outputfile_path)
        println("Please delete this file to start a new simulation.")
        error("Output File Error")
    end

    # Generate the general structure of the hdf5 file
    # Open the file
    hdf5_id = h5open(outputfile_path, "cw")

    # Add attributes to the file
    attributes(hdf5_id)["info"] = "HDF5 output file generated by RSA code"
    attributes(hdf5_id)["url"] = ""
    attributes(hdf5_id)["version"] = "0.1"
    attributes(hdf5_id)["date"] = string(now())
    attributes(hdf5_id)["threads"] = string(Threads.nthreads())
    attributes(hdf5_id)["cpuinfo"] = Base.Sys.cpu_info()[1].model
    attributes(hdf5_id)["runs"] = "1"

    # Generate groups
    Nlevel1 = 4
    level1_groups = ["inputfile-info", "preparation-info", "timings", "run-info"]
    for groups_id in 1:Nlevel1
        create_group(hdf5_id, level1_groups[groups_id])
    end

    # Add attributes to these groups
    group_id = hdf5_id["inputfile-info"]
    attributes(group_id)["info"] = "This group contains all information obtained and generated based on the user input filt. Every molecule and every grid is stored in a separate group."
    group_id = hdf5_id["preparation-info"]
    attributes(group_id)["info"] = "This group contains all distance matrices generated based on the input file. These matrices are used to define possible outcomes of any event. Every run the is stored as a group."
    group_id = hdf5_id["timings"]
    attributes(group_id)["info"] = "This group contains all information regarding execution time of different code routines. For every run the common struct is stored."
    group_id = hdf5_id["run-info"]
    attributes(group_id)["info"] = "This group contains the results of every RSA simulation. Every RSA simulation is stored in a single group."

    # Generate the subgroups for run-1
    for groups_id in 1:Nlevel1
        create_group(hdf5_id, level1_groups[groups_id]*"/run-1")
    end

    # Generate a placeholder for the number of RSA simulations
    group_id = hdf5_id["run-info/run-1"]
    attributes(group_id)["number-of-simulations"] = "0"
    
    # Close the file
    close(hdf5_id)

    # Return the file name
    return outputfile_path

end

# A function to write input file information to an hdf5 file
function write_hdf5_input_information(hdf5_file, Nrun, Nmolecules, molecules, Ngrids, grids, lattice, events)

    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "cw")

    # Store the number of runs in case it is needed later
    # hdf5_id["number-of-runs"] = Nrun

    # Get the group id 
    group_id = hdf5_id["inputfile-info/run-"*string(Nrun)]

    # Write the datasets without subgroups
    group_id["number-of-molecules"] = Nmolecules
    group_id["number-of-grids"] = Ngrids

    # Write the lattice
    write_hdf5_lattice(lattice, group_id)

    # Write the molecules
    write_hdf5_molecules(Nmolecules, molecules, group_id)

    # Write the grids
    write_hdf5_grids(Ngrids, grids, group_id)

    # Write the events
    write_hdf5_events(events, group_id)

    # Close the file
    close(hdf5_id)
    
end

# A function to write lattice information to a hdf5 file
# The file must be open prior to calling this function
function write_hdf5_lattice(lattice, group_id)

    # Generate the group
    create_group(group_id, "lattice")
    lattice_id = group_id["lattice"]

    # Write the information
    lattice_id["dimension"] = lattice.dimension
    lattice_id["vectors"] = lattice.vectors
    lattice_id["transx"] = lattice.transx
    lattice_id["transy"] = lattice.transy
    lattice_id["Ncellx"] = lattice.Ncellx
    lattice_id["Ncelly"] = lattice.Ncelly
    lattice_id["transvectors"] = lattice.transvectors
    lattice_id["inversevectors"] = lattice.inversevectors

end

# A function to read lattice information from a hdf5 file
# The file must be open prior to calling this function
function read_hdf5_lattice(group_id)

    # Select the group
    lattice_id = group_id["lattice"]

    # Generate the struct
    lattice = lattice_struct()

    # Read the information
    lattice.dimension = read(lattice_id, "dimension")
    lattice.vectors = read(lattice_id, "vectors")
    lattice.transx = read(lattice_id, "transx")
    lattice.transy = read(lattice_id, "transy")
    lattice.Ncellx = read(lattice_id, "Ncellx") 
    lattice.Ncelly = read(lattice_id, "Ncelly")
    lattice.transvectors = read(lattice_id, "transvectors")
    lattice.inversevectors = read(lattice_id, "inversevectors")

    # Return results
    return lattice
    
end

# A function to write molecule information to a hdf5 file
# The file must be open prior to calling this function
function write_hdf5_molecules(Nmolecules, molecules, group_id)
    
    # Loop over all molecules
    for molecule_id in 1:Nmolecules

        # Generate the group
        create_group(group_id, "molecule-"*string(molecule_id))
        molecule_group_id = group_id["molecule-"*string(molecule_id)]

        # Write the information
        molecule_group_id["label"] = molecules[molecule_id].label
        molecule_group_id["grids"] = molecules[molecule_id].grids
        molecule_group_id["fixpointtype"] = molecules[molecule_id].fixpointtype
        molecule_group_id["fixpointatoms"] = molecules[molecule_id].fixpointatoms
        molecule_group_id["fixpoint"] = molecules[molecule_id].fixpoint
        molecule_group_id["rotationmodus"] = molecules[molecule_id].rotationmodus
        molecule_group_id["rotationangle"] = molecules[molecule_id].rotationangle
        molecule_group_id["rotationvalues"] = molecules[molecule_id].rotationvalues
        molecule_group_id["Nrotations"] = molecules[molecule_id].Nrotations
        molecule_group_id["rotations"] = molecules[molecule_id].rotations
        molecule_group_id["coordinate_path"] = molecules[molecule_id].coordinate_path
        molecule_group_id["Natoms"] = molecules[molecule_id].Natoms
        molecule_group_id["elements"] = molecules[molecule_id].elements
        molecule_group_id["elements_sorted"] = molecules[molecule_id].elements_sorted
        molecule_group_id["dimension"] = molecules[molecule_id].dimension
        molecule_group_id["coordinates"] = molecules[molecule_id].coordinates
        molecule_group_id["coordinates_sorted"] = molecules[molecule_id].coordinates_sorted
        molecule_group_id["maxradius"] = molecules[molecule_id].maxradius
        molecule_group_id["radii"] = molecules[molecule_id].radii
        molecule_group_id["radii_sorted"] = molecules[molecule_id].radii_sorted

        # Write information needing their own group
        create_group(molecule_group_id, "coordinates_rotated")
        rotations_group_id = molecule_group_id["coordinates_rotated"]
        for rotation_id in 1:molecules[molecule_id].Nrotations
            rotations_group_id["rotation-"*string(rotation_id)] = molecules[molecule_id].coordinates_rotated[rotation_id]
        end

    end

end

# A function to read molecule information from a hdf5 file
# The file must be open prior to calling this function
function read_hdf5_molecules(Nmolecules, group_id)
    
    # Generate the struct
    molecules = Vector{molecule_struct}(undef, 0)

    # Loop over all molecules
    for molecule_id in 1:Nmolecules

        # Select the group
        molecule_group_id = group_id["molecule-"*string(molecule_id)]

        # Generate the struct
        molecule = molecule_struct()

        # Read the information
        molecule.label = read(molecule_group_id, "label")
        molecule.grids = read(molecule_group_id, "grids")
        molecule.fixpointtype = read(molecule_group_id, "fixpointtype")
        molecule.fixpointatoms = read(molecule_group_id, "fixpointatoms")
        molecule.fixpoint = read(molecule_group_id, "fixpoint")
        molecule.rotationmodus = read(molecule_group_id, "rotationmodus")
        molecule.rotationangle = read(molecule_group_id, "rotationangle")
        molecule.rotationvalues = read(molecule_group_id, "rotationvalues")
        molecule.Nrotations = read(molecule_group_id, "Nrotations")
        molecule.rotations = read(molecule_group_id, "rotations")
        molecule.coordinate_path = read(molecule_group_id, "coordinate_path")
        molecule.Natoms = read(molecule_group_id, "Natoms")
        molecule.elements = read(molecule_group_id, "elements")
        molecule.elements_sorted = read(molecule_group_id, "elements_sorted")
        molecule.dimension = read(molecule_group_id, "dimension")
        molecule.coordinates = read(molecule_group_id, "coordinates")
        molecule.coordinates_sorted = read(molecule_group_id, "coordinates_sorted")
        molecule.maxradius = read(molecule_group_id, "maxradius")
        molecule.radii = read(molecule_group_id, "radii")
        molecule.radii_sorted = read(molecule_group_id, "radii_sorted")

        # Read information needing their own group
        rotations_group_id = molecule_group_id["coordinates_rotated"]
        molecule.coordinates_rotated = Vector{Matrix{Float64}}(undef, molecule.Nrotations)
        for rotation_id in 1:molecule.Nrotations
            molecule.coordinates_rotated[rotation_id] = read(rotations_group_id, "rotation-"*string(rotation_id))
        end

        # Add the molecule struct to the vector of structs
        push!(molecules, molecule)

    end

    # Return the result
    return molecules

end

# A function to write grid information to a hdf5 file
# The file must be open prior to calling this function
function write_hdf5_grids(Ngrids, grids, group_id)
    
    # Loop over all grids
    for grid_id in 1:Ngrids

        # Generate the group
        create_group(group_id, "grid-"*string(grid_id))
        grid_group_id = group_id["grid-"*string(grid_id)]

        # Write the information
        grid_group_id["label"] = grids[grid_id].label
        grid_group_id["Nuniquepoints"] = grids[grid_id].Nuniquepoints
        grid_group_id["uniquepoints"] = grids[grid_id].uniquepoints
        grid_group_id["dimension"] = grids[grid_id].dimension
        grid_group_id["Npoints"] = grids[grid_id].Npoints
        grid_group_id["points"] = grids[grid_id].points
        grid_group_id["mapping"] = grids[grid_id].mapping

    end

end

# A function to read grid information from a hdf5 file
# The file must be open prior to calling this function
function read_hdf5_grids(Ngrids, group_id)
    
    # Generate the struct
    grids = Vector{grid_struct}(undef,0)

    # Loop over all grids
    for grid_id in 1:Ngrids

        # Generate the struct
        grid = grid_struct()

        # Select the group
        grid_group_id = group_id["grid-"*string(grid_id)]

        # Read the information
        grid.label = read(grid_group_id, "label")
        grid.Nuniquepoints = read(grid_group_id, "Nuniquepoints")
        grid.uniquepoints = read(grid_group_id, "uniquepoints")
        grid.dimension = read(grid_group_id, "dimension")
        grid.Npoints = read(grid_group_id, "Npoints")
        grid.points = read(grid_group_id, "points")
        grid.mapping = read(grid_group_id, "mapping")

        # Add the grid to the struct vector
        push!(grids, grid)

    end

    # Return the results
    return grids

end

# A function to write event information to a hdf5 file
# The file must be open prior to calling this function
function write_hdf5_events(events, group_id)

    # Generate the group
    create_group(group_id, "events")
    events_id = group_id["events"]

    # Store datasets
    events_id["steps"] = events.steps
    events_id["break_steps"] = events.break_steps
    events_id["coverage_convergence"] = events.coverage_convergence
    events_id["break_convergence"] = events.break_convergence
    events_id["Nforce_adsorption"] = events.Nforce_adsorption
    events_id["force_adsorption"] = events.force_adsorption
    events_id["weigth_scale"] = events.weigth_scale
    events_id["force_weigth_scale"] = events.force_weigth_scale
    events_id["overlap2d"] = events.overlap2d
    events_id["overlap3d"] = events.overlap3d
    events_id["Nevents"] = events.Nevents
    events_id["Nadsorptions"] = events.Nadsorptions
    events_id["Ndiffusions"] = events.Ndiffusions
    events_id["Nrotations"] = events.Nrotations
    events_id["Nconformers"] = events.Nconformers

    create_group(events_id, "adsorptions")
    sub_events_id = events_id["adsorptions"]
    for ads_id in 1:events.Nadsorptions    
        create_group(sub_events_id, "adsorption-"*string(ads_id))
        adsorption_events_id = sub_events_id["adsorption-"*string(ads_id)]
        adsorption_events_id["molecule"] = events.adsorptions[ads_id].molecule
        adsorption_events_id["grid"] = events.adsorptions[ads_id].grid
        adsorption_events_id["weigth"] = events.adsorptions[ads_id].weigth
    end

    create_group(events_id, "rotations")
    sub_events_id = events_id["rotations"]
    for rot_id in 1:events.Nrotations
        create_group(sub_events_id, "rotation-"*string(rot_id))
        rotation_events_id = sub_events_id["rotation-"*string(rot_id)]
        rotation_events_id["molecule"] = events.rotations[rot_id].molecule
        rotation_events_id["grid"] = events.rotations[rot_id].grid
        rotation_events_id["weigth"] = events.rotations[rot_id].weigth
    end

    create_group(events_id, "diffusions")
    sub_events_id = events_id["diffusions"]
    for dif_id in 1:events.Ndiffusions
        create_group(sub_events_id, "diffusion-"*string(dif_id))
        diffusion_events_id = sub_events_id["diffusion-"*string(dif_id)]
        diffusion_events_id["molecule"] = events.diffusions[dif_id].molecule
        diffusion_events_id["grid_start"] = events.diffusions[dif_id].grid_start
        diffusion_events_id["grid_end"] = events.diffusions[dif_id].grid_end
        diffusion_events_id["weigth"] = events.diffusions[dif_id].weigth
        diffusion_events_id["radius"] = events.diffusions[dif_id].radius
    end

    create_group(events_id, "conformers")
    sub_events_id = events_id["conformers"]
    for con_id in 1:events.Nconformers
        create_group(sub_events_id, "conformer-"*string(con_id))
        conformer_events_id = sub_events_id["conformer-"*string(con_id)]
        conformer_events_id["molecule_start"] = events.conformers[con_id].molecule_start
        conformer_events_id["molecule_end"] = events.conformers[con_id].molecule_end
        conformer_events_id["grid"] = events.conformers[con_id].grid
        conformer_events_id["weigth"] = events.conformers[con_id].weigth
    end
    
end

# A function to read event information from a hdf5 file
# The file must be open prior to calling this function
function read_hdf5_events(group_id)

    # Select the group
    events_id = group_id["events"]

    # Generate the events struct
    events = events_struct()

    # Read datasets
    events.steps = read(events_id, "steps")
    events.break_steps = read(events_id, "break_steps")
    events.coverage_convergence = read(events_id, "coverage_convergence")
    events.break_convergence = read(events_id, "break_convergence")
    events.Nforce_adsorption = read(events_id, "Nforce_adsorption")
    events.force_adsorption = read(events_id, "force_adsorption")
    events.weigth_scale = read(events_id, "weigth_scale")
    events.force_weigth_scale = read(events_id, "force_weigth_scale")
    events.overlap2d = read(events_id, "overlap2d")
    events.overlap3d = read(events_id, "overlap3d")
    events.Nevents = read(events_id, "Nevents")
    events.Nadsorptions = read(events_id, "Nadsorptions")
    events.Ndiffusions = read(events_id, "Ndiffusions")
    events.Nrotations = read(events_id, "Nrotations")
    events.Nconformers = read(events_id, "Nconformers")

    # Read subgroups with their datasets
    sub_events_id = events_id["adsorptions"]
    for ads_id in 1:events.Nadsorptions    
        adsorption_events_id = sub_events_id["adsorption-"*string(ads_id)]
        subevent = event_adsorption_struct()
        subevent.molecule = read(adsorption_events_id, "molecule")
        subevent.grid = read(adsorption_events_id, "grid")
        subevent.weigth = read(adsorption_events_id, "weigth")
        push!(events.adsorptions, subevent)
    end

    sub_events_id = events_id["rotations"]
    for rot_id in 1:events.Nrotations
        rotation_events_id = sub_events_id["rotation-"*string(rot_id)]
        subevent = event_rotation_struct()
        subevent.molecule = read(rotation_events_id, "molecule")
        subevent.grid = read(rotation_events_id, "grid")
        subevent.weigth = read(rotation_events_id, "weigth")
        push!(events.rotations, subevent)
    end

    sub_events_id = events_id["diffusions"]
    for dif_id in 1:events.Ndiffusions
        diffusion_events_id = sub_events_id["diffusion-"*string(dif_id)]
        subevent = event_diffusion_struct()
        subevent.molecule = read(diffusion_events_id, "molecule")
        subevent.grid_start = read(diffusion_events_id, "grid_start")
        subevent.grid_end = read(diffusion_events_id, "grid_end")
        subevent.weigth = read(diffusion_events_id, "weigth")
        subevent.radius = read(diffusion_events_id, "radius")
        push!(events.diffusions, subevent)
    end

    sub_events_id = events_id["conformers"]
    for con_id in 1:events.Nconformers
        conformer_events_id = sub_events_id["conformer-"*string(con_id)]
        subevent = event_conformer_change_struct()
        subevent.molecule_start = read(conformer_events_id, "molecule_start")
        subevent.molecule_end = read(conformer_events_id, "molecule_end")
        subevent.grid = read(conformer_events_id, "grid")
        subevent.weigth = read(conformer_events_id, "weigth")
        push!(events.conformers, subevent)
    end
    
    # Return the results
    return events

end

# A function to read an hdf5 file to get all input information
function read_hdf5_input_information(hdf5_file, Nrun)
    
    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Input Reading Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "r")

    # Get the group id 
    group_id = hdf5_id["inputfile-info/run-"*string(Nrun)]

    # Read the datasets without subgroups
    Nmolecules = read(group_id, "number-of-molecules") 
    Ngrids = read(group_id, "number-of-grids")

    # Read the lattice
    lattice = read_hdf5_lattice(group_id)

    # Read the molecules
    molecules = read_hdf5_molecules(Nmolecules, group_id)

    # Read the grids
    grids = read_hdf5_grids(Ngrids, group_id)

    # Read the events
    events = read_hdf5_events(group_id)

    # Close the file
    close(hdf5_id)

    # Return input information
    return Nmolecules, molecules, Ngrids, grids, lattice, events

end

# A function to write preparation information to an hdf5 file
function write_preparation_information(hdf5_file, Nrun, Ngrids, grids, lattice, Nmolecules, molecules, unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, rate_constants_info, neighbour_list)
    
    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "cw")

    # Get the group id 
    main_id = hdf5_id["preparation-info"]["run-"*string(Nrun)]

    # Write gridpoint difference vectors
    tmp_matrix = Matrix{Float64}(undef, 3, 0)
    for grid_B_id in 1:Ngrids
        for grid_A_id in 1:Ngrids

            # Reduce all submatrices to one large matrix
            tmp_matrix_2 = unit_cell_gridpoints_difference[grid_A_id, grid_B_id]
            reduced_matrix = reduce(hcat, tmp_matrix_2)
            tmp_matrix = hcat(tmp_matrix, reduced_matrix)

        end
    end
    main_id["gridpoints-difference"] = tmp_matrix

    # Write the translation distance vectors
    # Reduce to one matrix
    reduced_matrix = reduce(hcat, translation_distance_vectors)
    main_id["translation-vectors"] = reduced_matrix

    # Write the rotation difference vectors
    tmp_matrix = Matrix{Float64}(undef, 3, 0)
    # Loop over molecules
    for molecule_B_id in 1:Nmolecules
        for molecule_A_id in 1:Nmolecules

            # Reduce to one matrix
            tmp_matrix_2 = rotation_difference_matrices[molecule_A_id, molecule_B_id]
            reduced_matrix1 = reduce(hcat, tmp_matrix_2)
            reduced_matrix2 = reduce(hcat, reduced_matrix1)
            tmp_matrix = hcat(tmp_matrix, reduced_matrix2)


        end
    end
    main_id["rotation-differences"] = tmp_matrix

    # Write the affected gridpoints matrices
    create_group(main_id, "affected-points")
    group_id = main_id["affected-points"]
    levels = [false, false, false]
    # Loop over molecules
    molec_subgroup_id = ""
    grid_subgroup_id = ""
    point_subgroup_id = ""
    for molecule_id in 1:Nmolecules

        # Set level 1 to false
        levels[1] = false

        # Loop over gridtypes
        for grid_id in 1:Ngrids

            # Set level 2 to false
            levels[2] = false

            # Loop over gridpoints
            for point_id in 1:grids[grid_id].Nuniquepoints

                # Set level 3 to false
                levels[3] = false

                # Loop over rotations
                for rotation_id in 1:molecules[molecule_id].Nrotations

                    if levels[3] == false

                        # Check whether element is empty
                        if ! isempty(Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id])

                            if levels[1] == false
                                levels[1] = true
                                # Create subgroup
                                molec_subgroup_name = "molecule-vector-" * string(molecule_id)
                                create_group(group_id, molec_subgroup_name)
                                molec_subgroup_id = group_id[molec_subgroup_name]
                            end
                            if levels[2] == false
                                levels[2] = true
                                # Create subgroup
                                grid_subgroup_name = "grid-vector-" * string(grid_id)
                                create_group(molec_subgroup_id, grid_subgroup_name)
                                grid_subgroup_id = molec_subgroup_id[grid_subgroup_name]
                            end
                            if levels[3] == false
                                levels[3] = true
                                # Create subgroup
                                point_subgroup_name = "point-vector-" * string(point_id)
                                create_group(grid_subgroup_id, point_subgroup_name)
                                point_subgroup_id = grid_subgroup_id[point_subgroup_name]
                            end

                            # Create dataset
                            dataset_name = "rotation-vector-" * string(rotation_id)
                            point_subgroup_id[dataset_name] = Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id]
                            continue

                        end

                    end

                    if ! isempty(Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id])
                        
                        # Create dataset
                        dataset_name = "rotation-vector-" * string(rotation_id)
                        point_subgroup_id[dataset_name] = Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id]
                    
                    end

                end

            end

        end

    end

    # Write the neighbour list
    create_group(main_id, "neighbour-list")
    group_id = main_id["neighbour-list"]
    # Loop over starting grid
    levels = [false, false, false]
    sgrid_subgroup_id = ""
    spoint_subgroup_id =""
    molec_subgroup_id = ""
    for start_grid_id in 1:Ngrids

        # Set level 1 to false
        levels[1] = false

        # Loop over starting point
        for start_point_id in 1:grids[start_grid_id].Nuniquepoints

            # Set level 2 to false
            levels[2] = false

            # Loop over molecules
            for molecule_id in 1:Nmolecules

                # Set level 3 to false
                levels[3] = false

                # Loop over end grid
                for end_grid_id in 1:Ngrids

                    if levels[3] == false

                        # Check whether element is empty
                        if ! isempty(neighbour_list[start_grid_id][start_point_id][molecule_id][end_grid_id])

                            if levels[1] == false
                                levels[1] = true
                                # Create subgroup
                                sgrid_subgroup_name = "start-grid-" * string(start_grid_id)
                                create_group(group_id, sgrid_subgroup_name)
                                sgrid_subgroup_id = group_id[sgrid_subgroup_name]
                            end
                            if levels[2] == false
                                levels[2] = true
                                # Create subgroup
                                spoint_subgroup_name = "start-point-" * string(start_point_id)
                                create_group(sgrid_subgroup_id, spoint_subgroup_name)
                                spoint_subgroup_id = sgrid_subgroup_id[spoint_subgroup_name]
                            end
                            if levels[3] == false
                                levels[3] = true
                                # Create subgroup
                                molec_subgroup_name = "molecule-" * string(molecule_id)
                                create_group(spoint_subgroup_id, molec_subgroup_name)
                                molec_subgroup_id = spoint_subgroup_id[molec_subgroup_name]
                            end

                            # Create dataset
                            dataset_name = "end-grid-" * string(end_grid_id)
                            molec_subgroup_id[dataset_name] = neighbour_list[start_grid_id][start_point_id][molecule_id][end_grid_id]
                            continue

                        end

                    end

                    if ! isempty(neighbour_list[start_grid_id][start_point_id][molecule_id][end_grid_id])

                        # Create dataset
                        dataset_name = "end-grid-" * string(end_grid_id)
                        molec_subgroup_id[dataset_name] = neighbour_list[start_grid_id][start_point_id][molecule_id][end_grid_id]

                    end

                end

            end

        end

    end

    # Write rate constants subgroup
    create_group(main_id, "rate-constants")
    subgroup_id = main_id["rate-constants"]
    subgroup_id["ads"] = rate_constants_info.ads
    subgroup_id["rot"] = rate_constants_info.rot
    subgroup_id["dif"] = rate_constants_info.dif
    subgroup_id["dif_rad"] = rate_constants_info.dif_rad
    subgroup_id["con"] = rate_constants_info.con
    
    # Close the file
    close(hdf5_id)

end

# A function to read preparation information from an hdf5 file
function read_preparation_information(hdf5_file, Nrun, Ngrids, grids, lattice, Nmolecules, molecules)
    
    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "r")

    # Get the group id 
    main_id = hdf5_id["preparation-info/run-"*string(Nrun)]

    # Read gridpoint difference vectors
    tmp_read_matrix = read(main_id, "gridpoints-difference")
    unit_cell_gridpoints_difference = Matrix{Matrix{Vector{Float64}}}(undef, Ngrids, Ngrids)
    id_count = 0
    # Loop over gridpoints
    for grid_B_id in 1:Ngrids
        for grid_A_id in 1:Ngrids
            
            # Generate tmp matrix
            gridpoints_difference_matrix = Matrix{Vector{Float64}}(undef, grids[grid_A_id].Nuniquepoints, grids[grid_B_id].Nuniquepoints)

            # Read in the values
            for column_id in 1:grids[grid_B_id].Nuniquepoints
                for row_id in 1:grids[grid_A_id].Nuniquepoints
                    id_count += 1
                    gridpoints_difference_matrix[row_id, column_id] = tmp_read_matrix[:,id_count]
                end
            end

            # Store in the final matrix
            unit_cell_gridpoints_difference[grid_A_id, grid_B_id] = gridpoints_difference_matrix

        end
    end

    # Read translation vectors
    translation_distance_vectors = Matrix{Vector{Float64}}(undef, lattice.transx + 1, lattice.transy + 1)
    tmp_read_matrix = read(main_id, "translation-vectors") 
    # Read in the values
    id_count = 0
    for column_id in 1:lattice.transy + 1
        for row_id in 1:lattice.transx + 1
            id_count += 1
            translation_distance_vectors[row_id, column_id] = tmp_read_matrix[:,id_count]
        end
    end

    # Read rotation difference vectors
    # Create the initial matrix to store all results
    rotation_difference_matrices = Matrix{Matrix{Matrix{Vector{Float64}}}}(undef, Nmolecules, Nmolecules)
    id_count = 0
    tmp_read_matrix = read(main_id, "rotation-differences") 
    # First loop over the molecules
    for molecule_B_id in 1:Nmolecules
        for molecule_A_id in 1:Nmolecules

            # Generate tmp matrix
            tmp_molecule_difference_matrix = Matrix{Matrix{Vector{Float64}}}(undef, molecules[molecule_A_id].Nrotations, molecules[molecule_B_id].Nrotations)

            # Loop over rotations
            for rotation_B_id in 1:molecules[molecule_B_id].Nrotations
                for rotation_A_id in 1:molecules[molecule_A_id].Nrotations
                        
                    # Define blank matrix
                    tmp_matrix = Matrix{Vector{Float64}}(undef, molecules[molecule_A_id].Natoms, molecules[molecule_B_id].Natoms)
            
                    # Read its values
                    for atom_B_id in 1:molecules[molecule_B_id].Natoms
                        for atom_A_id in 1:molecules[molecule_A_id].Natoms
                            id_count += 1
                            tmp_matrix[atom_A_id, atom_B_id] = tmp_read_matrix[:,id_count]
                        end
                    end

                    # Store it in the final matrix
                    tmp_molecule_difference_matrix[rotation_A_id, rotation_B_id] = tmp_matrix

                end
            end

            # Store in the final matrix
            rotation_difference_matrices[molecule_A_id, molecule_B_id] = tmp_molecule_difference_matrix

        end
    end

    # Read the affected points
    # Create an empty array to store all values
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

    # Read in the values 
    group_id = main_id["affected-points"]
    # Loop over molecules
    for molecule_id in 1:Nmolecules

        # Check whether this groups exists
        molec_subgroup_name = "molecule-vector-" * string(molecule_id)
        if ! haskey(group_id, molec_subgroup_name)
            continue
        end
        molec_subgroup_id = group_id[molec_subgroup_name]

        # Loop over gridtypes
        for grid_id in 1:Ngrids

            # Check whether this groups exists
            grid_subgroup_name = "grid-vector-" * string(grid_id)
            if ! haskey(molec_subgroup_id, grid_subgroup_name)
                continue
            end
            grid_subgroup_id = molec_subgroup_id[grid_subgroup_name]

            # Loop over gridpoints
            for point_id in 1:grids[grid_id].Nuniquepoints

                # Check whether this groups exists
                point_subgroup_name = "point-vector-" * string(point_id)
                if ! haskey(grid_subgroup_id, point_subgroup_name)
                    continue
                end
                point_subgroup_id = grid_subgroup_id[point_subgroup_name]

                # Loop over rotations
                for rotation_id in 1:molecules[molecule_id].Nrotations

                    # Check whether this dataset exists
                    dataset_name = "rotation-vector-" * string(rotation_id)
                    if ! haskey(point_subgroup_id, dataset_name)
                        continue
                    end
                    Affected_Points_Rotations[molecule_id][grid_id][point_id][rotation_id] = read(point_subgroup_id, dataset_name)

                end

            end

        end

    end

    # Read the neighbour list
    group_id = main_id["neighbour-list"]
    neighbour_list = Vector{Vector{Vector{Vector{Vector{Int64}}}}}(undef, Ngrids)
    for start_grid_id in 1:Ngrids
        neighbour_list[start_grid_id] = Vector{Vector{Vector{Vector{Int64}}}}(undef, grids[start_grid_id].Nuniquepoints)
        for start_point_id in 1:grids[start_grid_id].Nuniquepoints
            neighbour_list[start_grid_id][start_point_id] = Vector{Vector{Vector{Int64}}}(undef, Nmolecules)
            for molecule_id in 1:Nmolecules
                neighbour_list[start_grid_id][start_point_id][molecule_id] = Vector{Vector{Int64}}(undef, Ngrids)
                for end_grid_id in 1:Ngrids
                    neighbour_list[start_grid_id][start_point_id][molecule_id][end_grid_id] = Vector{Int64}(undef, 0)
                end
            end
        end
    end
    # Loop over starting grid
    for start_grid_id in 1:Ngrids

        # Check whether this groups exists
        sgrid_subgroup_name = "start-grid-" * string(start_grid_id)
        if ! haskey(group_id, sgrid_subgroup_name)
            continue
        end
        sgrid_subgroup_id = group_id[sgrid_subgroup_name]

        # Loop over starting point
        for start_point_id in 1:grids[start_grid_id].Nuniquepoints

            # Check whether this groups exists
            spoint_subgroup_name = "start-point-" * string(start_point_id)
            if ! haskey(sgrid_subgroup_id, spoint_subgroup_name)
                continue
            end
            spoint_subgroup_id = sgrid_subgroup_id[spoint_subgroup_name]

            # Loop over molecules
            for molecule_id in 1:Nmolecules

                # Check whether this groups exists
                molec_subgroup_name = "molecule-" * string(molecule_id)
                if ! haskey(spoint_subgroup_id, molec_subgroup_name)
                    continue
                end
                molec_subgroup_id = spoint_subgroup_id[molec_subgroup_name]

                # Loop over end grid
                for end_grid_id in 1:Ngrids

                    # Check whether this dataset exists
                    dataset_name = "end-grid-" * string(end_grid_id)
                    if ! haskey(molec_subgroup_id, dataset_name)
                        continue
                    end
                    neighbour_list[start_grid_id][start_point_id][molecule_id][end_grid_id] = read(molec_subgroup_id, dataset_name)

                end

            end

        end

    end

    # Read rate constants subgroup
    group_id = main_id["rate-constants"]
    rate_constants_info = rate_constants_struct()
    rate_constants_info.ads = read(group_id, "ads")
    rate_constants_info.rot = read(group_id, "rot")
    rate_constants_info.dif = read(group_id, "dif")
    rate_constants_info.dif_rad = read(group_id, "dif_rad")
    rate_constants_info.con = read(group_id, "con")
    
    # Close the file
    close(hdf5_id)

    # Return the results
    return unit_cell_gridpoints_difference, translation_distance_vectors, rotation_difference_matrices, Affected_Points_Rotations, neighbour_list, rate_constants_info

end

# A function to write timing information
function write_hdf5_timings(hdf5_file, Nrun, timings)
    
    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "cw")

    # Get the group id 
    main_id = hdf5_id["timings"]["run-"*string(Nrun)]

    # Write general timings
    main_id["Processing_Input"] = timings.Processing_Input
    main_id["Calculation_gridpoint_distance"] = timings.Calculation_gridpoint_distance
    main_id["Calculation_translation_distances"] = timings.Calculation_translation_distances
    main_id["Calculation_rotation_distances"] = timings.Calculation_rotation_distances
    #main_id["Calculation_length_sums"] = timings.Calculation_length_sums Currently not used
    main_id["Calculation_affected_points"] = timings.Calculation_affected_points
    main_id["Calculation_overlap"] = timings.Calculation_overlap
    main_id["Calculation_rateconstants"] = timings.Calculation_rateconstants
    main_id["Calculation_neighbour_list"] = timings.Calculation_neighbour_list
    main_id["RSA_runs"] = timings.RSA_runs

    # Currently not used
    #=
    # Combine all timings in one matrix
    rsa_timings_matrix = Matrix{Float64}(undef,4,0)
    for rsa_run_id in eachindex(timings.RSA)
        rsa_timings_matrix = hcat(rsa_timings_matrix, [timings.RSA[rsa_run_id].RSA_update_rateconstants, timings.RSA[rsa_run_id].RSA_select_events, timings.RSA[rsa_run_id].RSA_update_event_list, timings.RSA[rsa_run_id].RSA_compile_information])
    end
    

    # Write the results of every RSA run
    main_id["RSA_runs"] = rsa_timings_matrix
    =#

    # Close the file
    close(hdf5_id)

end

# A funktion to read timing information
function hdf5_read_timings(hdf5_file, Nrun)

    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "r")

    # Get the group id 
    main_id = hdf5_id["timings"]["run-"*string(Nrun)]

    # Create the struct
    timings = timings_struct()

    # Read the general timings
    timings.Processing_Input = read(main_id, "Processing_Input")
    timings.Calculation_gridpoint_distance = read(main_id, "Calculation_gridpoint_distance")
    timings.Calculation_translation_distances = read(main_id, "Calculation_translation_distances")
    timings.Calculation_rotation_distances = read(main_id, "Calculation_rotation_distances")
    # timings.Calculation_length_sums = read(main_id, "Calculation_length_sums") currently not used
    timings.Calculation_affected_points = read(main_id, "Calculation_affected_points")
    timings.Calculation_overlap = read(main_id, "Calculation_overlap")
    timings.Calculation_rateconstants = read(main_id, "Calculation_rateconstants")
    timings.Calculation_neighbour_list = read(main_id, "Calculation_neighbour_list")
    timings.RSA_runs = read(main_id, "RSA_runs")

    # Currently not used
    #=
    # Read the rsa run timings
    rsa_timings_matrix = read(main_id, "RSA_runs")

    # Rearrange the matrix to the struct
    timings.RSA = Vector{rsa_timings_struct}(undef, size(rsa_timings_matrix, 2))
    for rsa_rund_id in 1:size(rsa_timings_matrix, 2)
        timings.RSA[rsa_rund_id] = rsa_timings_struct()
    end
    for colum_id in axes(rsa_timings_matrix, 2)
        timings.RSA[colum_id].RSA_update_rateconstants = rsa_timings_matrix[1, colum_id]
        timings.RSA[colum_id].RSA_select_events = rsa_timings_matrix[2, colum_id]
        timings.RSA[colum_id].RSA_update_event_list = rsa_timings_matrix[3, colum_id]
        timings.RSA[colum_id].RSA_compile_information = rsa_timings_matrix[4, colum_id]
    end
    =#

    # Close the file
    close(hdf5_id)

    # Return the results
    return timings
    
end

# A function to write the results of RSA runs
function write_rsa_results(hdf5_file, Nrun, NRuns, rsa_results)
    
    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "cw")

    # Write the number of performed runs
    main_id = hdf5_id["run-info"]["run-"*string(Nrun)]
    main_id["Nruns"] = NRuns

    # Close the file
    close(hdf5_id)


    # Write every rsa run to the hdf5 file
    for run_id in 1:NRuns
        write_hdf5_rsa_run(hdf5_file, Nrun, run_id, rsa_results[run_id])    
    end
    
end

function write_hdf5_rsa_run(hdf5_file, Nrun, run_id, rsa_result)

    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "cw")

    # Create the group to store all information
    main_id = hdf5_id["run-info"]["run-"*string(Nrun)]
    create_group(main_id, "rsa-run-"*string(run_id))
    group_id = main_id["rsa-run-"*string(run_id)]

    # Generate a tmp matrix to store all elements of the Nevents matrix
    tmp_Nevents = reduce(hcat, rsa_result.Nevents)

    # Write all information
    group_id["randomseed"] = rsa_result.randomseed
    group_id["status"] = rsa_result.status
    group_id["Nsteps"] = rsa_result.Nsteps
    group_id["stepinfo"] = rsa_result.stepinfo
    group_id["Nevents"] = tmp_Nevents
    group_id["stepsize"] = rsa_result.stepsize
    group_id["size"] = rsa_result.size

    # Close the file
    close(hdf5_id)
    
end

# A function to read RSA results information
function read_rsa_results(hdf5_file, Nrun, Nmolecules, Ngrids)
    
    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "r")

    # Read the number of performed runs
    main_id = hdf5_id["run-info"]["run-"*string(Nrun)]
    NRuns = read(main_id, "Nruns") 

    # Close the file
    close(hdf5_id)

    # Allocate a resutls vector
    rsa_results = Vector{rsa_run_results_struct}(undef, NRuns)

    # Read every rsa run from the hdf5 file
    for run_id in 1:NRuns
        rsa_results[run_id] = read_hdf5_rsa_run(hdf5_file, Nrun, run_id, Nmolecules, Ngrids)    
    end

    # Return the results
    return rsa_results

end

function read_hdf5_rsa_run(hdf5_file, Nrun, run_id, Nmolecules, Ngrids)

    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Open the file
    hdf5_id = h5open(hdf5_file, "r")

    # Get the correct group
    group_id = hdf5_id["run-info"]["run-"*string(Nrun)]["rsa-run-"*string(run_id)]

    # Generate the results struct
    rsa_result = rsa_run_results_struct()

    # Read all information
    rsa_result.randomseed = read(group_id, "randomseed")
    rsa_result.status = read(group_id, "status")
    rsa_result.Nsteps = read(group_id, "Nsteps") 
    rsa_result.stepinfo = read(group_id, "stepinfo")
    tmp_Nevents = read(group_id, "Nevents")
    rsa_result.stepsize = read(group_id, "stepsize")
    rsa_result.size = read(group_id, "size")

    # Transform tmp_Nevents
    rsa_result.Nevents = Matrix{Vector{Int64}}(undef, Nmolecules, Ngrids)
    counter = 0
    for molecule_id in 1:Nmolecules
        for grid_id in 1:Ngrids
            counter += 1
            rsa_result.Nevents[molecule_id, grid_id] = tmp_Nevents[1:4, counter]
        end
    end

    # Close the file
    close(hdf5_id)

    # Return results
    return rsa_result
    
end


"""

    read_hdf5_output_file(hdf5_file::String)

Reads and returns all information of RSA simulations from a HDF5 file found under the given path `hdf5_file`.

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
function read_hdf5_output_file(hdf5_file)

    # Check whether the file already exists
    if ! isfile(hdf5_file)
        println("The following HDF5 file does not exist: " * outputfile_path)
        error("Output File Error")
    end

    # Currently no restart features are implemented
    Nrun = 1

    # Read input information
    Nmolecules, molecules, Ngrids, grids, lattice, events = read_hdf5_input_information(hdf5_file, Nrun)

    # Read RSA results
    rsa_results = read_rsa_results(hdf5_file, Nrun, Nmolecules, Ngrids)

    # Read timings
    timings = hdf5_read_timings(hdf5_file, Nrun)

    # Return results
    return rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events, timings
    
end
