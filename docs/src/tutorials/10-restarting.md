# [Tutorial 10: Restarting a Simulation Based on a Precoved Surface](@id tutorial-10)

!!! info "Learning Goals"
    * Use already created HDF5 file to restart RSA simulations
    * Adsorb multiple molecules in sequence on a surface


## Input File Setup
While interactively working with RSA simulations is often a first step, you might want the results to be stored. For this you must use the [HDF5](https://www.hdfgroup.org/solutions/hdf5/) features of this package.
. In this way you enable yourself to evaluate the RSA simulations at a later point and in addition are able to upload your data for any publication. Within this tutorial you can use the coordinates of Pyridine from the [Molecule Library](@ref). The main input file is identical to previous tutorials: 
```
# Adsorption on Cu(111)
Molecule
    rotationmodus = angle
    rotationangle = 60.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../pyridine-upright.xyz
End

Molecule
    rotationmodus = angle
    rotationangle = 60.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../pyridine-tilted.xyz
End

# Lattice of the Cu(111) surface
Lattice
    transx = 30
    transy = 18
    
    vectors
        2.51883   0.00000   0.00000
        0.00000   4.36274   0.00000
        0.00000   0.00000   1.00000
    end
End

# On-top grid points of the Cu(111) surface
Grid
    points
        0.00000   0.00000   0.00000
        1.25942   2.18137   0.00000
    end
End

# General settings for adsorption of pyridine  
Events
    steps = 1000
    #coverageconvergence = 100
    #forceadsorption = 25

    eventlist
	    1 ads 1 1.0
        2 ads 1 1.0
        1 con 2 1 10.0
    end
End
```

## Running the Simulation
To start the simulations you still use the [`perform_multiple_rsa_runs`](@ref) function but this time with the optional hdf5 flag enabled:
```
NRuns = 1000
inputfile_path = "...ADJUST-YOUR-PATH.../input.inp"
rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events = perform_multiple_rsa_runs(NRuns, inputfile_path, hdf5 = true);
```
This time a "h5" file is created storing all information of the RSA simulations - for the naming the "inp" ending of your input file is replaced with "h5".