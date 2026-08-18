# [Tutorial 10: Restarting a Simulation Based on a Precovered Surface](@id tutorial-10)

!!! info "Learning Goals"
    * Use already created HDF5 file to restart RSA simulations
    * Adsorb multiple molecules in sequence on a surface


## Input File Setup
All previous tutorials used a single adsorbate or multiple adsorbates simultaneously. To add multiple adsorbates in sequence a restart run is used. Basic steps are to create a HDF5 file in a first simulation and use a result from that calculation as the starting point for a new calculation. Consequently, at least two input files (for two adsorbates in sequence) are used. Within this tutorial you can use the coordinates of aniline and pyridine (upright structure) from the [Molecule Library](@ref). The first input file is nearly identical to previous tutorials: 
```
# Adsorption on Cu(111)
Molecule
    rotationmodus = angle
    rotationangle = 60.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../aniline.xyz
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

# General settings for adsorption of aniline  
Events
    steps = 30

    eventlist
	    1 ads 1 1.0
    end
End
```
The important change is that the number of steps per simulation is limited to 30. In this way we artificially stop the simulations before the surface is completely covered. Keep in mind that this trick is only meaningful for this tutorial. 


The second input file contains in addition all information for pyridine. In addition three restart keywords are present, which state the path to the HDF5 file, the generation of RSA simulations used, and a list of the individual RSA simulations of this generation used as initial seed.
```
# Adsorption on Cu(111)
Molecule
    rotationmodus = angle
    rotationangle = 60.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../aniline.xyz
End

Molecule
    rotationmodus = angle
    rotationangle = 60.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../pyridine-upright.xyz
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

# General settings for adsorption of aniline & pyrrole 
Events
    steps = 1000
    
    restart = 1
    restartruns = 1 115 541 
    restartfile = ...ADJUST-YOUR-PATH.../input-1.h5

    eventlist
	    2 ads 1 1.0
        1 dif 1 1 10.0 2.6
        1 rot 1 10.0
    end
End
```

!!! warning
    New molecules and grids of a restart run must appear after the original molecules and grids. Furthermore, you must not delete any molecule or grid block used in the initial calculations even if they are no longer used. 

## Running the Simulation
To start the simulations you still use the [`perform_multiple_rsa_runs`](@ref) function. You must use the optional hdf5 flag in both calls:
```
NRuns = 1000
inputfile_path = "...ADJUST-YOUR-PATH.../input-1.inp"
rsa_results-1, Nmolecules-1, molecules-1, Ngrids-1, grids-1, lattice-1, events-1 = perform_multiple_rsa_runs(NRuns, inputfile_path, hdf5 = true);

inputfile_path = "...ADJUST-YOUR-PATH.../input-2.inp"
rsa_results-2, Nmolecules-2, molecules-2, Ngrids-2, grids-2, lattice-2, events-2 = perform_multiple_rsa_runs(NRuns, inputfile_path, hdf5 = true);
```
The "h5" file will contain all information of both runs (as generation 1 and 2) while you can use all common analysis functions on the result structs. Keep in mind that the number of runs specified by `NRuns` will be performed for every selected run (1, 115, and 541) resulting in 3000 RSA simulations in this example. 

!!! info
    The result structs of a restart run also include the information of the selected initial surface. For example, if you adsorb multiple adsorbates in a series of restart runs, every histogram counting the number of adsorbates will inlcude all adsorbates currently on the surface and not only the newly added ones.  