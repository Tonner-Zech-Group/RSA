# Tutorial 4: Working with Multiple Grids

!!! info "Learning Goals"
    * Set up simulations with two or more grids
    * Basic evalution of those simulations

!!! tip
    Keep in mind that you can speed up simulations by using multiple threads as described in the [Using Multithreading with VSCode](@ref) section.

## Input File Setup   
The previous tutorial was based on working with multiple adsorbates. However, different adsorbates do not always prefer the same adsorption sites. Therefore, multiple grids are necessasry to offer different adsorption sites. In this tutorial we will work with methanesulfonic acid (MSA) and pyrrole as our adsorbates as these molecules prefer to adsorb at the hollow and on-top site, respectively. Create xyz files with the provided coordinates:
```
8
xyz file for deprotonated MSA
S   7.55175   5.84738    8.82304
C   7.54062   5.87192   10.59489
H   6.63347   5.35940   10.93177
H   8.43980   5.35291   10.94313
H   7.54197   6.91971   10.91366
O   6.31967   6.55908    8.39407
O   7.55093   4.41761    8.42300
O   8.78986   6.55470    8.40759
```
```
10
xyz file for pyrrole
N   7.56794   8.33034   9.01478
C   8.70449   7.55773   8.93506
C   8.30387   6.22604   8.84668
C   6.86940   6.21305   8.85184
C   6.44501   7.53688   8.94245
H   7.55869   9.34894   8.96275
H   9.69526   7.99974   8.99005
H   8.96215   5.36132   8.87377
H   6.22759   5.33612   8.88089
H   5.44698   7.96111   9.00395
```

As we work with two adsorbates, again two molecule blocks are needed in the main input file.
```
# Adsorption on Cu(111)
Molecule
    rotationmodus = angle
    rotationangle = 120.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../msa.xyz
End

Molecule
    rotationmodus = angle
    rotationangle = 60.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../pyrrole.xyz
End
```

We use the lattice of the Cu(111) surface. It is important to note that the lattice must be the same for any grids as only one lattice can be defined in every simulation.
```
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
``` 

To introduce multiple adsorption sites, multiple grid blocks are defined: 
```
# On-top grid points of the Cu(111) surface
Grid
    points
        0.00000   0.00000   0.00000
        1.25942   2.18137   0.00000
    end
End

# Hollow grid points of the Cu(111) surface
Grid
    points
        0.00000   1.45518   0.00000
        0.00000   2.90908   0.00000
        1.25942   0.72771   0.00000
        1.25942   3.63655   0.00000  
    end
End
```

Finally, as MSA (adsorbate 1) is adsorbing on the hollow sites (grid 2) and pyrrole (adsorbate 2) is adsorbing on the on-top sites (grid 1) the following event list is defined:
```
# General settings for adsorption  
Events
    steps = 1000

    eventlist
        1 ads 2 1.0
        2 ads 1 1.0
    end
End
```

## Running the Simulation
To start the simulations use the [`perform_multiple_rsa_runs`](@ref) function. Again, you should get a message that the RSA simulations are done after some time. 
```
NRuns = 1000
inputfile_path = "...ADJUST-YOUR-PATH.../input.inp"
rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events, timings = perform_multiple_rsa_runs(NRuns, inputfile_path);
```

## Evaluation
The [`plot_count_area_histograms`](@ref) and the [`plot_RSA_run`](@ref) functions are our main tools to evaluate the RSA simulations.
```
histo, mean, variance, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id = plot_count_area_histograms(NRuns, rsa_results, Nmolecules, molecules, lattice, plotonly = false, area = 0.6);
```
```
avgvalue_id[6]
surfaceplot = plot_RSA_run(rsa_results[ID of most average run].status, Ngrids, grids, Nmolecules, molecules, lattice)
```

!!! tip
    Change the *weight* of the adsorption events to test the impact on the adsorbate distribution. This setting is impacting the simulations even if adsorbates are not competing for the same adsorption site.
