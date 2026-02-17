# Tutorial 3: Working with Multiple Adsorbates

!!! info "Learning Goals"
    * Set up simulations with two or more adsorbates
    * Basic evalution of those simulations

!!! tip
    Keep in mind that you can speed up simulations by using multiple threads as described in the [Using Multithreading with VSCode](@ref) section.

## Input File Setup
Some simulations work with several adsorbates *simultanously*. In this case, every adsorbate must be provided by its own xyz file. In this tutorial we will use pyridine in its upright and tilted conformation as an example. Use the following coordinates to create two xyz files:   
```
11
xyz file for pyridine (upright)
N   6.28358   6.56065    8.45355
C   6.25263   7.72424    9.13811
C   6.23112   7.77781   10.52736
C   6.24195   6.58463   11.25023
C   6.27507   5.37951   10.54820
C   6.29524   5.40861    9.15778
H   6.22511   6.59391   12.34275
H   6.24512   8.63611    8.52875
H   6.20608   8.74862   11.02618
H   6.28509   4.41787   11.06492
H   6.32100   4.48683    8.56448
```
```
11
xyz file for pyridine (tilted)
N   6.29357   6.02642    8.42824
C   7.43949   5.53364    8.95227
C   7.45754   4.55838    9.94549
C   6.24499   4.07395   10.43689
C   5.05858   4.62236    9.94926
C   5.12504   5.59566    8.95637
H   6.22516   3.28470   11.19206
H   8.36942   5.93498    8.53437
H   8.41447   4.17565   10.30539
H   4.08392   4.29174   10.31264
H   4.21642   6.04634    8.54133
```

As we want to work with multiple adsorbates, a molecule block is required in the main input file for every adsorbate.
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
```
Other settings regarding the lattice and grid remain unchanged in comparison to previous tutorials.
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

# On-top grid points of the Cu(111) surface
Grid
    points
        0.00000   0.00000   0.00000
        1.25942   2.18137   0.00000
    end
End
```

However, the eventlist must be adapted to enable the adsorption of both adsorbates. 
```
# General settings for adsorption of pyridine  
Events
    steps = 1000

    eventlist
        1 ads 1 1.0
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
With the [`plot_count_area_histograms`](@ref) function you can obtain the information how many molecules of each adsorbate are present. As we work with two adsorbates, every output vector has six elemtens: Element 1-2 for adsorbate 1, element 3-4 for asdorbate 2, and element 5-6 for the sum of both adsorbates.
```
histo, mean, variance, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id = plot_count_area_histograms(NRuns, rsa_results, Nmolecules, molecules, lattice, plotonly = false, area = 0.6);
```
For a visualization you can use the [`plot_RSA_run`](@ref) function.

```
avgvalue_id[6]
surfaceplot = plot_RSA_run(rsa_results[ID of most average run].status, Ngrids, grids, Nmolecules, molecules, lattice)
```

!!! tip
    Change the *weight* of the adsorption events to test the impact on the adsorbate distribution.