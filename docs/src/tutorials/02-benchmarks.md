# Tutorial 2: Benchmarking of Input Settings

!!! info "Learning Goals"
    * Benchmark simulation cell size
    * Benchmark impact of rotation steps

!!! tip
    Keep in mind that you can speed up simulations by using multiple threads as described in the [Using Multithreading with VSCode](@ref) section.

## General
Every RSA simulation contains two important choices: The size of the unit cell as well as the number of rotations for every adsorbate. While the cell size can be adjusted based on a requested accuracy, the number of rotations should be motivated by the physics of the modeled system. Examples are small rotation steps for adsorbates that are considered to be freely rotating, a step size of 60° or 120° if an adsorbate is adapting to the symmetry of a Cu(111) surface, or even values of 180° or 360° (0°) if the adsorbate is assumend to be static. The value and number of valid rotations is therefore benchmarked to gain information on the sensitivity of the RSA simulations to this value.

!!! info
    The larger the unic cell and the smaller the rotation steps the more demanding are the RSA simulations. So choose wisely! 

## Input File Setup
In this tutorial we will again use aniline as our adsorbate. Store the following coordinates in an xyz file:
```
14
xyz file for aniline
N   5.04118   8.64882   8.56544
C   6.26128   5.16607   9.00043
C   6.26490   6.56249   8.89207
C   5.04571   7.26376   8.83488
C   3.83123   6.55488   8.89628
C   3.84412   5.15850   9.00450
C   5.05496   4.45477   9.03852
H   7.21439   4.63227   9.04668
H   7.21136   7.11047   8.87828
H   2.88123   7.09678   8.88507
H   2.89468   4.61854   9.05367
H   5.05865   3.36626   9.11942
H   5.88382   9.14949   8.85915
H   4.19684   9.14457   8.86255
```

The structure of the main input file is identical to the previous tutorial. You must create different input files to benchmark the cell size and rotation values. 
```
# Adsorption on Cu(111)
Molecule
    rotationmodus = angle
    rotationangle = 60.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../aniline.xyz
End
```
For a benchmark of the rotation values change the *rotationangle* settings from 360° to 5° or even 2° in case you are patient (recommended: 360°, 180°, 90°, 60°, 45°, 10°, 5°). Keep the cell size fixed during this benchmark. Start with the largest value to get a feeling of the growing compute time with decreasing values.
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
For a benchmark of the cell size change the *transx* and *transy* settings from 50x30 to 10x6 (recommended: 50x30, 40x24, 30x18, 20x12, 10x6). Again, keep other settings fixed. Start with the smallest value to get a feeling of the compute time. All other settings are identical to the previous tutorial:
```
# On-top grid points of the Cu(111) surface
Grid
    points
        0.00000   0.00000   0.00000
        1.25942   2.18137   0.00000
    end
End

# General settings for adsorption of aniline  
Events
    steps = 1000

    eventlist
	    1 ads 1 1.0
    end
End
```

## Running the Simulation
To start the simulations use the [`perform_multiple_rsa_runs`](@ref) function. Again, you should get a message that the RSA simulations are done after some time based an your chosen benchmark settings.
```
NRuns = 1000
inputfile_path = "...ADJUST-YOUR-PATH.../input.inp"
rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events, timings = perform_multiple_rsa_runs(NRuns, inputfile_path);
```

## Evaluation
Use the [`plot_count_area_histograms`](@ref) function to get the mean and variance of your simulations with changing inputs.
```
histo, mean, variance, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id = plot_count_area_histograms(NRuns, rsa_results, Nmolecules, molecules, lattice, plotonly = false, area = 0.6);
mean[4]
variance[4]
```

!!! info
    A reasonable setting should show a converged mean value and a variance small enough to be acceptable for your desired accuracy. Keep in mind that these settings are system dependend!