# Tutorial 5: Rotation of Adsorbates

!!! info "Learning Goals"
    * Include adsorbate rotations in simulations
    * Adjust input settings for "infinite" simulations
    * Evaluate rotation events 

!!! tip
    Keep in mind that you can speed up simulations by using multiple threads as described in the [Using Multithreading with VSCode](@ref) section.

## Input File Setup  

For this tutorial you can use the coordinates of aniline from the [Molecule Library](@ref). The largest part of the main input file is identical to previous tutorials: 
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
```

The major difference in the input file is the presence of [rotation events](@ref events) for the adsorbate in the eventlist. Furthermore *coverageconvergence* and *forceadsorption* are added. The first keyword will stop the RSA simulation once the given number of steps were performed without the occurence of an adsorption event - thereby the surface coverage is assumed to be converged. The second keyword will try to force an adsorption event once the given the given number of steps were performed without the occurence of an adsorption event. 
```
# General settings for adsorption of aniline  
Events
    steps = 1000
    coverageconvergence = 100
    forceadsorption = 25

    eventlist
        1 ads 1 1.0
        1 rot 1 1.0
    end
End
```
!!! info
    Carefully test the impact of the *steps*, *coverageconvergence*, and *forceadsorption* settings on your calculations! Also change the *weight* of the events to test the impact on the simulations.

## Running the Simulation
To start the simulations use the [`perform_multiple_rsa_runs`](@ref) function:
```
NRuns = 1000
inputfile_path = "...ADJUST-YOUR-PATH.../input.inp"
rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events, timings = perform_multiple_rsa_runs(NRuns, inputfile_path);
```

## Evaluation
As a first step of the evaluation you should check how your input keywords affected the simulations. Here, you can check how many steps were performed by using the *Nsteps* field, how often a certain event was performed by using the *Nevents* field, or which event was still possible at the end of the simulation by using the *stepinfo* field.
- *Nsteps*: Check whether the maximum allowed number of steps or the surface coverage convergence was reached.
- *Nevents*: Check how many adsorption events (firt value) and how many rotation events (second number) were performed.
- *stepinfo*: Check the last columns of this matrix. The first value within each column tells you how many events in total were possible while the second and third value tell you how many adsorption and rotations events were possible, respectively. What type of events were performed at the end of the simulation?
```
histo, mean, variance, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id = plot_count_area_histograms(NRuns, rsa_results, Nmolecules, molecules, lattice, plotonly = false, area = 0.6);
avgvalue_id[4]

rsa_results[ID of most average run].Nsteps
rsa_results[ID of most average run].Nevents
rsa_results[ID of most average run].stepinfo
```
!!! tip
    More information on the rsa\_results structure can be found in the documentation: [`RSA.rsa_run_results_struct`](@ref).

---

As usual the [`plot_RSA_run`](@ref) function can be used to visualize the final state of a simulation. In addition, it might be interesting to use the [`animate_RSA_run`](@ref) function to create an animation of the RSA simulation: 
```
avgvalue_id[4]
surfaceplot = plot_RSA_run(rsa_results[ID of most average run].status, Ngrids, grids, Nmolecules, molecules, lattice)

anim = animate_RSA_run(rsa_results[ID of most average run].stepinfo, Ngrids, grids, Nmolecules, molecules, lattice)
gif(anim, "...ADJUST-YOUR-PATH.../rsa_animation.gif", fps = 10, loop = -1)
```

!!! warning
    Be aware that the generation of animations is an slow process. It is strongly advised to create animations only with a small (less than 500) number of images. If you want to play around with animations, have a look at the [Plots](https://docs.juliaplots.org/stable/) package, which is used here to create these animations.
