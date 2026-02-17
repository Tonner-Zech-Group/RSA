# Tutorial 8: Mimic Reactions

!!! info "Learning Goals"
    * Combine conversion and diffusion events to mimic reactions 

!!! tip
    Keep in mind that you can speed up simulations by using multiple threads as described in the [Using Multithreading with VSCode](@ref) section.

## Input File Setup 

For this tutorial you can use the coordinates of diethylsulfite and its thiolate from the [Molecule Library](@ref). The largest part of the main input file is identical to previous tutorials: 
```
# Adsorption on Cu(111)
Molecule
    rotationmodus = angle
    rotationangle = 60.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../des.xyz
End

Molecule
    rotationmodus = angle
    rotationangle = 60.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = ...ADJUST-YOUR-PATH.../thiolate.xyz
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

# Bridge grid points of the Cu(111) surface
Grid
    points
        0.00000   2.18137   0.00000
        0.62971   1.09068   0.00000
        0.62971   3.27205   0.00000
        1.25942   0.00000   0.00000
        1.88912   1.09068   0.00000
        1.88912   3.27205   0.00000
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

This time, the eventlist contains several "chained" events: Adsorption events are only possible for molecule 1 - which should be DES. DES can "decompose" and form a thiolate fragment - that should be molecule 2. This conversion takes place only on grid points of the first grid (bridge grid points). The third event contains the diffusion of the thiolate fragment from grid 1 to grid 2 (hollow grid points). Usually a realistic description would always combine the decomposition of DES and the diffsuion of its fragments to their prefered adsorption site. Within our RSA simulation we have to separate these steps and use the *weigths* to make sure that both events usually happen in succession. That is why the weigth for the diffusion step is so large. As a final event the diffusion of molecule 2 on grid 2 is added as we assume the thiolate fragment to be mobile. 
```
# General settings for adsorption of DES  
Events
    steps = 2500
    coverageconvergence = 100
    forceadsorption = 25

    eventlist
	    1 ads 1 1.0
        1 con 2 1 5.0
        2 dif 1 2 1000.0 2.6
        2 dif 2 2 10.0 2.6
    end
End
```

!!! info
    In addition to the current events, you could also add the *backreaction* of the conversion to ensure that it is less likely that a thiolate remains on grid points of grid 1 when neighbouring grid points of grid 2 are currently blocked for the diffusion step. 

!!! info
    Carefully test the impact of the *steps*, *coverageconvergence*, and *forceadsorption* settings on your calculations! Also change the *weight* of the events to test the impact on the simulations. Keep in mind that you only mimic reactivity to reach a limit for surface coverage of a certain adsorbate.

## Running the Simulation
To start the simulations use the [`perform_multiple_rsa_runs`](@ref) function:
```
NRuns = 100
inputfile_path = "...ADJUST-YOUR-PATH.../input.inp"
rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events, timings = perform_multiple_rsa_runs(NRuns, inputfile_path);
```

!!! info
    As this type of simulation is more expensive than simulations of previous tutorials (we use multiple grids, multiple adsorbates, 4 events, and with DES a stretched adsorbate), reduce the number of RSA simulations to 100. 

## Evaluation
As a first step of the evaluation you should check how your input keywords affected the simulations. Here, you can check how many steps were performed by using the *Nsteps* field, how often a certain event was performed by using the *Nevents* field, or which event was still possible at the end of the simulation by using the *stepinfo* field.
- *Nsteps*: Check whether the maximum allowed number of steps or the surface coverage convergence was reached.
- *Nevents*: Check how many adsorption events (firt value) and how many conversion events (fourth number) were performed.
- *stepinfo*: Check the last columns of this matrix. The first value within each column tells you how many events in total were possible while the second and fifth value tell you how many adsorption and conversion events were possible, respectively. What type of events were performed at the end of the simulation?
```
histo, mean, variance, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id = plot_count_area_histograms(NRuns, rsa_results, Nmolecules, molecules, lattice, plotonly = false, area = 0.6);
avgvalue_id[6]

rsa_results[ID of most average run].Nsteps
rsa_results[ID of most average run].Nevents
rsa_results[ID of most average run].stepinfo
```
!!! tip
    More information on the rsa\_results structure can be found in the documentation: [`RSA.rsa_run_results_struct`](@ref).

---

As usual the [`plot_RSA_run`](@ref) function can be used to visualize the final state of a simulation. In addition, it might be interesting to use the [`animate_RSA_run`](@ref) function to create an animation of the RSA simulation: 
```
avgvalue_id[6]
surfaceplot = plot_RSA_run(rsa_results[ID of most average run].status, Ngrids, grids, Nmolecules, molecules, lattice)

anim = animate_RSA_run(rsa_results[ID of most average run].stepinfo, Ngrids, grids, Nmolecules, molecules, lattice)
gif(anim, "...ADJUST-YOUR-PATH.../rsa_animation.gif", fps = 10, loop = -1)
```

!!! warning
    Be aware that the generation of animations is an slow process. It is strongly advised to create animations only with a small (less than 500) number of images. If you want to play around with animations, have a look at the [Plots](https://docs.juliaplots.org/stable/) package, which is used here to create these animations.