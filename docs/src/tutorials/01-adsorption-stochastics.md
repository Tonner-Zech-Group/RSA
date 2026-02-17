# Tutorial 1: Adsorption of a single adsorbate

!!! info "Learning Goals"
    * Run simple RSA simulations: Aniline on Cu(111)
    * Basic evaluations

!!! tip
    Keep in mind that you can speed up simulations by using multiple threads as described in the [Using Multithreading with VSCode](@ref) section.

## Input File Setup    
For the most simple RSA simulations two input files are needed: The coordinates of the adsorbate and the main input file of the RSA simulation. We will work with aniline as our adsorbate in this tutorial, so please save the following coordinates in an xyz file:
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


Our main input files must contain all information regarding the used [adsorbate](@ref molecules), the [lattice](@ref lattice), the [grid points](@ref grids) as well as possible [events](@ref events). Please combine the following blocks in one file:
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
Within our simulations we will stick aniline by its N atom (first atom in xyz file) to the grid points. In addition, we sample the rotation of aniline around its z-axis in steps of 60°. Keep in mind to adjust the path to your xyz file.

!!! tip
    You can work with absolute paths or use the [cd](https://docs.julialang.org/en/v1/base/file/#Base.Filesystem.cd-Tuple{AbstractString}) command in the notebook to change you working directory and the use relative paths (keep in mind that the relative paths should look like this:"./...RELATIVE-PATH.../aniline.xyz". You can check your current working directory with the [pwd](https://docs.julialang.org/en/v1/base/file/#Base.Filesystem.pwd) command.)

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
As we model the adsorption to a Cu(111) surface, the lattice vectors are those of a Cu(111) unit cell. The z-vector is only needed as we work with 3 dimensions for all our coordinates. As the unit cell is too small to hold many adsorbates, we create a super cell by adding 30 unit cells in x-direction and 18 in y-direction. Consequently, the surface we are modelling is 31x19 unit cells in size or 76Åx79Å.
```
# On-top grid points of the Cu(111) surface
Grid
    points
        0.00000   0.00000   0.00000
        1.25942   2.18137   0.00000
    end
End
```
Based on prior investigatins we assume that aniline only adsorbs to the on-top position of the Cu(111) surface. Therefore, only these adsorption sites are part of our grid points.
```
# General settings for adsorption of aniline  
Events
    steps = 1000

    eventlist
	    1 ads 1 1.0
    end
End
```
To enable the adsorption of aniline the event is added to the eventlist. The "1 ads 1 1.0" line reads as follows: The molecule with the label "1" is "adsorbing" on grid points with the label "1" whereby this event is weigthed to all other events with the factor "1.0". All molecules and grids are independently and automatically labeled based on their order in your input file. As a final setting, we will stop the RSA simulations if no further event is possible or latest when 1000 steps are reached.   

## Running the Simulation
Once both input files are prepared, running the RSA simulations is quite simple: You have to load the RSA module and then start the simulations.
Loading the module is only necessary once for the whole notebook, by adding the path to the source code to the environment:
```
push!(LOAD_PATH,"...ADJUST-YOUR-PATH.../RSA")
using RSA
```

To start the simulations, define a number of RSA simulations you want to perform (here we want to perform 1000 simulations), state the path to the main input file, and call the [`perform_multiple_rsa_runs`](@ref) function. After some time you should get a message that the RSA simulations are done and you are ready for the evaluation steps.
```
NRuns = 1000
inputfile_path = "...ADJUST-YOUR-PATH.../input.inp"
rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events, timings = perform_multiple_rsa_runs(NRuns, inputfile_path);
```

!!! info
    The names ("rsa_results", "Nmolecules", etc) of the variables holding the output of the function call can freely be choosen. The only important part is that 8 variables are provided as the function will return 8 objects.

## Evaluation
For the evaluation we will use the [`plot_count_area_histograms`](@ref), the [`plot_RSA_run`](@ref), and the [`plot_effective_gap_size`](@ref) functions. The first function is used to obtain the histograms showing the distribution of the number of adsorbates and covered area (in %) over all RSA simulations. The output returned by the function is a vector containing adsorbate count and covered area pairs: One pair for every type of adsorbate and one final pair for the sum of all adsorbate types. As we currently use only a single adsorbates, these pairs are identical:
```
histo = plot_count_area_histograms(NRuns, rsa_results, Nmolecules, molecules, lattice, area = 0.6);
histo[3] # identical to histo[1] 
histo[4] # identical to histo[2]
```
As our RSA results are also returned as a vector (one element for one RSA simulation), we can also easily check how the histograms would look in case we take a reduced set of runs. For example we take a subset of 100 runs, namely run 201 to 300: 
```
histo = plot_count_area_histograms(100, rsa_results[201:300], Nmolecules, molecules, lattice, area = 0.6);
histo[3] # identical to histo[1] 
histo[4] # identical to histo[2]
```

!!! info
    In Julia vector elements are always accessed by writing a number or a range (start and end element separated by a colon) in square brackets behind the name of the variable storing the vector. 

---

If you need more information, you can use the optional keyword `plotonly = false` to request more output information:
```
histo, mean, variance, minvalue, minvalue_id, maxvalue, maxvalue_id, avgvalue, avgvalue_id = plot_count_area_histograms(NRuns, rsa_results, Nmolecules, molecules, lattice, plotonly = false, area = 0.6);
histo[4]
mean[4]
variance[4]
```
Now, we obtain the histograms and in addition information about the mean, variance, min and max values, as well as the most average value in terms of adsorbate count and coverage. Check how the mean and variance change, with the number of RSA simulations by evaluating different subsets!

!!! info
    The number of RSA simulations you perform is an important setting. You should test for your system how strongly the mean and variance change with the number of RSA simulations. 

---

You can obtain a visual representation of the covered surface by using the [`plot_RSA_run`](@ref) function. The best representation of all runs is usually the most average RSA simulation:
```
avgvalue_id[4]
surfaceplot = plot_RSA_run(rsa_results[ID of most average run].status, Ngrids, grids, Nmolecules, molecules, lattice)
```
In addition, you can obtain the effective gap sizes - the largest possible circle fitting on any grid point - with the [`plot_effective_gap_size`](@ref) function. Again you can use the most average RSA simulation. As you might be interested in the size range of gaps, it is also advised to have a look at other RSA simulation for example the simulations with the smallest or largest coverage.
```
avgvalue_id[4]
gaphisto, gapplot = plot_effective_gap_size(rsa_results[ID of most average run].status, Ngrids, grids, Nmolecules, molecules, lattice, stepsize = 0.1)
gaphisto
gapplot
```

!!! info
    You might realized that we used the "status" of a RSA simulation in previous function calls. The reason for this is that the results of a RSA simulation are not stored in one large matrix but in multiple numbers, vectors and matrices, which are combined in a so called "struct" (other programming languages call it an "object"). To access any part - in our cases the "status" matrix - of an struct you have to add the name of the subset to the struct name separated by a dot.

