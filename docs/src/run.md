# Running Simulations

After completing the [Installation](@ref) of the package, you shoud have access to several high level functions to run and analyze RSA simulations. A complete lists of high level functions is shown in the [Exported Functions](@ref) section.

To start RSA simulations only two inputs, namely the number of simulations you want to perform and the path to the input file, are needed:
```
NRuns = 1000
inputfile_path = "/home/fabian/Sync/Coding/Julia/Scripts/lattice.inp"
rsa_results, Nmolecules, molecules, Ngrids, grids, lattice, events, timings = perform_multiple_rsa_runs(NRuns, inputfile_path)
```
As a result you will obtain all information found in the input file (*Nmolecules*, *molecules*, *Ngrids*, *grids*, *lattice*, *events*), the timings of the individual steps (*timings*), as well as all results of the rsa simulations (*rsa_results*). Details on the obtained objects are explained in the [Exported Functions](@ref) section as well as the [Datastructure](@ref) section.


For simple analysis tasks, the obtained objects are simply passed to another function:
```
myplot = plot_RSA_run(rsa_results[352].status, Ngrids, grids, Nmolecules, molecules, lattice)
savefig(myplot, "RSA_run_352.png")
```
Default functions for plotting RSA runs, evaluating the covered area, and calculating the effective gap size are provided in the [Analysis](@ref) section. 

## Using Multithreading with VSCode
Running RSA simulations over VSCode in parallel is possible by changing the "julia.NumThreads" setting. Simply search for "num thread" in the settings searchbar and change to your needs. 
To test the new settings use the following command in your notebook:
```
Threads.nthreads()
```

## Exported Functions

### RSA Runs
```@docs
perform_multiple_rsa_runs
read_hdf5_output_file
```

### Analysis
```@docs
plot_RSA_run
plot_single_molecule
animate_RSA_run
plot_count_area_histograms
plot_effective_gap_size
```

## Datastructure
```@docs; canonical=false
RSA.rsa_run_results_struct
RSA.molecule_struct
RSA.grid_struct
RSA.lattice_struct
RSA.events_struct
RSA.timings_struct
```