# RSA.jl

A *Random Sequential Adsorption (RSA)* package for the modelling of molecular adsorption events within area-selective deposition. 

## Features
* Random adsorption of adsorbates on a surface grid
* Support of multiple adsorbates and surface grids
* Support of diffusion and rotation events
* Support of adsorbate conversion events
* Analysis of surface coverage and effective gap size
* Multithreading support for running multiple simulations simultaneously

All keywords of the input file are discussed in detail in the corresponding sections ([Molecules](@ref molecules), [Lattice](@ref lattice), [Grids](@ref grids), [Events](@ref events)) while an easier introduction into RSA simulations is provided by the [Tutorials](@ref). Furthermore, all functions to run or evaluate RSA simulations are provided in the [Running Simulations](@ref) section.

## Documentation
```@contents
Pages = ["install.md",
         "keywords.md",
         "molecules.md",
         "lattice.md",
         "grids.md",
         "events.md",
         "run.md",
         "tutorials.md",
         "systems/molecules.md",
         "systems/surfaces.md",
         "developers.md"
        ]
Depth = 1
```
