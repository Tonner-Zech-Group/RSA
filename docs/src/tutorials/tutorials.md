# Tutorials

!!! warning
    All tutorials use simplified input settings to speed up the RSA simulations. For your product runs, carefully benchmark and motivate your settings!

!!! info
    A complete and more extensive explanation of each input keyword is found in the "Input Keywords" part of this documentation starting with the introduction of [Keyword Types](@ref). 

!!! tip
    As common for simulation softwares, it is highly advised to sort your simulations by several folders. If you want to run RSA simulations on a computing cluster, have a look at [Tutorial 9](@ref tutorial-9) to learn how to write all information to a file for later evaluation. 

## Basics
These basic tutorials are designed to be performed in sequence. Each tutorial will build on the knowledge gained in previous tutorials.
```@contents
Pages = [
         "01-adsorption-stochastics.md",
         "02-benchmarks.md",
         "03-multiple-adsorbates.md",
         "04-multiple-grids.md"
        ]
Depth = 1
```

## Advanced
These advanced tutorials are designed to be performed independently. Still, the knowledge gained by the basic tutorials is mandatory.
```@contents
Pages = [
         "05-events-rotations.md",
         "06-events-diffusions.md",
         "07-events-conformer-changes.md",
         "08-events-conversions.md"
        ]
Depth = 1
```

## Situational
```@contents
Pages = [
         "09-hdf5.md"
        ]
Depth = 1
```