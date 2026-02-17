# [Events Keywords](@id events)

## Overview
The events [Block-Keyword](@ref) contains all information of the events taking place in the RSA simulation. In addition, more general settings as for example the maximum number of steps are defined here. 
```
Events [Block-Keyword]

  steps [Integer-Keyword]
  coverageconvergence [Integer-Keyword]
  forceadsorption [Integer-Keyword]
  overlap [Text-Keyword]

  eventlist [Block-Keyword]
    ... one event per line
  end

End
``` 

## Details
The following list states all keywords of the events block with their default value:
* `steps = 0`

    This keyword defines the maximum number of steps to be performed in the RSA simulations. A value of "0" indicates that there is no maximum number of steps. Either this or the *coverageconvergence* keyword must be set by the user.

* `coverageconvergence = 0`

    This keyword defines the number of consecutive non-adsorption events at which the surface coverage is consider to be converged. Once this threshold is reached the simulation stops. A value of "0" indicates that the convergence of coverage is not used. Either this or the *steps* keyword must be set by the user.

* `forceadsorption = 0`

    To enforce the selection of adsorption events this keyword can be used. It defines the number of non-adsorption events after which an adsorption should be enforced - only if an adsorption event is currently possible. A value of "0" deactivates this behavior.

* `overlap = 2D`

    This keyword defines how the overlap between two molecules is defined:
    + `2D` : The van der Waals spheres of the atoms are projected in two dimensions to judge the overlap between adsorbates. This is the default approach in most RSA simulations.
    + `3D` : The overlap between the van der Waals spheres of the atoms is tested in three dimensions. This setting can only be used in case all coordinates (molecules, lattice, grids) are stated with three dimensions. 

* `eventlist ... end`

    This block keyword is used to define all possible events in the RSA simulation. To define events the labels of the molecules and grids are used. In addition, every event contains a "weight" to balance the likelyhood of certain events to each other. Every line states one event.
    + `[molecule label] ads [grid label] [weigth]` : For an adsorption event only the molecule, the grid and the weight of this event are needed.
    + `[molecule label] rot [grid label] [weigth]` : For an rotation event again only the molecule, the grid and the weight of this event are needed.
    + `[molecule label] dif [grid label] [grid label] [weigth] [diffusion radius]` : Diffusion events are defined by stating the molecule, the starting grid, the destination grid, the weigth, and the maximum diffusion distance (in angstrom) in order.
    + `[molecule label] con [molecule label] [grid label] [weigth]` : Conversion events state the starting molecule, the destination molecule, the grid type - on which the conversion takes place -, and the weight of this event.

    An example for a diffusion event is shown in the following:
    + `2 dif 1 2 10.0 3.2`

