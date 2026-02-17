# [Grids Keywords](@id grids)

!!! info
    While all keywords can start with a small or capital letter it is recommended to stick to only small letters for simplicity. 
    All keywords belong to one of the following keyword types: [Block-Keyword](@ref), [Integer-Keyword](@ref integer-keyword), [Float-Keyword](@ref integer-keyword), [Text-Keyword](@ref).

## Overview
The grid [Block-Keyword](@ref) contains all information for a **single grid type** of the surface model. In case you are working with several grids, you have to specify one block per grid.
```
Grid

  label [Integer-Keyword]

  points [Block-Keyword]
    ... one grid point per line
  end

End
```

## Details
The following list states all keywords of the grid block with their default value:
* `label = #`

    Gives every grid a unique label. The default label is based on the occurence order, i.e. the first defined grid will get the label "1" etc.

* `points ... end`

    Within the points block all grid points are defined. Every line states the cartesian coordinates of a single grid point.

!!! warning
    Atomic coordinates provided in the xyz file of any molecule, the coordinates of specified grid points, and the stated lattice vectors must adhere to the same dimensionality!