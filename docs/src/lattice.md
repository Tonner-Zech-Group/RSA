# [Lattice Keywords](@id lattice)

!!! info
    While all keywords can start with a small or capital letter it is recommended to stick to only small letters for simplicity. 
    All keywords belong to one of the following keyword types: [Block-Keyword](@ref), [Integer-Keyword](@ref integer-keyword), [Float-Keyword](@ref integer-keyword), [Text-Keyword](@ref).

## Overview
The lattice [Block-Keyword](@ref) contains all information for the lattice of the surface model. Periodic boundary conditions are always included. Importantly, the final surface cell is defined by a unit cell, which is repeated in x and y-direction. For performance reasons: Choose your unit cell as small as possible and adjust the number of repetitions to obtain your desired surface model. 
```
Lattice

  transx [Integer-Keyword]
  transy [Integer-Keyword]

  vectors [Block-Keyword]
    ... one vector per line	
  end

End
```

## Details
The following list states all keywords of the lattice block with their default value:
* `transx = 0`

    Defines the number of repetitions in x-direction. A value of "0" indicates that only the original unit cell is used.

* `transy = 0`

    Defines the number of repetitions in y-direction. A value of "0" indicates that only the original unit cell is used.

* `vectors ... end`

    This block-keyword is used to define the lattice vectors of the unit cell. Every line specifies one vector. Importantly, the number of vectors has to match the dimension of each vector, i.e. if you work with **3** coordinates per vector you have to specify **3** unit vectors.

!!! tip
    For the best performance of the RSA simulations you should choose a unit cell as small as possible to represent the symmetry of the surface. Here, it is not relevant whether your adsorbate could fit in this unit cell.

!!! warning
    Atomic coordinates provided in the xyz file of any molecule, the coordinates of specified grid points, and the stated lattice vectors must adhere to the same dimensionality!

