# [Molecule Keywords](@id molecules)

!!! info
    While all keywords can start with a small or capital letter it is recommended to stick to only small letters for simplicity. 
    All keywords belong to one of the following keyword types: [Block-Keyword](@ref), [Integer-Keyword](@ref integer-keyword), [Float-Keyword](@ref integer-keyword), [Text-Keyword](@ref).

## Overview
The molecule [Block-Keyword](@ref) contains all information for a **single molecule**. In case you are working with several molecules, you have to specify one block per molecule.
```
Molecule [Block-Keyword]

  label [Integer-Keyword]

  rotationmodus [Text-Keyword]
  rotationangle [Float-Keyword]
  rotationvalues [Block-Keyword]
    ... one value per line
  end

  fixpointtype [Text-Keyword]
  fixpointatoms [Text-Keyword]

  structure [Text-Keyword]

End
``` 

## Details
The following list states all keywords of the molecule block with their default value:
* `label = #`

    Gives every molecule a unique label. The default label is based on the occurence order, i.e. the first defined molecule will get the label "1" etc.

* `rotationmodus = angle`

    The rotation modus defines how the rotation angles of the molecule are specified.
    + `angle` : Only a single angle is defined by which the molecule is stepwise rotated.
    + `values` : All values of the rotation angles are explicitly stated. The usage of this keywords demands the usage of the *rotationvalues* keyword.

* `rotationangle = 360.0`

    This value defines the angle for a stepwise rotation of the molecule. A rotation value of 0.0° is always included. The molecule is stepwise rotated by the input value until an angle of 360° would be exceeded. Therefore, a value of 360.0 or larger results in using only the initial orientation. 

* `rotationvalues ... end`

    No default value defined. This keyword must be activated by the *rotationmodus* keyword. A list with one rotation value per input line must be specified by the user.

* `fixpointtyp = centroid`

    The fixpoint of every molecule defines its position on a grid point of the surface model. Two options are available to define a fixpoint:
    + `centroid`: The centroid of the molecule is used to define the position on every grid point.
    + `atoms` : The coordinates of selected atoms are used to define a centroid and thereby the position on grid points. The usage of this keywords demands the usage of the *fixpointatom* keyword.

* `fixpointatom = ...`

    No default value defined. This keyword must be activated by the *fixpointtyp* keyword. All atoms, which should be used to define a centroid, are specified separated by spaces.

* `structure = ...`

    No default value defined. For every molecule a xyz formated file must be provided with the atomic elements and coordinates. The absolut path to this file is stated with this keyword.
    
!!! warning
    Atomic coordinates provided in the xyz file of any molecule, the coordinates of specified grid points, and the stated lattice vectors must adhere to the same dimensionality!

