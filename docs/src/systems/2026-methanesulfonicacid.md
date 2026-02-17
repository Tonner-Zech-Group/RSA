# Methanesulfonic Acid (MSA)

- Origin: DFT optimized structure
- Dimension: 3

## Input File Settings
```
# Adsorption on Cu(111)
Molecule
    rotationmodus = angle
    rotationangle = 120.0
    fixpointtype = atoms
    fixpointatoms = 1
    structure = .../structure.xyz
End
```

## Coordinate File, Cu(111)
```
9
xyz file for MSA
S   7.65786   5.99763    9.27688
C   7.89196   5.94626   11.03037
H   6.96905   5.57390   11.48738
H   8.72280   5.25342   11.20492
H   8.13802   6.95457   11.37892
H   6.20477   7.25944    8.29920
O   6.43256   7.01418    9.25446
O   7.25538   4.68725    8.80429
O   8.82600   6.60652    8.65417
```

## Coordinate File, Cu(111), without proton
```
8
xyz file for deprotonated MSA
S   7.55175   5.84738    8.82304
C   7.54062   5.87192   10.59489
H   6.63347   5.35940   10.93177
H   8.43980   5.35291   10.94313
H   7.54197   6.91971   10.91366
O   6.31967   6.55908    8.39407
O   7.55093   4.41761    8.42300
O   8.78986   6.55470    8.40759
```