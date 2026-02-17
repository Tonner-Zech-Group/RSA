# Keyword Types

All keywords of this package belong to one of the following types. The type of keyword defines how a given input must be specified. While some keywords have a default value, it is **highly recommended to not rely** on default settings. Instead, any relevant keyword should be explicitly stated in the simulation.

!!! info
    Every line containing a "!" or "#" sign is considered to be a comment and therefore ignored.

### Block-Keyword
Block-keywords are used to group a list of keywords or input parameters. Every block is started with the keyword - "Molecule" in the following example - and closed by an "End" statement. Every line between the keyword and end statement is considered as a single input. 
```
Block-Keywords:
Molecule
 
 ... more keywords and inputs
 ... one input per line

End
```

### [Integer-Keyword and Float-Keyword](@id integer-keyword)
Integer and float-keywords expect a single number as an input. The keyword and its value are seperated by an equal "=" sign.
```
Integer-Keywords:
steps = 10000

Float-Keywords:
rotationangle = 10.0
```

### Text-Keyword
Non numeric values and lists of values are handled by text-keywords. The keyword and input is again seperated by an equal "=" sign. The input text should not be surrounded by any quotes. In case of a single input, the value should not contain spaces. For lists, values are separated by a space (and not by any special sign). 
```
Text-Keywords:
rotationmodus = angle
fixpointatoms = 1 2 5 8

```

