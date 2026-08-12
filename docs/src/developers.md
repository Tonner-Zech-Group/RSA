# Developers


## Testsuit
The RSA package is shipped with some tests for the most important parts of the code. Tests can be run via Pkg or by running the `runtest.jl` file within the test folder. It is also possible to run some test individually.
 
* `rsa_tests.jl`

    Running some of the tutorials to ensure that the output of the RSA simulations is correct.

* `io_tests.jl`

    Reading input files (normal and hdf5) to check that the main input information (molecules, grids, lattice, events) are generated correctly.

!!! info
    The test suit is fixing the random seed to ensure that an identical result is obtained in every run. Any changes to how random numbers are generated - especially by Random.jl - might results in all tests failing. In this case tests have to be performed with an older version of Random.jl or manually updated to a newer Version.


## Complete List of Documented Functions and Types
```@autodocs
Modules = [RSA]
Order   = [:function, :type]
```

