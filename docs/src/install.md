# Installation

!!! info
    To install and mange your version of [Julia](http://www.julialang.org) we recommend the usage of [Juliaup](https://github.com/JuliaLang/juliaup).
    Furthermore, we recommend the usage of [VSCode](https://code.visualstudio.com/docs/languages/julia) and its [Julialang](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia) extension to work interactively with the RSA package, while the [IJulia](https://julialang.github.io/IJulia.jl/dev/manual/installation/) package can be used to run Julia in Jupyter notebooks.

## Get the Code
As our RSA package is not registered yet you have to download the source code and make it available by changing the *LOAD_PATH* variable later on. Clone the repository to get the code:
```
git clone https://github.com/Tonner-Zech-Group/RSA
```

## Install Dependencies
Within a terminal, install all dependencies by switching in the *RSA* root folder (or use the full path in the following commands). Start the Julia commnad line (usually by typing `julia`) and activate the RSA environment and install all dependencies with the following commands:
```
using Pkg
```
```
Pkg.activate("/PATH/to/source/code/RSA")
```
```
Pkg.instantiate()
```

After the installation of all dependencies you can switch back to the default environment:
```
Pkg.activate()
```

## Using the RSA package
You can use the RSA package in any of your scripts without activating the project (this would only be necessary if you want to contribute to the development of the package). Simply add the package directory to the *LOAD_PATH*:
```
push!(LOAD_PATH,"/PATH/to/source/code/RSA")
using RSA
```
To test whether the module was loaded you can use the following command in the Julia REPL:
```
?RSA
```