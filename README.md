# Development | Team notes


💡 Feel free to create a new folder under `TeamNotes` with your github handle to store scripts, notes, etc.

## Contents

- [Tools](https://github.com/JuliaReach/dev/blob/master/Tools.md) -- Description of the tools that we use and collection of best practices.

## Recommended setup

First create a new folder `JuliaReach` and clone the repositories of interest, here we assume `LazySets` and `ReachabilityAnalysis`.

```bash
$ mkdir JuliaReach
$ cd JuliaReach
$ git clone git@github.com:JuliaReach/LazySets.jl.git
$ git clone git@github.com:JuliaReach/ReachabilityAnalysis.jl.git
```

Then open a new environment and add the JuliaReach packages for development, as well as several dependencies that are needed to build tests and documentation.

```julia
julia> using Pkg

julia> Pkg.develop(path="/home/mforets/JuliaReach/LazySets.jl")

julia> Pkg.develop(path="/home/mforets/JuliaReach/ReachabilityAnalysis.jl")

julia> Pkg.add(["AbstractTrees",
                "BenchmarkTools",
                "CDDLib",
                "DifferentialEquations",
                "DisplayAs",
                "Distributions",
                "Documenter",
                "Expokit",
                "ExponentialUtilities",
                "IntervalMatrices",
                "JLD2",
                "LaTeXStrings",
                "Literate",
                "Makie",
                "Optim",
                "Plots",
                "Polyhedra",
                "RecipesBase",
                "StaticArrays",
                "Symbolics",
                "TaylorModels"])
```
