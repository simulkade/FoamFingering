# FoamFingering
## The models and results of our SPE-173193 paper
To run these files, you need to install [Julia](http://julialang.org/downloads/) and [Anaconda](https://www.continuum.io/downloads) and the following Julia packages:
  + IJulia
  + JFVM
  + PyPlot
  + Dierckx
  + DataFrames
  + ProgressMeter
  + Roots
  + JLD

Just open `Julia` and run `Pkg.add("package_name_from_above")`.  
You can open the notebooks in `Jupyter` and run the code. It solves water-alternating-gas (WAG) and surfactant-alternating-gas (SAG) processes in 1D analytically and numerically, and in 2D numerically.
