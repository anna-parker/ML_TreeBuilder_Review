Simulations are performed using julia. As I use non registered packages I would suggest installing packages from scratch. 
To do this run `julia --project=.` in this directory and executing the following commands
precompile the simulation script for faster simulations:
```
using Pkg
Pkg.add("TreeTools")
Pkg.add("ArgParse")
Pkg.add("BioSequences")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Combinatorics")
Pkg.add("Distributions")
Pkg.add("Random")
Pkg.add(url="https://github.com/diegozea/ROC.jl", rev="master")
Pkg.add(url="https://github.com/PierreBarrat/ARGTools")
Pkg.add(url="https://github.com/anna-parker/TreeAlgs-Fork")
```
For faster simulations I additionally create a system image (note this has to be rerun if any modifications are made to the script or packages).
```
using PackageCompiler

create_sysimage(
    [:ARGTools, :TreeAlgs,
    :TreeTools, :CSV, :ArgParse,
    :Random, :Combinatorics, :Distributions],
    sysimage_path="tree_simulations.so",
    precompile_execution_file="Performance_Influenza/PhyloFormerSimulations.jl" # the new line
)
```
Note I also include a text file with c values for the influenza model of coalescence - these have been experimentally determined 
to result in trees of the the desired resolution, the file is used in the simulation script. 
