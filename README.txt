### Review of PhyloFormer - a deep learning framework to build phylogenetic trees

As PhyloFormer does not come with an executable, I clone their package from: https://github.com/lucanest/Phyloformer 
and then insert a link to this folder (where the `predict.py` script is located) into my jupyter notebook.

In linux this can be done as follows:
```
ln -s path_to_Phyloformer/predict.py PhyloFormer_predict
```
I then initialize the conda environment (this is the environment for PhyloFormer with the additional packages Bio and augur):

```
conda install -n base -c conda-forge mamba
mamba env create -f environment.yml
conda activate phylo
```
I perform simulations in [julia] (https://julialang.org/downloads/), to set this up see the instructions in the TreeSimulations folder.

A detailed description of the simulation method, as well as the test framework is given in the jupyter notebook.