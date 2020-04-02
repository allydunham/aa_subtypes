# Characterising amino acid positions and their subtypes using deep mutational landscapes
This repo contains the code for the publication XXX

## Abstract

## Running the pipeline
## Overview
1. Clone this repo - `git clone https://github.com/allydunham/aa_subtypes.git`
2. Install dependancies
4. Setup local environment settings
5. Download raw data files

I used a combination of local python virtual environments (managed with pipenv) and a 
`conda` environment on our computer cluster for running the heavier programs.
Conda can be used for all but a few packages, which were just installed manually.

## Dependancies

Python 3.7 and R 3.6.3 were used

### Python packages
* Biopython
* colorama
* numpy
* matplotlib
* pandas
* ruamel.yaml

In addition the python modules in src/ must be able to be imported.
These are also listed in the Pipfile.

### R packages
* Biostrings
* broom
* cluster
* dbscan
* dynamicTreeCut
* GGally
* ggdendro
* ggpubr
* grid
* ggtext
* magrittr
* mclust
* multipanelfigure
* plotlistr (available at github.com/allydunham/plotlistr)
* png
* readxl
* rlang
* tblhelpr (available at github.com/allydunham/tblhelpr)
* tools
* tidyverse
* uwot
* yaml

## Other Tools
* SIFT4G
* FoldX 
* Naccess
* Porter5
* Pymol (requires pymol2 to be importable)

## 4. Local Environment
1. Add this repos `src` directory to your PYTHONPATH
2. Update the Porter5 and uniref90 paths in snakemake.yaml
3. Update cluster.yaml if you have different preferences for running snakemake in
   cluster mode. I used [snakemake-lsf](https://github.com/Snakemake-Profiles/snakemake-lsf)
   to manage cluster mode, but the required setup will depend on your computer cluster setup.

## 5. Raw Data Files
