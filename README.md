# Characterising amino acid positions and their subtypes using deep mutational landscapes
This repo contains the code for the publication XXX

## Abstract

## Project Overview

```
subtypes
|-- Snakefile                                    - Main pipeline
|-- bin                                          - Scripts (mainly R & Python), split into various sections
|   `-- pipeline                                 - Snakemake rules for subsections of the pipeline
|-- data                                         - Data files
|   |-- studies                                  - Directory for each input study with raw data, processing
|   |   |-- study                                  script and processed data tsv
|   |   |   `-- raw
|   |-- subtypes                                 - Output of position subtype clustering procedures
|   |-- tools                                    - Various folders with the output of each tool on each gene
|   `-- (long_)combined_mutational_scans.tsv     - Processed, combined data tables, with full details for each mutation or averaged over the position
|-- docs                                         - Document describing final subtypes
|-- figures                                      - Directories of figures for each stage of the project
|   |-- 0\_data
|   |-- 1\_landscape
|   |-- 2\_subtypes
|   |-- 3\_continuous
|   `-- 4\_figures                               - Final paper figures (including subparts)
|-- logs                                         - Logs from snakemake
|-- meta                                         - Various metadata files (uniprot domains, residue hydrophobicity etc.)
|   `-- subtypes                                 - Configuration files for subtyping algorithms
`-- src                                          - Python and R modules with shared functions
```

## Running the pipeline
### Overview
1. Clone this repo - `git clone https://github.com/allydunham/aa_subtypes.git`
2. Download raw data files
3. Install dependancies
4. Setup local environment settings

I used a combination of local python virtual environments (managed with pipenv) and a 
`conda` environment on our computer cluster for running the heavier programs.
Conda can be used for all but a few packages, which were just installed manually.

### 2. Raw Data Files
To run the full pipeline the required files need to be downloaded from each study website, 
as indicated in the study YAML files, and placed in the appropriate `raw/` folder.
However, the data analysis pipeline can be run using the supplied aggrgated data files, 
continaing processed data from each study and associated results from various tools (SIFT,
FoldX etc.).

### 3. Dependancies

Python 3.7 and R 3.6.3 were used

#### Python packages
* Biopython
* colorama
* numpy
* matplotlib
* pandas
* ruamel.yaml

In addition the python modules in src/ must be able to be imported.
These are also listed in the Pipfile.

#### R packages
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

#### Other Tools
* SIFT4G (patched to output 4 decimal places)
* FoldX 5 
* Naccess 2.1.1
* Porter5 
* Pymol 2.3.0 open source (and requires pymol2 to be importable)

Of these only Pymol is required for the analysis part of the pipeline, with the results 
of the other tool being included in the supplied data files.

### 4. Local Environment
1. Add this repos `src` directory to your PYTHONPATH
2. Update the Porter5 and uniref90 paths in snakemake.yaml
3. Update cluster.yaml if you have different preferences for running snakemake in
   cluster mode. I used [snakemake-lsf](https://github.com/Snakemake-Profiles/snakemake-lsf)
   to manage cluster mode, but the required setup will depend on your setup.


