# Installation instructions

## Overview
1. Clone this repo - `git clone https://github.com/allydunham/aa_subtypes.git`
2. Install neccessary R & python packages
3. Install other dependancies
4. Setup local environment settings
5. Download raw data files

The dependancies can be installed with your prefered package manager. I used a combination
of local python virtual environments (managed with pipenv) and a `conda` environment on
our computer cluster for running the heavier programs. Conda can be used for all but a few
packages, which were just installed manually.

## 2. Required Python Packages


## 3. Other Dependancies


## 4. Local Environment
1. Add this repos `src` directory to your PYTHONPATH
2. Update the Porter5 and uniref90 paths in snakemake.yaml
3. Update cluster.yaml if you have different preferences for running snakemake in 
   cluster mode. I used [snakemake-lsf](https://github.com/Snakemake-Profiles/snakemake-lsf)
   to manage cluster mode, but the required setup will depend on your computer cluster setup.

## 5. Raw Data Files


