# Exploring amino acid functions in a deep mutational landscape 
Alistair Dunham & Pedro Beltrao, Molecular Systems Biology, 2021

This repo contains the code used to generate the results of [Dunham & Beltrao (2021)](https://doi.org/10.15252/msb.202110305). The accompanying [DeepScanScape](https://github.com/allydunham/deepscanscape) package allows users to apply a similar analysis to their own data.

## Abstract

Amino acids fulfil a diverse range of roles in proteins, each utilising its chemical properties in different ways in different contexts to create required functions.
For example, cysteines form disulphide or hydrogen bonds in different circumstances and charged amino acids do not always make use of their charge. The repertoire of amino acid functions and the frequency at which they occur in proteins remains understudied.
Measuring large numbers of mutational consequences, which can elucidate the role an amino acid plays, was prohibitively time-consuming until recent developments in deep mutational scanning.
In this study, we gathered data from 28 deep mutational scanning studies, covering 6,291 positions in 30 proteins, and used the consequences of mutation at each position to define a mutational landscape.
We demonstrated rich relationships between this landscape and biophysical or evolutionary properties.
Finally, we identified 100 functional amino acid subtypes with a data-driven clustering analysis and studied their features, including their frequencies and chemical properties such as tolerating polarity, hydrophobicity or being intolerant of charge or specific amino acids.
The mutational landscape and amino acid subtypes provide a foundational catalogue of amino acid functional diversity, which will be refined as the number of studied protein positions increases.

## Project Overview

```plaintext
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
|   |-- 0_data
|   |-- 1_landscape
|   |-- 2_subtypes
|   |-- 3_continuous
|   `-- 4_figures                                - Final paper figures (including subparts)
|-- logs                                         - Logs from snakemake
|-- meta                                         - Various metadata files (uniprot domains, residue hydrophobicity etc.)
|   `-- subtypes                                 - Configuration files for subtyping algorithms
`-- src                                          - Python and R modules with shared functions
```

## Running the pipeline

The pipeline is managed by Snakemake, with various convenience rules supplied:

* all - Generate the subtypes and figures for the paper along with a few additional
analyses.
* full\_analysis - Perform all these and various additional analyses, primarily generating
subtypes using other tested methods and generating many more summary plots for the various
stages of the pipeline. This includes mapping various factors onto the protein structures
and longer clustering algorithms, so takes significantly more time.
* full\_clean - Remove all generated files, including the main data files (combined\_mutational\_scans.tsv).
This is mainly useful when running the pipeline from raw data downloaded from individual studies.
* quick\_clean - Remove all generated files apart from those that take a long time to
generate (FoldX results, SIFT results etc.) and the main data files.

### Setup Overview

1. Clone this repo - `git clone https://github.com/allydunham/aa_subtypes.git`
2. Download data files (optional)
3. Install dependancies
4. Setup local environment settings

I used a combination of local python virtual environments and a
`conda` environment on our computer cluster for running the heavier programs.
Conda can be used for all but a few packages, which were just installed manually.

### 2. Data Files

To run the full pipeline the required files need to be downloaded from each study website,
as indicated in the study YAML files, and placed in the appropriate `raw/` folder.
However, the data analysis pipeline can be run using the supplied aggregated data files,
containing processed data from each study and associated results from various tools (SIFT,
FoldX etc.).

For analyses using the clinvar dataset clinvar\_20200106.vcf.gz and
clinvar\_20200106.vcf.gz.tbi must be downloaded from the Clinvar website and placed in
data/clinvar.
Newer releases can also be used if the appropriate script is updated
(bin/analysis/0\_data/clinvar\_classification.R).
If these files are not available Snakemake will proceed without this analysis, which
is supplementary to the main focus.

### 3. Dependancies

Python 3.7 and R 3.6.3 were used

#### Python packages

* biopython
* colorama
* numpy
* matplotlib
* pandas
* ruamel.yaml

In addition the python modules in src/ must be able to be imported.

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

Clinvar analyses also requires:

* BSgenome.Hsapiens.UCSC.hg38
* Homo.sapiens
* TxDb.Hsapiens.UCSC.hg38.knownGene
* VariantAnnotation

#### Other Tools

* SIFT4G (patched to output 4 decimal places)
* FoldX 5
* Naccess 2.1.1
* Porter5
* Pymol 2.3.0 open source (and requires pymol2 to be importable)

Of these only Pymol is required for the analysis part of the pipeline, and then
only when running the complete pipeline which generates various plots of features mapped
to protein structures.
The standard `Snakemake` command doesn't produce these extra plots by default so doesn't require
pymol.
The results of the other tool are all pregenerated and included in the supplied data files,
so are only required to run the pipeline from raw data.

### 4. Local Environment

1. Add this repos `src` directory to your $PYTHONPATH
2. Update the Porter5 and uniref90 paths in snakemake.yaml
3. Update cluster.yaml if you have different preferences for running snakemake in
   cluster mode. I used [snakemake-lsf](https://github.com/Snakemake-Profiles/snakemake-lsf)
   to manage cluster mode, but the required setup will depend on your setup.
