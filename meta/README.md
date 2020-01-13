# Meta Data folder

This folder contains various bits of metadata required to execute the project
and a few summary tables output by the pipeline.

Metadata:
* **structures.yaml** - details the SwissMODEL structures chosen for the studied proteins. 
  Each gene is linked to a SwissMODEL template ID, the type of structure it is (e.g. 
  x-ray, homology model etc.) and a list of sections to use. Each section is defined by
  a chain and region, with an optional offset for sequence numbering compared to Uniprot.
  This manual system could theoretically be replaced by automated model selection via the
  SwissMODEL API, but the complexity of the choice and the small number of proteins made 
  manual selection easier. This file is used downstream to select regions of the PDB to
  process (e.g. with FoldX and Naccess) and to convert results to Uniprot sequence numbering.
  If a new study is added for a new protein a model must be chosen and added here.
* **study\_template.yaml** - Template YAML file for adding new studies to the project
* **subtypes/** - folder containing YAML config files for the various clustering approaches 
  to extracting amino acid subtypes.
* **residue\_hydrophobicity.tsv** - Average amino acid hydrophobicity table, sourced from 
  [Bandyopadhyay & Mehler (2008)](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.21958)
* **fasta/** - can contain .fa files of the form {species/strain}\_{gene}.fa that act as 
  master copies when validating study config files. Not required for normal execution.
  In general the fasta files are sourced from Uniprot, although some come from the 
  studies directly or have manual edits based on the mutations used in a study.

Generated Files:
* **overall\_summary** - Summary of the project overall, giving the total number of 
  studies, genes, etc. processed.
* **gene\_summary.tsv** - Table summarising the properties of the genes included in the
  project, including stats such as the number of mutants and mutant coverage.
* **study\_summary.tsv** - Table summarising the properties of the studies included in 
  the project.
