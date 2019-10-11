#!/usr/bin/env bash
# Setup a folder for a study with the given ID

# Setup dir with yaml meta file and raw subdir
mkdir "$1"
touch "$1/$1.yaml"
mkdir "$1/raw"

cp "../../bin/standardise_study_template.R" "$1/standardise_$1.R" || echo "Couldn't copy template R script"
