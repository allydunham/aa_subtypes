#!/usr/bin/env bash
# Setup a folder for a study with the given ID

# Config
proj_root="$HOME/phd/subtypes"
template_yaml="$proj_root/meta/study_template.yaml"
template_script="$proj_root/bin/standardise_study_template.R"

# Setup dir with yaml meta file and raw subdir
mkdir "$1" || exit
mkdir "$1/raw"
cp "$template_yaml" "$1/$1.yaml" || echo "Couldn't copy template YAML"
cp "$template_script" "$1/standardise_$1.R" || echo "Couldn't copy template R script"
