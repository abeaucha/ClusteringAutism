#!/bin/zsh

# Make conda tools visible
source $(conda info --base)/etc/profile.d/conda.sh

conda remove --name clustering-autism-minc-env --all -y
conda remove --name clustering-autism-env --all -y
conda clean --all -y
