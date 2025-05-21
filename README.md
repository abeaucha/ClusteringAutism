
# Overview

This repository contains all the source code pertaining to the manuscript:

Ellegood J, Beauchamp A, Yee Y et al. 2025. "Assigning Targetable Molecular Pathways to Transdiagnostic Subgroups Across Autism and
Related Neurodevelopmental Disorders".

The manuscript is currently openly accessible on [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.03.04.641443v1).

# System requirements

## Hardware requirements

This project is supported for Linux x86_64 computers and Apple ARM64 computers.

## Software requirements

### Operating systems

The following operating systems are currently supported:
- macOS Sequoia 15.4.1.
- Ubuntu 22.04.5 LTS

### Software

The following software are required to build the project environment:
- conda
- GFortran

Suggested installation instructions can be found below. 

The project environment is constructed using conda and installs the following software:
- R v4.4.2
- Python v3.12.10
- minc-toolkit-v2 v1.9.19

### R dependencies

```text
devtools v2.4.5
tidyverse v2.0.0
lme4 v1.1-37
visNetwork v2.1.2
rjson v0.2.23
dt v0.33
doParallel v1.0.17
parallelly v1.43.0
cowplot v1.1.3
RSpectra v0.16.2
Rcpp v1.0.14
RcppEigen v0.3.4.0.2
RcppArmadillo v14.4.1-1
optparse v1.7.5
SNFtool v2.3.1
doSNOW v1.0.20
rcartocolor v2.1.1
gridExtra v2.3
patchwork v1.3.0
viridis v0.6.5
pROC v1.18.5
RMINC v1.5.3.0
MRIcrotome v0.2.0
```

R dependencies are stored in the file `R_packages_conda.txt` and are installed automatically as part of the environment setup process.

### Python dependencies

```text
NumPy v2.2.5
SciPy v1.15.2
pandas v2.2.3
tqdm v4.67.1
scikit-learn v1.6.1
matplotlib v3.10.1
seaborn v0.13.2
utils v1.0.2
pyminc v0.57
```

Python dependencies are stored in the files `python_packages_conda.txt` and `python_packages_pip.txt` and are installed automatically as part of the environment setup process.

# Installation guide

## Requirement: conda

This project requires a conda distribution for environment and package management. A minimal installation of conda can be obtained via [Miniforge](https://github.com/conda-forge/miniforge). 

Download the installation script from Github:
- Apple ARM64: Miniforge3-MacOSx-arm64.sh
- Linux x86_64: Miniforge3-Linux-x86_64.sh

```shell
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSx-arm64.sh"
```

Install Miniforge into your preferred destination (e.g. `~/opt/miniforge3`):
```shell
bash Miniforge3-MacOSx-arm64.sh -b -p ~/opt/miniforge3
```

Initialize conda for `zsh` or `bash`:
```shell
~/opt/miniforge3/bin/conda init zsh
```

Disable auto-activation of base environment on startup:
```shell
~/opt/miniforge3/bin/conda config --set auto_activate_base false
```

Restart Terminal or re-source the startup files. 
```shell
# zsh startup files
source ~/.zprofile
source ~/.zshrc

# bash startup files
source ~/.bashrc
```


## Requirement: GFortran

Installing `RMINC` requires a Fortran compiler. Since Apple uses clang by default, most users will not have a Fortran
compiler installed. 

On macOS, GFortran can be installed as part of the GNU Compiler Collection via Homebrew: 

```shell
brew install gcc
```

On Linux systems, GFortran can be installed using the default package manager:

```shell
apt install gfortran
```

Regardless of how you install GFortran you will need to tell
R where to find it. If you installed via Homebrew you can run the
following in Terminal:

```shell
mkdir ~/.R
touch ~/.R/Makevars
echo "FLIBS=-L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14/" >> ~/.R/Makevars
```

This assumes that Homebrew is installed in the standard directory `/opt/homebrew/`. Update the path as needed if you installed Homebrew elsewhere.

## Build the project environment

The typical installation time for the project environment is approximatley 10 minutes.

The project (conda) environment is constructed by sourcing the setup script:

```shell
./setup.sh
```

This creates two conda environments:
1. The main project environment named `clustering-autism-env` containing R and Python.  
2. A support environment named `clustering-autism-minc-env` containing the minc-toolkit-v2 installation

## Activate the project environment

To run project files (e.g. pipelines), first activate the project environment:

```shell
conda activate clustering-autism-env
```

In addition to specifying the versions of R and Python, this defines a number of environment variables
- `PROJECTPATH`: The path to the local copy of the project repository.
- `SRCPATH`: The path to the project source files (i.e. subdirectory `src`)
- `MINC_TOOLKIT`: The path to the minc-toolkit-v2 installation in the `clustering-autism-minc-env` environment.


# Demo

A demo data set (~ 5GB) is available upon request. The pipeline demo can be found in the shell script `demo.sh`. Execute the demo by sourcing the script from the project root:

```shell
./demo.sh
```

The expected output are the following pipeline output sub-directories and files:

```text
data
|-- human
|   |-- derivatives
|       |-- metadata.csv
|       |-- 001
|           |-- demographics.csv
|           |-- centroids
|           |   |-- resolution_3.0
|           |       |-- absolute
|           |       |   |-- centroid_nk_2_k_*_.mnc
|           |       |-- relative
|           |           |-- centroid_nk_2_k_*_.mnc
|           |-- clusters
|           |   |-- resolution_3.0
|           |       |-- affinity.csv
|           |       |-- clusters.csv
|           |-- effect_sizes
|           |   |-- resolution_3.0
|           |       |-- absolute
|           |       |   |-- effect_sizes.csv
|           |       |   |-- participant_*.mnc
|           |       |-- relative
|           |           |-- effect_sizes.csv
|           |           |-- participant_*.mnc
|           |-- jacobians
|               |-- absolute
|               |   |-- participant_*.mnc
|               |-- relative
|                   |-- participant_*.mnc
|-- mouse
|   |-- derivatives
|       |-- metadata.csv
|       |-- 001
|           |-- model_list_Feb2023.RData
|           |-- scanbase_scans_autism_filtered_Feb2023.csv
|           |-- centroids
|           |   |-- resolution_0.2
|           |       |-- absolute
|           |       |   |-- centroid_nk_2_k_*_.mnc
|           |       |-- relative
|           |           |-- centroid_nk_2_k_*_.mnc
|           |-- clusters
|           |   |-- resolution_0.2
|           |       |-- affinity.RData
|           |       |-- clusters.csv
|           |-- effect_sizes
|               |-- 200
|                   |-- *_ES_Absolute_200.mnc
|                   |-- *_ES_Relative_200.mnc
|-- cross_species
    |-- metadata.csv
    |-- 001
        |-- similarity
            |-- centroid_pairs.csv
            |-- similarity.csv
```

Notes:
- The human Jacobian and effect size MINC files run from "participant_1.mnc" to "participant_120.mnc". 
- The mouse effect size MINC image files are prefixed by the study name and genotype code, and the directory should contain 50 images in total.
- The human and mouse centroid sub-directories should contain two MINC files: "centroid_nk_2_k_1.mnc" and "centroid_nk_2_k_2.mnc".

The demo pipeline should execute in approximately 150 seconds.


# Instructions for use

This repository provides the pipelines and source code used to generate the results in our manuscript. There are four main pipelines in the processing workflow:
1. `process_human_images.py`
2. `process_mouse_images.py`
3. `compute_cluster_similarity.py`
4. `permute_cluster_similarity.py`

All pipelines can be found in the sub-directory `src/pipelines`.

### 1. Processing human images

Human images are processed using the `process_human_images.py` pipeline, which implements the following:
1. Generates absolute and relative effect size images for the participants
2. Generates clusters of participants based on the effect size images
3. Evaluates absolute and relative centroid images for each of the clusters

```shell
process_human_images.py \
  --params-id 001 \
  --pipeline-dir data/human/derivatives/ \
  --input-dir data/human/registration/jacobians/ \
  --demographics data/human/registration/subject_info/demographics.csv \
  --mask data/human/registration/reference_files/mask_0.8mm.mnc \
  --datasets POND SickKids \
  --es-group patients \
  --es-method normative-growth \
  --es-batch Site Scanner \
  --es-df 3 \
  --cluster-resolution 3.0 \
  --cluster-nk-max 10 \
  --cluster-metric correlation \
  --cluster-K 10 \
  --cluster-sigma 0.5 \
  --cluster-t 20 \
  --centroid-method mean \
  --execution local \
  --nproc 8
```

### 2. Processing mouse images

Mouse images are processed using the `process_mouse_images.R` pipeline, which implements the same steps as the human processing pipeline on the mouse data:
1. Generates absolute and relative effect size images for the mouse models
2. Generates clusters of mouse models based on the effect size images
3. Evaluates absolute and relative centroid images for each of the clusters

```shell
process_mouse_images.R
```

### 3. Evaluating the similarity of mouse and human clusters

Once human and mouse clusters have been identified, the cross-species transcriptomic similarity between the clusters is evaluated using the `computer_cluster_similarity.py` pipeline.

```shell
compute_cluster_similarity.py \
  --pipeline-dir data/cross_species/ \
  --params-id 001 \
  --species human mouse \
  --input-dirs data/human/derivatives/ data/mouse/derivatives/ \
  --input-params-ids 001 001 \
  --expr-dirs data/human/expression data/mouse/expression \
  --masks data/human/registration/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
  --microarray-coords data/human/expression/AHBA_microarray_coordinates_study.csv \
  --gene-space average-latent-space \
  --n-latent-spaces 100 \
  --metric correlation \
  --signed true \
  --threshold top_n \
  --threshold-value 0.2 \
  --threshold-symmetric true \
  --jacobians absolute relative \
  --execution local \
  --nproc 8
```

This pipeline can also be used to evaluate the transcriptomic similarity between clusters within the same species (e.g. human-human).

### 4. Evaluating permutation resampling of the cluster similarity for empiricial null distributions

Once the true cluster similarities have been evaluated, permutation resampling can be applied to build up empirical null distributions for the similarity of cluster pairs. This is done using the `permute_cluster_similarity.py` pipeline.

```shell
permute_cluster_similarity.py \
--pipeline-dir data/cross_species/ \
--param-id 001 \
--input-dirs data/human/derivatives/ data/mouse/derivatives/ \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/AHBA_microarray_coordinates_study.csv \
--permutations-start 1 \
--permutations-n 500 \
--off-diagonal 1 \
--execution local \
--nproc 8
```



## Reproduction instructions

The quantitative results in our manuscript were generated using the following pipeline configurations.

### 1. Processing human and mouse images

Processing the human images in the POND+ dataset:
```shell
process_human_images.py \
  --params-id 700
  --pipeline-dir data/human/derivatives/ \
  --input-dir data/human/registration/jacobians_resampled/resolution_0.8/ \
  --demographics data/human/registration/subject_info/demographics.csv \
  --mask data/human/registration/reference_files/mask_0.8mm.mnc \
  --datasets POND SickKids \
  --es-group patients \
  --es-method normative-growth \
  --es-batch Site Scanner \
  --es-df 3 \
  --cluster-resolution 3.0 \
  --cluster-nk-max 10 \
  --cluster-metric correlation \
  --cluster-K 10 \
  --cluster-sigma 0.5 \
  --cluster-t 20 \
  --centroid-method mean \
  --execution slurm \
  --nproc 8 \
  --slurm-njobs 300 \
  --slurm-time 120 \
  --slurm-mem 16G
```

Processing the human images in the HBN dataset:
```shell
process_human_images.py \
  --params-id 013 \
  --pipeline-dir data/human/derivatives/ \
  --input-dir data/human/registration/jacobians_resampled/resolution_0.8/ \
  --demographics data/human/registration/subject_info/demographics.csv \
  --mask data/human/registration/reference_files/mask_0.8mm.mnc \
  --datasets HBN \
  --es-group patients \
  --es-method normative-growth \
  --es-batch Site Scanner \
  --es-df 3 \
  --cluster-resolution 3.0 \
  --cluster-nk-max 10 \
  --cluster-metric correlation \
  --cluster-K 10 \
  --cluster-sigma 0.5 \
  --cluster-t 20 \
  --centroid-method mean \
  --execution slurm \
  --nproc 8 \
  --slurm-njobs 300 \
  --slurm-time 120 \
  --slurm-mem 16G
```

Processing the mouse images:
```shell
process_mouse_images.R
```

### 2. Evaluating cluster similarity

Evaluating the similarity between POND+ human clusters and mouse clusters:

```shell
compute_cluster_similarity.py \
  --pipeline-dir data/cross_species/ \
  --params-id 375 \
  --species human mouse \
  --input-dirs data/human/derivatives/ data/mouse/derivatives/ \
  --input-params-ids 700 001 \
  --expr-dirs data/human/expression data/mouse/expression \
  --masks data/human/registration/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
  --microarray-coords data/human/expression/AHBA_microarray_coordinates_study.csv \
  --gene-space average-latent-space \
  --n-latent-spaces 50 \
  --metric correlation \
  --signed true \
  --threshold top_n \
  --threshold-value 0.2 \
  --threshold-symmetric true \
  --jacobians absolute relative \
  --execution slurm \
  --slurm-njobs 300 \
  --slurm-mem 16G \
  --slurm-time 8:00:00
```

Evaluating the similarity between HBN human clusters and mouse clusters:

```shell
compute_cluster_similarity.py \
  --pipeline-dir data/cross_species/ \
  --params-id 861 \
  --species human mouse \
  --input-dirs data/human/derivatives/ data/mouse/derivatives/ \
  --input-params-ids 013 001 \
  --expr-dirs data/human/expression data/mouse/expression \
  --masks data/human/registration/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
  --microarray-coords data/human/expression/AHBA_microarray_coordinates_study.csv \
  --gene-space average-latent-space \
  --n-latent-spaces 50 \
  --metric correlation \
  --signed true \
  --threshold top_n \
  --threshold-value 0.2 \
  --threshold-symmetric true \
  --jacobians absolute relative \
  --execution slurm \
  --slurm-njobs 300 \
  --slurm-mem 16G \
  --slurm-time 8:00:00
```

Evaluating the similarity between POND+ clusters and HBN clusters:

```shell
compute_cluster_similarity.py \
  --pipeline-dir data/cross_species/ \
  --params-id 779 \
  --species human human \
  --input-dirs data/human/derivatives/ data/human/derivatives/ \
  --input-params-ids 700 013 \
  --expr-dirs data/human/expression data/human/expression \
  --masks data/human/registration/reference_files/mask_0.8mm.mnc data/human/registration/reference_files/mask_0.8mm.mnc \
  --microarray-coords data/human/expression/AHBA_microarray_coordinates_study.csv \
  --gene-space average-latent-space \
  --n-latent-spaces 50 \
  --metric correlation \
  --signed true \
  --threshold top_n \
  --threshold-value 0.2 \
  --threshold-symmetric true \
  --jacobians absolute relative \
  --execution slurm \
  --slurm-njobs 300 \
  --slurm-mem 16G \
  --slurm-time 8:00:00
```

### 3. Permuting the cluster similarity 

Permuting the similarity between POND+ human clusters and mouse clusters:

```shell
permute_cluster_similarity.py \
--pipeline-dir data/cross_species/ \
--param-id 375 \
--input-dirs data/human/derivatives/ data/mouse/derivatives/ \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/AHBA_microarray_coordinates_study.csv \
--permutations-start 1 \
--permutations-n 500 \
--off-diagonal 1 \
--execution slurm \
--slurm-njobs 300 \
--slurm-mem 16G \
--slurm-time 8:00:00
```

Permuting the similarity between HBN human clusters and mouse clusters:

```shell
permute_cluster_similarity.py \
--pipeline-dir data/cross_species/ \
--param-id 861 \
--input-dirs data/human/derivatives/ data/mouse/derivatives/ \
--expr-dirs data/human/expression data/mouse/expression \
--masks data/human/registration/reference_files/mask_0.8mm.mnc data/mouse/atlas/coronal_200um_coverage_bin0.8.mnc \
--microarray-coords data/human/expression/AHBA_microarray_coordinates_study.csv \
--permutations-start 1 \
--permutations-n 500 \
--off-diagonal 1 \
--execution slurm \
--slurm-njobs 300 \
--slurm-mem 16G \
--slurm-time 8:00:00
```

Permuting the similarity between POND+ clusters and HBN clusters:

```shell
permute_cluster_similarity.py \
--pipeline-dir data/cross_species/ \
--param-id 779 \
--input-dirs data/human/derivatives/ data/human/derivatives/ \
--expr-dirs data/human/expression data/human/expression \
--masks data/human/registration/reference_files/mask_0.8mm.mnc data/human/registration/reference_files/mask_0.8mm.mnc \ 
--microarray-coords data/human/expression/AHBA_microarray_coordinates_study.csv \
--permutations-start 1 \
--permutations-n 500 \
--off-diagonal 1 \
--execution slurm \
--slurm-njobs 300 \
--slurm-mem 16G \
--slurm-time 8:00:00
```

The source code for all downstream analyses and visualizations related to figures in the manuscript can be found in the subdirectory `figures/v3`.