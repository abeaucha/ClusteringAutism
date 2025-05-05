
# System requirements

## Hardware requirements

This project is supported for ARM64 Apple Silicon computers.

## Software requirements

### Operating systems

The following operating systems are supported:
- macOS Sequoia 15.4.1.

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

R dependencies are installed automatically as part of the environment setup process.

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

Python dependencies are installed automatically as part of the environment setup process.

---

# Installation guide

## Requirement: conda

This project requires a conda distribution for environment and package management. A minimal installation of conda can be obtained via [Miniforge](https://github.com/conda-forge/miniforge). 

Download the installation script from Github:
```shell
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSx-arm64.sh"
```

Install Miniforge into your preferred destination (e.g. `~/opt/miniforge3`):
```shell
bash Miniforge3-MacOSx-arm64.sh -b -p ~/opt/miniforge3
```

Initialize conda for `zsh`:
```shell
~/opt/miniforge3/bin/conda init zsh
```

Disable auto-activation of base environment on startup:
```shell
~/opt/miniforge3/bin/conda config --set auto_activate_base false
```

Restart Terminal or re-source the startup files:
```shell
source ~/.zprofile
source ~/.zshrc
```


## Requirement: GFortran

Installing `RMINC` requires a Fortran compiler. Since Apple uses clang by default, most users will not have a Fortran
compiler installed. GFortran can be installed as part of the GNU Compiler Collection via Homebrew. 

```shell
brew install gcc
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

Construct the project (conda) environment by sourcing the setup script:

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

---

# Demo

A demo pipeline is available in the shell script `pipeline_demo.sh`. Execute the demo by sourcing the script.

```shell
./pipeline_demo.sh
```

The expected output of the demo are the following pipeline output sub-directories and files:

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

