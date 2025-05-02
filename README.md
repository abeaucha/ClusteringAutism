

# Installation

Currently, this project can only be installed on ARM64 Apple Silicon computers. 

## Pre-requisite: conda

This project requires a conda distribution for environment and package management. A minimal installation of conda can be obtained via [Miniforge](https://github.com/conda-forge/miniforge). 

Download the installation script from Github:
```{zsh}
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSx-arm64.sh"
```

Install Miniforge into your preferred destination (e.g. `~/opt/miniforge3`):
```{zsh}
bash Miniforge3-MacOSx-arm64.sh -b -p ~/opt/miniforge3
```

Initialize conda for `zsh`:
```{zsh}
~/opt/miniforge3/bin/conda init zsh
```

Disable auto-activation of base environment on startup:
```{zsh}
~/opt/miniforge3/bin/conda config --set auto_activate_base false
```

Restart Terminal or re-source the startup files:
```{zsh}
source ~/.zprofile
source ~/.zshrc
```


## Pre-requisite: GFortran

Installing `RMINC` requires a fortran compiler. Since Apple uses clang by default, most users will not have a fortran
compiler installed. GFortran can be installed as part of the Gnu Compiler Collection via Homebrew. 

```{zsh}
brew install gcc
```

Regardless of how you install GFortran you will need to tell
R where to find it. If you installed via Homebrew you can run the
following in Terminal:

```{zsh}
mkdir ~/.R
touch ~/.R/Makevars
echo "FLIBS=-L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14/" >> ~/.R/Makevars
```

This assumes that Homebrew is installed in the standard directory `/opt/homebrew/`. Update the path as needed if you installed Homebrew elsewhere.

## Build the project environment

Install the project (conda) environment by sourcing the setup script:

```{zsh}
./setup.sh
```

This creates two conda environments:
1. The main project environment named `clustering-autism-env` containing R and Python.  
2. A support environment named `clustering-autism-minc-env` containing the minc-toolkit-v2 installation

## Activate the project environment

To run project files (e.g. pipelines), first activate the project environment:

```{zsh}
conda activate clustering-autism-env
```

In addition to specifying the versions of R and Python, this defines a number of environment variables
- `PROJECTPATH`: The path to the local copy of the project repository.
- `SRCPATH`: The path to the project source files (i.e. subdirectory `src`)
- `MINC_TOOLKIT`: The path to the minc-toolkit-v2 installation in the `clustering-autism-minc-env` environment.


# Demonstration

