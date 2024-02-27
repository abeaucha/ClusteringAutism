#!/bin/bash

# Use MICe module files
module purge
module use /hpf/largeprojects/MICe/tools/modulefiles

# Import the MICe environment
module load mice-env/1.1.0

# If venv does not exist, create it
if [ ! -d ".venv_hpc" ]; then
	echo "Initializing python virtual environment..."
	python3 -m venv .venv_hpc
	. .venv_hpc/bin/activate
	echo "Upgrading pip..."
	pip install pip --upgrade
	echo "Installing python packages..."
	pip3 install -r python_reqs_hpc.txt
	deactivate
fi

# Activate the venv
source .venv_hpc/bin/activate

# Set environment variables
SRCPATH="$PWD/src"
PYTHONPATH="$SRCPATH:$PYTHONPATH"
PATH="$SRCPATH:$SRCPATH/pipelines:$PATH"
RMINC_BATCH_CONF=/hpf/largeprojects/MICe/tools/RMINC/1.5.2.2_2/packrat/lib/x86_64-pc-linux-gnu/3.6.1/RMINC/parallel/slurm_batchtools.R

export SRCPATH
export PYTHONPATH
export PATH
export RMINC_BATCH_CONF
