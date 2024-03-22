#!/bin/bash

# Use MICe module files
module purge
module use /hpf/largeprojects/MICe/tools/modulefiles

# Import the MICe environment
module load mice-env/1.1.0

# If venv does not exist, create it
if [ ! -d ".venv_hpc" ]; then
	echo "Initializing python virtual environment..."

	# Create the venv
	python3 -m venv .venv_hpc

	# Activate the venv
	. .venv_hpc/bin/activate

	# Upgrade pip
	echo "Upgrading pip..."
	pip install pip --upgrade

	# Install required python packages
	echo "Installing python packages..."
	pip3 install -r python_reqs_hpc.txt

	# Deactivate the venv
	deactivate
fi

# Activate the venv
source .venv_hpc/bin/activate

# Set environment variables
SRCPATH="$PWD/src"
PYTHONPATH="$SRCPATH:$PYTHONPATH"
PATH="$SRCPATH:$SRCPATH/drivers:$SRCPATH/pipelines:$PATH"
RMINC_BATCH_CONF=/hpf/largeprojects/MICe/tools/RMINC/1.5.2.2_2/packrat/lib/x86_64-pc-linux-gnu/3.6.1/RMINC/parallel/slurm_batchtools.R

export SRCPATH
export PYTHONPATH
export PATH
export RMINC_BATCH_CONF
