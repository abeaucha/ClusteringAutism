#!/bin/bash

#On MICe machines: Remove all modules
module purge
module use /hpf/largeprojects/MICe/tools/modulefiles/linux-centos7-sandybridge/
module load python-3.9.13-gcc-8.2.0-kagensa

#If venv does not exist, create it
if [ ! -d ".venv" ]; then
	echo "Initializing python virtual environment..."
	python3 -m venv .venv
	. .venv/bin/activate
	echo "Upgrading pip..."
	pip install pip --upgrade
	echo "Installing python packages..."
	pip3 install -r python_reqs_hpc.txt
	deactivate
fi

#Load necessary modules
#NOTE: minc-stuffs module break python virtual environment
#module load \
#	minc-toolkit \
#	r/3.6.3 \
#	r-packages \
#	ants
module load minc-toolkit-1.9.18.2-gcc-8.2.0-p4ac6lj
module load r-3.6.3-gcc-8.2.0-gypogri
module load r-packages/20220704
module load ants-20220609-gcc-8.2.0-lnn77l7 ants-2.4.0-gcc-8.2.0-z6w27gh

#Activate the python venv
source .venv/bin/activate
SRCPATH="$PWD/src"
PYTHONPATH="$SRCPATH:$PYTHONPATH"
PATH="$SRCPATH:$SRCPATH/pipelines:$PATH"
RMINC_BATCH_CONF=/hpf/largeprojects/MICe/tools/r-packages/20220704/packratting/packrat/lib/x86_64-pc-linux-gnu/3.6.3/RMINC/parallel/slurm_batchtools.R

export SRCPATH
export PYTHONPATH
export PATH
export RMINC_BATCH_CONF