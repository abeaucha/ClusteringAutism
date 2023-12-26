#!/bin/bash

#On MICe machines: Remove all modules
module purge

#If venv does not exist, create it
if [ ! -d ".venv" ]; then
	echo "Initializing python virtual environment..."

	#Create the venv
	python3 -m venv .venv

	#Activate the venv
	source .venv/bin/activate

	#Upgrade pip
	echo "Upgrading pip..."
	pip install pip --upgrade

	#Install necessary python packages
	echo "Installing python packages..."
	pip3 install -r python_reqs.txt

	#Deactivate the venv
	deactivate
fi

#Load necessary modules 
#NOTE: minc-stuffs module break python virtual environment
module load \
	minc-toolkit \
	r/3.6.3 \
	r-packages \
	ants

#Activate the python venv
source .venv/bin/activate
export SRCPATH=$PWD/src
export PATH=$PATH:$PWD/src
export PATH=$PATH:$PWD/src/pipelines
export PYTHONPATH=$PYTHONPATH:$PWD/src
