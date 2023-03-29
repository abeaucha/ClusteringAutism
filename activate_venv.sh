#!/bin/bash
#Project environment on Graham

module purge
module use /project/def-jlerch/tools/modulefiles/
module load python/3.10.2

if [ ! -d ".venv" ]; then
	echo "Initializing python virtual environment..."
	python3 -m venv .venv
	source .venv/bin/activate
	pip install pip --upgrade
	pip3 install -r python_reqs.txt
	deactivate
fi

module load mice-env
source .venv/bin/activate
export PYTHONPATH="$PWD/src"
