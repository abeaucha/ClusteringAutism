#!/bin/bash

module purge
module use /project/def-jlerch/tools/modulefiles/
module load python/3.10.2

if [ ! -d ".ccenv" ]; then
	echo "Initializing python virtual environment..."
	python3 -m venv .ccenv
	source .ccenv/bin/activate
	pip install pip --upgrade
	pip3 install -r python_reqs_cc.txt
	deactivate
fi

module load mice-env
source .ccenv/bin/activate
export PYTHONPATH="$PWD/src"
