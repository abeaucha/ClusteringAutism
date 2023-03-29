#!/bin/bash

module use /project/def-jlerch/tools/modulefiles/
module load mice-env

if [ ! -d ".ccenv" ]; then
	echo "Initializing python virtual environment..."

	python3 -m venv .ccenv

	source .ccenv/bin/activate

fi
