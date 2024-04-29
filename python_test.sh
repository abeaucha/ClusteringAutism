#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=test_%j.out

echo "PATH before venv:"
echo $PATH
source activate_venv_hpc.sh
echo "PATH after venv:"
echo $PATH
python_test.py