#!/bin/bash
#SBATCH --job-name=metric_uncertainty_POND
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --chdir=/scratch/abeaucha/ClusteringAutism/main/
#SBATCH --output=metric_uncertainty_POND_%j.out
#SBATCH --error=metric_uncertainty_POND_%j.err

# Activate virtual environment
source activate_venv.sh

# (Optional but usually good hygiene)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

B=1000
# B=96
NPROCS=24
DATASET=POND

# echo "Launching job 1..."
# srun --exclusive -N1 -n1 -c 8 Rscript figures/v3/evaluate_cluster_metric_uncertainty_copy.R 1 8 8 ${DATASET} figures/v3/resources/chunk_1.csv &

# echo "Launching job 2..."
# srun --exclusive -N1 -n1 -c 8 Rscript figures/v3/evaluate_cluster_metric_uncertainty_copy.R 9 16 8 ${DATASET} figures/v3/resources/chunk_2.csv &


# Launch 24 job steps, each step gets 1 task (=1 R process) and 8 CPUs
for t in $(seq 1 $NPROCS); do
  start=$(( (t-1)*B/NPROCS + 1 ))
  end=$(( t*B/NPROCS ))

  echo "Launching task $t for iterations $start to $end"

  srun --overlap -N1 -n1 -c ${SLURM_CPUS_PER_TASK} \
    Rscript figures/v3/evaluate_cluster_metric_uncertainty.R ${start} ${end} ${SLURM_CPUS_PER_TASK} ${DATASET} figures/v3/resources/cluster_metrics_uncertainty_${DATASET}_${t}.csv &
done

wait


deactivate
