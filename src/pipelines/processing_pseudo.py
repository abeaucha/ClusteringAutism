
pipelines.processing.initialize()
pipelines.processing.effect_sizes()
pipelines.processing.clustering()
pipelines.processing.centroids()
pipelines.processing.main() # ???

def initialize(...):
    return

def effect_sizes(...):
    # This function computes the effect sizes for both absolute and
    # relative Jacobians.
    # Can be executed either locally in series on distributed on a cluster

    if environment is 'local':
        # Iterate over Jacobians
        for j in enumerate(jacobians):
            # Execute the subprocess compute_effect_sizes.R
            execute_local('compute_effect_sizes.R')

    if environment is 'slurm':
        # Deploy a job array with 2 jobs
        # Each job will execute a script compute_effect_sizes.R
        # that computes the effect sizes for one set of Jacobians
        execute_slurm('compute_effect_sizes.R')

def clustering():
    return


def centroids():
    return



def execute_local(script, args):
    # Executes a subprocess for the script with the given
    # command line args
    # 1. Determine if script is Python or R
    # 2. Deploy the appropriate subprocess

def execute_slurm(script, args, args to map over?):
    # Submits a job array executing the script over the parameters
    # How do I submit a job array like this?


# What does compute_effect_sizes.R do?
# Take a set of images and compute effect sizes with the desired method
# Should be able to run propensity matching or normative growth modelling.
# Ideally should be able to distribute the voxel computations across a cluster
# as well.

