
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

    # This function also has to generate the effect size matrices.
    # Is that separate from compute_effect_sizes.R?

def clustering():
    # This function clusters the absolute and relative effect sizes
    # This is the part I'm unsure about. If I run a job array for the
    # effect sizes, then I can submit this as a job with dependencies.
    if environment is 'local':
        execute_local('clustering.R')
    if environment is 'slurm':
        execute_slurm('clustering.R')
        # This needs to get the jobs from the effect_sizes() module to use
        # as dependencies
    return


def centroids():
    # This function generates the cluster centroid images
    if environment is 'local':
        execute_local('create_centroids.R')
    if environment is 'slurm':
        execute_slurm('create_centroids.R')
        # This needs to get the job from the clustering() module to use
        # as a dependency
    return



def execute_local(script, args):
    # Executes a subprocess for the script with the given
    # command line args
    # 1. Determine if script is Python or R
    # 2. Deploy the appropriate subprocess
    # If I make the scripts executable then I can just source them, no?
    # I don't have to specify the interpreter

def execute_slurm(script, args, args to map over?):
    # Submits a job array executing the script over the parameters
    # How do I submit a job array like this?
    # This also needs to be able to submit a single job, rather than just a
    # job array


# What does compute_effect_sizes.R do?
# Take a set of images and compute effect sizes with the desired method
# Should be able to run propensity matching or normative growth modelling.
# Ideally should be able to distribute the voxel computations across a cluster
# as well.

