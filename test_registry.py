#!/usr/bin/env python3

import pandas as pd
import utils
import os

resources = dict(
    nodes = 1,
    mem = "8G",
    time = "00:10:00"
)

registry = utils.Registry(resources = resources)


registry.create_batches(x = pd.DataFrame([range(100), range(100)]),
                        nbatches = 2)


registry.create_jobs(script = 'sleep 5m')

import subprocess

output = subprocess.run(['sbatch', registry.jobs[0]], capture_output=True)

# To get job status
# squeue --me --name=tmpiqg_by6t --format="%t" --noheader

job = registry.jobs[0]
jobname = os.path.basename(job)
jobname = os.path.splitext(jobname)[0]
