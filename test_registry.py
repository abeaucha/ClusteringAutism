import pandas as pd
import utils
import os
import subprocess
from importlib import reload 
from time import sleep

reload(utils)

resources = dict(
    nodes = 1,
    mem = "8G",
    time = "00:10:00"
)

registry = utils.Registry(resources = resources)

registry.create_batches(x = pd.DataFrame([range(100), range(100)]),
                        nbatches = 2)

registry.create_jobs(script = 'sleep 60m')

registry.submit_jobs()

cmd = 'sacct --jobs=9958837 --format=JobID,STATE --noheader'
tmp = subprocess.run(cmd.split(' '), capture_output=True)
tmp = tmp.stdout.decode('UTF-8')
tmp = tmp.splitlines()


# To get job status
# squeue --me --name=tmpiqg_by6t --format="%t" --noheader

cmds = ['squeue --me --name={} --format=%i --noheader'.format(job) for job in registry.jobnames]
output = [subprocess.run(cmd.split(' '),capture_output=True) for cmd in cmds]
output = [out.stdout.decode('UTF-8').replace('\n', '') for out in output]

# To get job ID
# squeue --me --name=tmpiqg_by6t --format="%i" --noheader


job = registry.jobs[0]
jobname = os.path.basename(job)
jobname = os.path.splitext(jobname)[0]

registry.jobnames = [os.path.basename(job) for job in registry.jobs]
registry.jobnames = [os.path.splitext(job)[0] for job in registry.jobnames]
