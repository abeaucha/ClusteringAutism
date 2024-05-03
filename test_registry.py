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
    time = "00:60:00"
)

registry = utils.Registry(resources = resources)

registry.create_batches(x = pd.DataFrame([range(100), range(100)]),
                        nbatches = 2)

registry.create_jobs(script = 'sleep 60m')

registry.submit_jobs()



# for job in zip(registry.jobids, registry.jobnames):
# cmd = 'sacct --jobs={} --name={} --format=JobID%8,JobName%11,STATE --noheader'.format(','.join(registry.jobids), ','.join(registry.jobnames))

sacct --jobs=9971796,9971797 --name=registry_002_batch_0__vncqagco,registry_002_batch_1__fonns_vd --format=JobID%8,JobName%30,STATE

cmd = ('sacct --jobs={} --name={} --format=JobID%8,JobName{},STATE'
               .format(','.join(registry.jobids), ','.join(registry.jobnames), 
                       len(registry.jobnames[0])))
output = (subprocess.run(cmd.split(' '), capture_output = True)
            .stdout.decode('UTF-8').splitlines())
output = [out.split() for out in output]
keys = output[0]
rows = [r for r in output[2:] 
        if r[0] in registry.jobids and r[1] in registry.jobnames]
states = {x[1]:[r[x[0]] for r in rows] for x in enumerate(keys)}


wait = True
# while wait:

# def _check_status(self):


status = registry._fetch_status()

codes_active = ['RUNNING', 'PENDING']
n_active = 0
for code in codes_active:
    n_active += status['State'].count(code)
completed = True if n_active == 0 else False
return completed






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
