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

registry.create_batches(x = pd.DataFrame(dict(x=range(5832),y=range(5832))),
                        nbatches = 400)

registry.create_jobs(script = 'sleep 60m')

registry.submit_jobs(wait = True, cleanup = False)


cmd = ('sacct --jobs={} --name={} --format=JobID%20,JobName%{},STATE%20'
        .format(','.join(registry.jobids), ','.join(registry.jobnames), 
                len(registry.jobnames[0])))

output = (subprocess.run(cmd.split(' '), capture_output = True)
            .stdout.decode('UTF-8').splitlines())

output = [out.split() for out in output]
keys = output[0]
rows = [r for r in output[2:] if r[0] in registry.jobids and r[1] in registry.jobnames]

test = [r for r in output[2:] if r[0] in registry.jobids]
test = [r for r in output[2:] if r[1] not in registry.jobnames]


import pandas as pd
import utils
import os
from importlib import reload 
reload(utils)
name = 'compute_cluster_similarity_registry_002'
registry = utils.Registry(name = name, attach = True)
