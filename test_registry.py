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

registry.submit_jobs(wait = False, cleanup = False)


batches = [os.path.join(registry.paths['batches'], batch) for batch in os.listdir(registry.paths['batches'])]

jobs = [os.path.join(registry.paths['jobs'], job) for job in os.listdir(registry.paths['jobs'])]

[os.path.splitext(os.path.basename(job))[0] for job in registry.jobs]




import pandas as pd
import utils
import os
from importlib import reload 
reload(utils)
name = 'compute_cluster_similarity_registry_002'
registry = utils.Registry(name = name, attach = True)
