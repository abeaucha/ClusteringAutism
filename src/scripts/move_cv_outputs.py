

import os
from shutil import copytree, rmtree

path_hpf = '/hpf/largeprojects/MICe/abeauchamp/Paper_ClusteringAutism/main/data/human/derivatives/v3/916/cross_validation/'

path_zfs = '/projects/abeauchamp/Projects/MouseHumanMapping/Paper_ClusteringAutism/main/data/human/derivatives/v3/916/cross_validation'

start = 50
end = 70

for i in range(start, end+1):
    print("Copying sample {}".format(i))
    subdir = 'sample_{}'.format(i)
    copytree(os.path.join(path_hpf, subdir, ''), os.path.join(path_zfs, subdir, ''), dirs_exist_ok = True)

for i in range(start, end+1):
    print("Removing sample {}".format(i))
    subdir = 'sample_{}'.format(i)
    rmtree(os.path.join(path_hpf, subdir, ''))
