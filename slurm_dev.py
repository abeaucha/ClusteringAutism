import processing
import os
import pandas as pd
from glob import glob
import utils

if __name__ == '__main__':

    pipeline_dir = 'data/tmp/'
    input_dir = 'data/human/registration/jacobians_resampled/'
    resolution = 3.0
    demographics = 'data/human/registration/DBM_input_demo_passedqc_wfile.csv'
    datasets = ['POND', 'SickKids']
    mask = 'data/human/registration/reference_files/mask_3.0mm.mnc'
    es_matrix_file = 'effect_sizes.csv'
    parallel = True
    nproc = 8
    es_method = 'normative-growth'
    es_df = 3
    es_combat = True
    es_combat_batch = ['Site', 'Scanner']
    es_ncontrols = None


    df_demographics = pd.read_csv(demographics)

    df_demographics = (df_demographics
                       .loc[df_demographics['Dataset'].isin(datasets)]
                       .copy())

    demographics = os.path.join(pipeline_dir, os.path.basename(demographics))
    df_demographics.to_csv(demographics, index = False)

    input_dir = os.path.join(input_dir,
                             'resolution_{}'.format(resolution),
                             '')

    # Paths to pipeline image directory
    imgdir = os.path.join(pipeline_dir,
                          'jacobians',
                          'resolution_{}'.format(resolution),
                          '')

    es_dir = os.path.join(pipeline_dir, 'effect_sizes')

    jac = 'absolute'
    input_files = glob(os.path.join(input_dir, jac, '') + '*.mnc')
    input_files_in_dataset = [[f for f in input_files if g in f][0] for g in
                              df_demographics['file'].to_list()]
    imgfiles = utils.mk_symlinks(src = input_files_in_dataset,
                                 dst = os.path.join(imgdir, jac, ''))

    es_kwargs = dict(imgdir = os.path.join(imgdir, jac, ''),
                     demographics = demographics,
                     mask = mask,
                     outdir = os.path.join(es_dir, jac, ''),
                     matrix_file = es_matrix_file,
                     parallel = parallel,
                     nproc = nproc,
                     method = es_method)

    if es_method == 'normative-growth':
        es_kwargs.update(
            dict(df = es_df,
                 combat = es_combat,
                 combat_batch = es_combat_batch)
        )
    else:
        es_kwargs.update(
            dict(ncontrols = es_ncontrols)
        )

    es_files = processing.calculate_human_effect_sizes(**es_kwargs)