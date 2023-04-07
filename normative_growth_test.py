
import processing

imgdir = 'data/human/derivatives/POND_SickKids/jacobians/resolution_3.0/absolute/'
demographics = 'data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc_wfile.csv'
mask = 'data/human/registration/v1/reference_files/mask_3.0mm.mnc'
outdir = 'data/human/test/'
key = 'file'
df = 3
batch = ['Site', 'Scanner']
nbatches = 2
nproc = 8

tmp = processing.normative_growth_norm(imgdir = imgdir,
                     demographics = demographics,
                     mask = mask,
                     outdir = outdir,
                     key = key,
                     df = df,
                     batch = batch,
                     nbatches = nbatches,
                     nproc = nproc)