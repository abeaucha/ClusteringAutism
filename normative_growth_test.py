
import processing

imgdir = 'data/human/registration/v1/jacobians_resampled/resolution_0.8/absolute/'
demographics = 'data/human/derivatives/POND_SickKids/DBM_input_demo_passedqc_wfile.csv'
mask = 'data/human/registration/v1/reference_files/mask_0.8mm.mnc'
outdir = '/scratch/abeaucha/data/human/test/normative_growth/'
key = 'file'
df = 3
batch = ['Site', 'Scanner']
nbatches = 4
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
