import os
import subprocess
from re import sub
from glob import glob
from tqdm                   import tqdm
import multiprocessing      as mp
from functools              import partial
from src.utils import execute_R

def gunzip_file(gzfile, keep = True, outdir = None):
    
    """
    Unzip a compressed file.
    
    Arguments
    ---------
    gzfile: str
        Path to the file to unzip.
    keep: bool
        Option to keep the compressed file after extraction.
    outdir: str
        Path to the output directory. If None, the file will
        be unzipped in its native directory.
        
    Returns
    -------
    outfile: str
        Path to the unzipped file.
    """
    
    os.system('gunzip -f -k {}'.format(gzfile))
    outfile = sub(r'.gz', '', gzfile)
    if not keep:
        os.system('rm {}'.format(gzfile))
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        os.system('mv {} {}'.format(outfile, outdir))
        outfile = os.path.join(outdir, os.path.basename(outfile))

    return outfile


def gunzip_files(infiles, keep = True, outdir = None, parallel = False, nproc = None):
    
    """
    Unzip a set of compressed files.
    
    Arguments
    ---------
    infiles: list
        List of paths to files to unzip.
    keep: bool
        Option to keep the compressed files after extraction.
    outdir: None
        Path to the output directory. If None, the files will
        be unzipped in their native directories.
    parallel: bool
        Option to run in parallel.
    nprpc: int
        Number of processors to use in parallel.
        
    Returns
    -------
    outfiles: list
        List of paths to the unzipped files.
    """
    
    gunzip_file_partial = partial(gunzip_file, 
                                  keep = keep,
                                  outdir = outdir)
    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the number of processors to use")
        pool = mp.Pool(nproc)
        outfiles = []
        for outfile in tqdm(pool.imap(gunzip_file_partial, infiles), total = len(infiles)):
            outfiles.append(outfile)
        pool.close()
        pool.join()
    else:
        outfiles = list(map(gunzip_file_partial, tqdm(infiles)))
        
    return outfiles


def nii_to_mnc(infile, keep = True, outdir = None):
    
    """
    Convert a NIFTY file to MINC format.
    
    Arguments
    ---------
    infile: str
        Path to the NIFTY file to convert.
    outdir: str
        Path to the output directory. If None, the output
        MINC file will be stored in the native directory.
    keep: bool
        Option to keep the input file. 
        
    Returns
    -------
    outfile: str
        Path to the converted MINC file.
    """
    
    os.system('nii2mnc -quiet -clobber {}'.format(infile))
    outfile = sub(r'.nii', '.mnc', infile)
    if not keep:
        os.system('rm {}'.format(infile))
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        os.system('mv {} {}'.format(outfile, outdir))
        outfile = os.path.join(outdir, os.path.basename(outfile))
        
    return outfile


def mnc_to_nii(infile, keep = True, outdir = None):
    
    """
    Convert a MINC file to NIFTY format.
    
    Arguments
    ---------
    infile: str
        Path to the MINC file to convert.
    outdir: str
        Path to the output directory. If None, the output
        MINC file will be stored in the native directory.
    keep: bool
        Option to keep the input file. 
        
    Returns
    -------
    outfile: str
        Path to the converted NIFTY file.
    """
    
    os.system('mnc2nii -quiet {}'.format(infile))
    outfile = sub(r'.mnc', '.nii', infile)
    if not keep:
        os.system('rm {}'.format(infile))
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        os.system('mv {} {}'.format(outfile, outdir))
        outfile = os.path.join(outdir, os.path.basename(outfile))
        
    return outfile
    
    
def convert_images(infiles, input_format = 'nifty', output_format = 'minc', keep = True, outdir = None, parallel = False, nproc = None):
    
    """
    Convert a set of images from NIFTY to MINC or vice versa.
    
    Arguments
    ---------
    infiles: list
        List of paths to images to convert.
    input_format: str
        One of {'nifty', 'minc'} indicating the input format.
    output_format: str
        One of {'nifty', 'minc'} indicating the output format. 
        Must be different from the input format.
    keep: bool
        Option to keep the input file. 
    outdir: str
        Path to the output directory. If None, the converted
        images will be stored in the native directory.
    parallel: bool
        Option to run in parallel.
    nprpc: int
        Number of processors to use in parallel.
        
    Returns
    -------
    outfiles: list
        List of paths to the converted images.
    """
        
    if (input_format == 'nifty') & (output_format == 'minc'):
        converter = partial(nii_to_mnc, keep = keep, outdir = outdir)
    elif (input_format == 'minc') & (output_format == 'nifty'):
        converter = partial(mnc_to_nii, keep = keep, outdir = outdir)
    else:
        raise ValueError 
        
    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the number of processors to use")
        pool = mp.Pool(nproc)
        outfiles = []
        for outfile in tqdm(pool.imap(converter, infiles), total = len(infiles)):
            outfiles.append(outfile)
        pool.close()
        pool.join()
    else:
        outfiles = list(map(converter, tqdm(infiles)))
        
    return outfiles
    
    
def calculate_human_effect_sizes(demographics, imgdir, maskfile, outdir, ncontrols = 10, threshold = 10, dataset = 1, parallel = False, nproc = None):

    """
    Calculate human effect size images.
    
    Arguments
    ---------
    demographics: str
        Path to the CSV file containing the human demographics data.
    imgdir: str
        Path to the directory containing the MINC images to use to compute the effect sizes.
    maskfile: str
        Path to the mask MINC file for the images.
    outdir: str
        Path to the directory in which to save the effect size MINC images.
    ncontrols: int
        Number of propensity-matched controls to use when computing the effect sizes.
    threshold: int
    dataset: int
        A numeric flag indicating which data to use. Set to 1 for all data. Set to 2 for POND and SickKids. Set to 3 for POND only.
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
    
    Returns
    -------
    outfiles: list
        List of paths to the effect size images.
    """

    script = 'calculate_human_effect_sizes.R'
    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the number of processors to use")
    else: 
        nproc = 1
    parallel = 'true' if parallel else 'false'
    script_args = {'demographics':demographics,
                   'imgdir':imgdir,
                   'maskfile':maskfile,
                   'outdir':outdir,
                   'ncontrols':ncontrols,
                   'threshold':threshold,
                   'dataset':dataset,
                   'parallel':parallel,
                   'nproc':nproc}
    execute_R(script = script, args = script_args)
    outfiles = os.listdir(outdir)
    return outfiles


def resample_image(infile, isostep, outdir = None):
    
    """
    Resample a MINC image
    
    Arguments
    ---------
    infile: str
        Path to the image file to resample.
    isostep: float
        Isotropic dimension of voxels in the resampled image (millimeters).
    outdir: str
        Path to the directory in which to save the resampled image. If None, resampled image will be written to the directory of the input file.
    
    
    Returns
    -------
    outfile: str
        Path to the resampled image. 
    """
    
    outfile = sub(r'.mnc', '_resampled_{}.mnc'.format(isostep), infile)
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outfile = os.path.basename(outfile)
        outfile = os.path.join(outdir, outfile)
    cmd_autocrop = ('autocrop -quiet -clobber -isostep {} {} {}'
                    .format(isostep, infile, outfile))
    os.system(cmd_autocrop)

    return outfile


def resample_images(infiles, isostep, outdir = None, parallel = False, nproc = None):
    
    """
    Resample a set of MINC images.
     
    Arguments
    ---------
    infiles: list
        List of paths to images to resample.
    outdir: str
        Path to the output directory. If None, the resampled 
        images will be stored in the native directory.
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
        
    Returns
    -------
    outfiles: list
        List of paths to the resampled images.
    """

    resampler = partial(resample_image,
                        isostep = isostep,
                        outdir = outdir)
    
    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the number of processors to use in parallel.")
        pool = mp.Pool(nproc)
        outfiles = []
        for outfile in tqdm(pool.imap(resampler, infiles), total = len(infiles)):
            outfiles.append(outfile)
        pool.close()
        pool.join()
    else:
        outfiles = list(map(resampler, tqdm(infiles)))
        
    return outfiles
        
