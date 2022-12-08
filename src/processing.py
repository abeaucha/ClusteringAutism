import os
import subprocess
from re import sub
from glob import glob
from tqdm                   import tqdm
import multiprocessing      as mp
from functools              import partial
from src.utils import execute_R
from pyminc.volumes.factory import volumeFromFile
import numpy as np
from warnings import warn
import pandas as pd

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
        

def import_image(img, mask = None, flatten = True):
    
    """
    Import a MINC image.
    
    Arguments
    ---------
    img: str
        Path to the MINC image to import.
    mask: str
        Optional path to the MINC mask. 
    flatten: bool
        Option to flatten image into a 1-dimensional array. 
        If True and a mask is provided, only the voxels in the mask will be returned.
    
    Returns
    -------
    img: numpy.ndarray
        A NumPy array containing the (masked) image.
    """
    
    img_vol = volumeFromFile(img)
    img_dims = img_vol.getSizes()
    img_seps = img_vol.getSeparations()
    img = np.array(img_vol.data)
    img_vol.closeVolume()

    if flatten:
        img = img.flatten()

    if mask is not None:

        mask_vol = volumeFromFile(mask)
        mask_dims = mask_vol.getSizes()
        mask_seps = mask_vol.getSeparations()
        mask = np.array(mask_vol.data)
        mask_vol.closeVolume()

        if mask_seps != img_seps:
            raise Exception("Input image and mask have different resolutions.")
        if mask_dims != img_dims:
            raise Exception("Input image and mask have different dimensions.")

        if flatten:
            mask = mask.flatten()
            img = img[mask == 1]
        else:
            img[mask == 0] = 0
    
    return img


def import_images(infiles, mask = None, output_format = 'list', flatten = True, parallel = False, nproc = None):

    """
    Import a set of MINC images.
    
    Arguments
    ---------
    infiles: list
        List of paths to images to import.
    mask: str
        Optional path to a mask image. 
    output_format: str
        One of {'list', 'numpy', 'pandas'} indicating what format to return.
    flatten: bool
        Option to flatten images into a 1-dimensional array. 
        If True and a mask is provided, only the voxels in the mask will be returned.
        If False and output_format is not 'list', images will be flattened regardless.
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
    
    Returns
    -------
    imgs
        A list, NumPy array, or Pandas DataFrame containing the (masked) images.
    """

    format_opts = ['list', 'numpy', 'pandas']
    format_test = sum([output_format == opt for opt in format_opts])
    format_err = ("Argument output_format must be one of {}: {}"
                   .format(format_opts, output_format))
    if format_test != 1:
        raise ValueError(format_err)

    if not flatten:
        if output_format != 'list':
            msg_warn = ("flatten = False is only valid when output_format = 'list'. "
                        "Proceeding with flattened images.")
            warn(msg_warn)
            flatten = True

    importer = partial(import_image,
                       mask = mask,
                       flatten = flatten)

    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the number of processors to use in parallel.")
        pool = mp.Pool(nproc)
        imgs = []
        for img in tqdm(pool.imap(importer, infiles), total = len(infiles)):
            imgs.append(img)
        pool.close()
        pool.join()
    else:
        imgs = list(map(importer, tqdm(infiles)))

    imgsize_test = [len(img) for img in imgs]
    imgsize_err = "Images provided contain different numbers of voxels."
    if len(set(imgsize_test)) != 1:
        raise Exception(imgsize_err)

    if output_format == 'numpy':
        return np.asarray(imgs)
    elif output_format == 'pandas':
        return pd.DataFrame(np.asarray(imgs))
    else: 
        return imgs
    
    

def build_voxel_matrix(infiles, mask = None, file_col = False, sort = False, save = False, outfile = 'voxel_matrix.csv', parallel = False, nproc = None):
    
    """
    Create a data frame of voxels from a set of images.
    
    Arguments
    ---------
    infiles: list
        List of paths to images.
    mask: str
        Optional path to a mask image. 
    file_col: bool
        Option to store input files in a column. If true, the paths in infiles are stored in a column called 'file'.
    sort: bool
        Option to sort rows based on file names.
    save: bool
        Option to save to CSV.
    outfile: str
        Path to the output CSV file. Ignored if save = False.
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
    
    Returns
    -------
    df_imgs
        A Pandas DataFrame containing the (masked) images.
    """
    
    df_imgs = import_images(infiles = infiles,
                           mask = mask,
                           output_format = 'pandas',
                           flatten = True,
                           parallel = parallel,
                           nproc = nproc)
    df_imgs['file'] = infiles

    if sort:
        df_imgs = df_imgs.sort_values(by = 'file')

    if not file_col:
        df_imgs = df_imgs.drop('file', axis = 1)
    
    if save:
        df_imgs.to_csv(outfile, index = False)
    
    return df_imgs


def cluster_human_data(infiles, rownames = None, nk_max = 10, metric = 'correlation', K = 10, sigma = 0.5, t = 20, outfile = 'clusters.csv'):

    """
    Cluster human effect size images using similarity network fusion.
    
    Arguments
    ---------
    infiles: list
        List of paths to the CSV files containing effect size data.
    rownames: str
        Column in CSV files containing row names.
    nk_max: int
        Maximum number of clusters to identify. Solutions will be obtained for nk = 2 to nk = nk_max.
    metric: str
        Distance metric used to compute the SNF affinity matrices.
    K: int
        Number of nearest-neighbours used to compute the SNF affinity matrices.
    sigma: float
        Variance for the local model in the SNF affinity matrices.
    t: int
        Number of iterations for the diffusion process in SNF.
    outfile: str
        Path to the CSV file in which to write the cluster assignments.
        
    Returns
    -------
    A Pandas DataFrame containing cluster assignments.
    """
    
    # TO-DO
    # - Make this function more general so that it can run SNF on any set of voxel matrices
    #   rather than just the specific human effect sizes. 
    # - Option to run SNF on more than two matrices? 
    # Make it so that this function exits if an error happens in the R script
    
    
    if type(infiles) is not list:
        raise ValueError("Argument infiles must be a list.")
    
    if len(infiles) != 2:
        raise Exception("Argument infiles must have length 2.")
    
    script = 'cluster_human_data.R'
    script_args = {'file1':infiles[0],
                   'file2':infiles[1],
                   'rownames':rownames,
                   'nclusters':nk_max,
                   'metric':metric,
                   'K':K,
                   'sigma':sigma,
                   't':t,
                   'outfile':outfile}
    
    if rownames is None:
        del script_args['rownames']
    
    execute_R(script = script, args = script_args)
    
    return pd.read_csv(outfile)