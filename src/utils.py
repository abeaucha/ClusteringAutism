import os
import subprocess
import multiprocessing as mp
import pandas as pd
from functools import partial
from random import randint
from re import sub
from tqdm import tqdm


def execute_R(script, args):
    """
    Execute an R script.
    
    Arguments
    ---------
    :param script: str
        String containing the name of the R script to execute.
    :param args: dict
        Dictionary of key-value pairs containing command line arguments 
        to pass to the script.
    
    Returns
    -------
    :return: None
    """

    args = [['--' + str(key), str(val)] for key, val in args.items()]
    args = sum(args, [])
    cmd = ['Rscript'] + [script] + args
    subprocess.run(cmd)
    return


def mkdir_from_list(inlist, basedir = './', sep = '_'):
    """
    Create a directory from a list of strings.
    
    Arguments
    ---------
    inlist: list
        List of strings to use in directory name.
    basedir: str
        Path to the base directory in which to create the new directory.
    sep: str
        Separator used to concatenate the strings in `inlist`.
    
    Returns
    -------
    outdir: str
        Path to the new directory.
    """

    inlist = sep.join(inlist)
    outdir = os.path.join(basedir, inlist)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir


def get_params_id(params, metadata):
    df_params = pd.DataFrame(params, index = [0], dtype = str)
    if os.path.exists(metadata):
        df_metadata = pd.read_csv(metadata, dtype = str)
        df_match = pd.merge(df_metadata, df_params, how = 'inner')
        nmatch = df_match.shape[0]
        if nmatch == 1:
            return df_match['id'].values[0]
        else:
            return None
    else:
        return None

def fetch_metadata(metadata, **kwargs):

    if not os.path.exists(metadata):
        raise OSError("Metadata file not found: {}".format(metadata))

    df_metadata = pd.read_csv(metadata, dtype = str)
    if any(kwargs.keys()):
        df_params = pd.DataFrame(kwargs, index = [0], dtype = str)
    else:
        print("No parameters provided. Returning all metadata.")
        return df_metadata

    df_match = pd.merge(df_metadata, df_params, how = 'inner')
    nmatch = df_match.shape[0]
    if nmatch == 0:
        return None
    else:
        return df_match

def set_params_id(params, metadata, params_id = None):
    df_params = pd.DataFrame(params, index = [0], dtype = str)
    if os.path.exists(metadata):
        df_metadata = pd.read_csv(metadata, dtype = str)
        df_match = pd.merge(df_metadata, df_params, how = 'inner')
        nmatch = df_match.shape[0]
        if nmatch == 1:
            df_params['id'] = df_match['id'].values[0]
        elif nmatch == 0:
            df_params['id'] = (random_id(3) if params_id is None
                               else str(params_id))
            df_metadata = pd.concat([df_metadata, df_params])
            df_metadata.to_csv(metadata, index = False)
        else:
            raise Exception
    else:
        df_params['id'] = random_id(3) if params_id is None else str(params_id)
        df_params.to_csv(metadata, index = False)

    return df_params['id'].values[0]


def mkdir_from_params(params, outdir, metadata = None, params_id = None):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if metadata is None:
        metadata = os.path.join(outdir, 'metadata.csv')

    params_id_check = get_params_id(params = params,
                                    metadata = metadata)

    if params_id_check is None:
        params_id = set_params_id(params = params,
                                  metadata = metadata,
                                  params_id = params_id)
    else:
        if params_id is not None:
            print("Parameters already identified: {}".format(params_id_check))
        params_id = params_id_check

    outdir = os.path.join(outdir, params_id, '')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    return outdir


def mk_symlinks(src, dst):
    """
    Create a set of symbolic links.
    
    Arguments
    ---------
    src: list
        List of paths to the source files.
    dst: str
        Path to the directory in which to save the links.
    
    Returns
    -------
    dstfiles: list
        List of paths to the symbolic links. 
    """

    # Create output dir if needed
    dst = os.path.join(dst, '')
    if not os.path.exists(dst):
        os.makedirs(dst)

    # Create symlinks
    dstfiles = []
    for i, srcfile in enumerate(src):
        dstfile = os.path.join(dst, os.path.basename(srcfile))
        if not os.path.exists(dstfile):
            os.symlink(os.path.relpath(srcfile, dst), dstfile)
        dstfiles.append(dstfile)

    return dstfiles


def random_id(n):
    """
    Create a random ID with a specific number of digits.
    
    Arguments
    ---------
    n: int
        Number of digits in the random ID.
        
    Returns
    -------
    randid: str
        Random ID as string.
    """

    num = randint(0, (10 ** n) - 1)
    fstring = '{:0' + str(n) + '}'
    randid = fstring.format(num)
    return randid


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

    subprocess.run(['gunzip', '-f', '-k', gzfile])
    outfile = sub(r'.gz', '', gzfile)
    if not keep:
        subprocess.run(['rm', gzfile])
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        subprocess.run(['mv', outfile, outdir])
        outfile = os.path.join(outdir, os.path.basename(outfile))

    return outfile


def gunzip_files(infiles, keep = True, outdir = None,
                 parallel = False, nproc = None):
    """
    Unzip a set of compressed files.
    
    Arguments
    ---------
    infiles: list
        List of paths to files to unzip.
    keep: bool
        Option to keep the compressed files after extraction.
    outdir: None
        Path to the output directory. If None, the files will be 
        unzipped in their native directories.
    parallel: bool
        Option to run in parallel.
    nproc: int
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
            raise ValueError("Set the nproc argument to specify the number ",
                             "of processors to use")
        pool = mp.Pool(nproc)
        outfiles = []
        for outfile in tqdm(pool.imap(gunzip_file_partial, infiles),
                            total = len(infiles)):
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
        Path to the output directory. If None, the output MINC file 
        will be stored in the native directory.
    keep: bool
        Option to keep the input file. 
        
    Returns
    -------
    outfile: str
        Path to the converted MINC file.
    """

    subprocess.run(['nii2mnc', '-quiet', 'clobber', infile])
    outfile = sub(r'.nii', '.mnc', infile)
    if not keep:
        subprocess.run(['rm', infile])
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        subprocess.run(['mv', outfile, outdir])
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
        Path to the output directory. If None, the output MINC file 
        will be stored in the native directory.
    keep: bool
        Option to keep the input file. 
        
    Returns
    -------
    outfile: str
        Path to the converted NIFTY file.
    """

    subprocess.run(['mnc2nii', '-quiet', infile])
    outfile = sub(r'.mnc', '.nii', infile)
    if not keep:
        subprocess.run(['rm', infile])
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        subprocess.run(['mv', outfile, outdir])
        outfile = os.path.join(outdir, os.path.basename(outfile))

    return outfile


def convert_images(infiles, input_format = 'nifty', output_format = 'minc',
                   keep = True, outdir = None, parallel = False, nproc = None):
    """
    Convert a set of images between NIFTY and MINC formats.
    
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
    nproc: int
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
            raise ValueError("Set the nproc argument to specify the number ",
                             "of processors to use")
        pool = mp.Pool(nproc)
        outfiles = []
        for outfile in tqdm(pool.imap(converter, infiles),
                            total = len(infiles)):
            outfiles.append(outfile)
        pool.close()
        pool.join()
    else:
        outfiles = list(map(converter, tqdm(infiles)))

    return outfiles


def transform_image(infile, like, transform, outdir = None, suffix = None):
    # Append suffix if specified
    if suffix is not None:
        outfile = (os.path.splitext(infile)[0] +
                   suffix +
                   os.path.splitext(infile)[1])
    else:
        outfile = infile

    # Create output directory if needed
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outfile = os.path.basename(outfile)
        outfile = os.path.join(outdir, outfile)
    else:
        if suffix is None:
            raise Exception("Arguments outdir and suffix "
                            "cannot both be None.")

    # Command line mincresample command
    cmd_mincresample = ['mincresample', '-quiet', '-clobber',
                        '-like', like, '-transform', transform,
                        infile, outfile]
    log = subprocess.run(cmd_mincresample, stdout = subprocess.PIPE)

    return outfile


def transform_images(infiles, like, transform, outdir = None, suffix = None,
                     parallel = False, nproc = None):
    # Create output directory if needed
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    # Partial transforming function
    transformer = partial(transform_image,
                          like = like,
                          transform = transform,
                          outdir = outdir,
                          suffix = suffix)

    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the "
                             "number of processors to use in parallel.")
        pool = mp.Pool(nproc)
        outfiles = []
        for outfile in tqdm(pool.imap(transformer, infiles),
                            total = len(infiles)):
            outfiles.append(outfile)
        pool.close()
        pool.join()
    else:
        outfiles = list(map(transformer, tqdm(infiles)))

    return outfiles


def resample_image(infile, isostep, outdir = None, suffix = None):
    """
    Resample a MINC image
    
    Arguments
    ---------
    infile: str
        Path to the image file to resample.
    isostep: float
        Isotropic dimension of voxels in the resampled image (mm).
    outdir: str
        Path to the directory in which to save the resampled image. 
        If None, resampled image will be written to the directory 
        of the input file.
    suffix: str
        Suffix to append to output filename before the file extension.
    
    Returns
    -------
    outfile: str
        Path to the resampled image. 
    """

    # Append suffix if specified
    if suffix is not None:
        outfile = (os.path.splitext(infile)[0] +
                   suffix +
                   os.path.splitext(infile)[1])
    else:
        outfile = infile

    # Create output directory if needed
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outfile = os.path.basename(outfile)
        outfile = os.path.join(outdir, outfile)
    else:
        if suffix is None:
            raise Exception("Arguments outdir and suffix "
                            "cannot both be None.")

    # Autocrop command
    cmd_autocrop = ['autocrop', '-quiet', '-clobber',
                    '-isostep', isostep, infile, outfile]
    log = subprocess.run(cmd_autocrop, stdout = subprocess.PIPE)

    return outfile


def resample_images(infiles, isostep, outdir = None, suffix = None,
                    parallel = False, nproc = None):
    """
    Resample a set of MINC images.
     
    Arguments
    ---------
    infiles: list
        List of paths to images to resample.
    isostep: float
        Isotropic dimension of voxels in the resampled image (mm).
    outdir: str
        Path to the output directory. If None, the resampled 
        images will be stored in the native directory.
    suffix: str
        Suffix to append to output filename before the file extension.
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
        
    Returns
    -------
    outfiles: list
        List of paths to the resampled images.
    """

    # Create output directory if needed
    if outdir is not None:
        outdir = os.path.join(outdir, '')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

    # Partial resampling function
    resampler = partial(resample_image,
                        isostep = isostep,
                        outdir = outdir,
                        suffix = suffix)

    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the "
                             "number of processors to use in parallel.")
        pool = mp.Pool(nproc)
        outfiles = []
        for outfile in tqdm(pool.imap(resampler, infiles),
                            total = len(infiles)):
            outfiles.append(outfile)
        pool.close()
        pool.join()
    else:
        outfiles = list(map(resampler, tqdm(infiles)))

    return outfiles
