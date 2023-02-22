import os
import subprocess
import numpy   as np
import pandas  as pd
from functools import partial
from random    import randint
from re        import sub
from tqdm      import tqdm
from warnings  import warn


def execute_R(script, args):
    
    """
    Execute an R script.
    
    Arguments
    ---------
    script: str
        String containing the name of the R script to execute.
    args: dict
        Dictionary of key-value pairs containing command line arguments 
        to pass to the script.
    
    Returns
    -------
    None
    """
    
    args = [['--'+str(key), str(val)] for key, val in args.items()]
    args = sum(args, [])
    cmd = ['Rscript']+[script]+args
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


def mkdir_from_params(params, basedir = './', metadata = 'metadata.csv', 
                      params_id = None, ndigits = 3, return_metadata = False):

    """
    Create a random directory based on a parameter set.
    
    Arguments
    ---------
    params: dict
        Dictionary of key-value pairs containing the parameter set.
    basedir: str
        Path to the base directory in which to create the directory
        for the parameter set.
    metadata: str
        Name of the CSV file in which to store the parameter set 
        metadata.
    params_id: str
        Option to specify the ID associated with the parameter set.
    ndigits: int
        Number of digits to use when creating the random ID of the 
        parameter directory.
    return_metadata: bool
        Option to return the path to the metadata file.
    
    Returns
    -------
    params_dir: str
        Path to the random directory associated with the parameter set.
    metadata: str (optional)
        Path to the metadata file. 
    """
    
    params_cp = params.copy()
    
    #Create base directory if needed
    basedir = os.path.join(basedir, '')
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    #If metadata file already exists, import it. 
    metadata = os.path.join(basedir, metadata)
    if os.path.exists(metadata):
        
        #Check existence of parameter set in metadata
        df_metadata = pd.read_csv(metadata, dtype = 'str')
        df_match = pd.merge(df_metadata, 
                            pd.DataFrame(params_cp, index = [0], dtype = 'str'),
                            how = 'inner')
        
        #If parameter set doesn't exist, add it to the metadata
        if df_match.shape[0] == 0:

            params_cp['id'] = random_id(ndigits) if params_id is None else params_id
            flag = 1
            while flag == 1:
                id_exists = np.isin(params_cp['id'], df_metadata['id'])
                if id_exists:
                    warn("Parameter ID exists: {}. Randomly generating a new ", 
                         "one.".format(params_id))
                    params_cp['id'] = random_id(ndigits)
                else: 
                    flag = 0
            df_metadata = pd.concat([df_metadata, 
                                     pd.DataFrame(params_cp, index = [0])])
            df_metadata.to_csv(metadata, index = False)
            
        elif df_match.shape[0] == 1:
            params_cp['id'] = df_match['id'].values[0]
            
        else:
            raise Exception
            
    else: 
        #If metadata file doesn't exist, create it with the parameter set
        params_cp['id'] = random_id(ndigits) if params_id is None else params_id
        pd.DataFrame(params_cp, index = [0]).to_csv(metadata, index = False)

    #Create parameter directory if needed
    params_dir = os.path.join(basedir, params_cp['id'], '')
    if not os.path.exists(params_dir):
        os.makedirs(params_dir)

    #Return path to parameter directory and metadata file as desired
    if return_metadata: 
        return params_dir, metadata
    else:
        return params_dir

    
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
    
    #Create output dir if needed
    dst = os.path.join(dst, '')
    if not os.path.exists(dst):
        os.makedirs(dst)
    
    #Create symlinks
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
    
    num = randint(0, (10**n)-1)
    fstring = '{:0'+str(n)+'}'
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
            raise ValueError("Set the nproc argument to specify the number ",
                             "of processors to use")
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
        Path to the output directory. If None, the output MINC file 
        will be stored in the native directory.
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
        Path to the output directory. If None, the output MINC file 
        will be stored in the native directory.
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
            raise ValueError("Set the nproc argument to specify the number ",
                             "of processors to use")
        pool = mp.Pool(nproc)
        outfiles = []
        for outfile in tqdm(pool.imap(converter, infiles), total = len(infiles)):
            outfiles.append(outfile)
        pool.close()
        pool.join()
    else:
        outfiles = list(map(converter, tqdm(infiles)))
        
    return outfiles