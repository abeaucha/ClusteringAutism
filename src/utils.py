import subprocess
import os
import pandas as pd
import numpy as np
from warnings import warn
from random import randint

def execute_R(script, args):
    
    """
    Execute an R script.
    
    Arguments
    ---------
    script: str
        String containing the name of the R script to execute.
    args: dict
        Dictionary of key-value pairs containing command line arguments to pass to the script
    
    Returns
    -------
    None
    """
    
    args = [['--'+str(key), str(val)] for key, val in args.items()]
    args = sum(args, [])
    cmd = ['Rscript']+[script]+args
    subprocess.run(cmd)
    return


def random_id(n):
    
    """
    Create a random ID with a specific number of digits
    
    Arguments
    ---------
    n: int
        Number of digits in the random ID
        
    Returns
    -------
    randid: str
        Random ID as string
    """
    
    num = randint(0, (10**n)-1)
    fstring = '{:0'+str(n)+'}'
    randid = fstring.format(num)
    return randid


def mkdir_from_params(params, basedir = './', metadata = 'metadata.csv', params_id = None, ndigits = 3, return_metadata = False):

    """
    Create a random directory based on a parameter set.
    
    Arguments
    ---------
    params: dict
    basedir: str
    metadata: str
    params_id: str
    ndigits: int
    return_metadata: bool
    
    
    Returns
    -------
    params_dir: str
    metadata: str (optional)
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
                    warn('Parameter ID exists: {}. Randomly generating a new one.'.format(params_id))
                    params_cp['id'] = random_id(ndigits)
                else: 
                    flag = 0
            df_metadata = pd.concat([df_metadata, pd.DataFrame(params_cp, index = [0])])
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

    
def mkdir_from_list(inlist, basedir = './', sep = '_'):
    
    """
    
    Arguments
    ---------
    inlist: list
    basedir: str
    sep: str
    
    Returns
    -------
    outdir: str
    
    """
    
    inlist = sep.join(inlist)
    outdir = os.path.join(basedir, inlist)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir


def mk_symlinks(src, dst):

    """
    Create a set of symbolic links
    
    Arguments
    ---------
    src: list
    dst: str
    
    Returns
    -------
    dstfiles: list
    
    """
    
    dst = os.path.join(dst, '')
    if not os.path.exists(dst):
        os.makedirs(dst)
    
    dstfiles = []
    for i, srcfile in enumerate(src):
        dstfile = os.path.join(dst, os.path.basename(srcfile))
        if not os.path.exists(dstfile):
            os.symlink(os.path.relpath(srcfile, dst), dstfile)
        dstfiles.append(dstfile)
            
    return dstfiles