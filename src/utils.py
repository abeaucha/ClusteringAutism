import math
import os
import shutil
import subprocess
import sys
import multiprocessing as mp
import pandas as pd
import tempfile
from contextlib import contextmanager
from functools import partial
from functools import wraps
from random import randint
from re import sub
from time import time, perf_counter, sleep
from tqdm import tqdm
from warnings import warn


# Classes ---------------------------------------------------------------------


class Registry:

    def _create_registry(self):

        if self.resources is None:
            raise ValueError("Argument 'resources' must be specified "
                             "when creating a registry.")        
        
        i = 1; iterate = True
        while iterate:
            id = f'{i:03d}'
            name_test = '{}_{}'.format(self.name, id)
            if os.path.exists(name_test):
                i += 1
            else:
                self.name = name_test
                iterate = False

        if self._verbose:
            print("Creating registry: {}".format(self.name), flush = True)

        self.paths = dict(
            jobs = os.path.join(self.name, 'jobs'),
            logs = os.path.join(self.name, 'logs'),
            batches = os.path.join(self.name, 'batches'),
            outputs = os.path.join(self.name, 'outputs')
        )
        self.batches = []
        self.jobs = []
        self.jobnames = []
        self.outputs = []

        for path in self.paths.values():
            if not os.path.exists(path):
                os.makedirs(path)

        return

    
    def _attach_registry(self):

        if not os.path.exists(self.name):
            raise OSError("Registry not found: {}".format(self.name))

        if self._verbose:
            print("Attaching registry: {}".format(self.name), flush = True)

        self.paths = dict(
            jobs = os.path.join(self.name, 'jobs'),
            logs = os.path.join(self.name, 'logs'),
            batches = os.path.join(self.name, 'batches'),
            outputs = os.path.join(self.name, 'outputs')
        )
        self.batches = [os.path.join(self.paths['batches'], batch) 
                        for batch in os.listdir(self.paths['batches'])]
        self.jobs = [os.path.join(self.paths['jobs'], job) 
                        for job in os.listdir(self.paths['jobs'])]
        self.jobnames = [os.path.splitext(os.path.basename(job))[0] 
                         for job in self.jobs]
        self.outputs = [os.path.join(self.paths['outputs'], out) 
                        for out in os.listdir(self.paths['outputs'])]
        
        return


    def __init__(self, attach = False, resources = None, name = 'registry', verbose = True):

        self.name = name
        self._verbose = verbose
        self.resources = resources
        if attach:
            self._attach_registry()
        else:
            self._create_registry()
        return


    def create_batches(self, x, nbatches = 2, prefix = 'batch'):

        if self._verbose:
            print("Batching data...", flush = True)

        if type(x) is pd.core.frame.DataFrame:
            batch_size = math.ceil(len(x) / nbatches)
            batch_start_ind = range(0, len(x), batch_size)
            self.nbatches = len(batch_start_ind)
            self.batches = []
            for i in range(self.nbatches):
                batch_file = '{}_{}.csv'.format(prefix, i)
                batch_file = os.path.join(self.paths['batches'], batch_file)
                batch_start = batch_start_ind[i]
                batch = x[batch_start:batch_start + batch_size]
                batch.to_csv(batch_file, index = False)
                self.batches.append(batch_file)
        else:
            raise ValueError("No method for type {}".format(type(x)))


    def _create_job(self, script, kwargs = None, prefix = None):

        # Remove null arguments
        if kwargs is not None:
            kwargs = {key:val for key, val in kwargs.items() if val is not None}

        # Create empty temp file
        file = tempfile.NamedTemporaryFile(mode = 'w',
                                           dir = self.paths['jobs'],
                                           prefix = prefix,
                                           suffix = '.sh',
                                           delete = False,
                                           newline = '\n')

        # Job name
        job_resources = self.resources.copy()
        job_resources['job_name'] = os.path.basename(file.name).replace('.sh', '')

        # Job log file
        logfile = os.path.basename(file.name)
        logfile = logfile.replace('.sh', '.out')
        logfile = os.path.join(self.paths['logs'], logfile)
        job_resources['output'] = logfile

        # Open file for editing
        with file as f:

            # Initial shebang
            f.write('#!/bin/bash\n')

            # SBATCH arguments
            for key, val in job_resources.items():
                line = '#SBATCH --{}={}\n'.format(key.replace('_', '-'), val)
                f.write(line)

            # Activate project virtual environment
            # f.write('source activate_venv_hpc.sh\n')

            # Script and arguments
            f.write(script + ' ')
            if kwargs is not None:
                for key, val in kwargs.items():
                    if type(val) is tuple:
                        val = ' '.join(val)
                    line = '--{} {} '.format(key.replace('_', '-'), val)
                    f.write(line)

            f.close()

        return file.name

    def create_jobs(self, script, kwargs = None):
        for i in range(self.nbatches):
            if kwargs is not None:
                kwargs_i = kwargs.copy()
                kwargs_i['input-file'] = self.batches[i]
                kwargs_i['output-file'] = os.path.join(self.paths['outputs'], 'batch_{}.csv'.format(i))
            else:
                kwargs_i = kwargs
            prefix = os.path.splitext(os.path.basename(self.batches[i]))[0]
            prefix = '_'.join([self.name, prefix, '_']) 
            jobfile_i = self._create_job(script = script,
                                         kwargs = kwargs_i,
                                         prefix = prefix)
            self.jobs.append(jobfile_i)
            self.jobnames.append(os.path.basename(jobfile_i))
        self.jobnames = [os.path.splitext(job)[0] for job in self.jobnames]
        

    def _fetch_job_ids(self):
        ids = []
        for job in self.jobnames:
            cmd = 'squeue --me --name={} --format=%i --noheader'.format(job) 
            id = subprocess.run(cmd.split(' '), capture_output = True)
            id = id.stdout.decode('UTF-8').replace('\n', '')
            if not id:
                warn('No job named: {}'.format(job), UserWarning)
            ids.append(id) # This line is being printed in the warning? Don't know why
        return ids


    def _fetch_status(self):
        cmd = ('sacct --jobs={} --name={} --format=JobID%20,JobName%{},STATE%20'
               .format(','.join(self.jobids), ','.join(self.jobnames), 
                       len(self.jobnames[0]))) # This jobnames thing is breaking when double digit batches
        # Fix maybe max([len(job) for job in self.jobnames])
        
        
        output = (subprocess.run(cmd.split(' '), capture_output = True)
                  .stdout.decode('UTF-8').splitlines())
        output = [out.split() for out in output]
        keys = output[0]
        rows = [r for r in output[2:] 
                if r[0] in self.jobids and r[1] in self.jobnames]
        status = {x[1]:[r[x[0]] for r in rows] for x in enumerate(keys)}
        pd.DataFrame(status).to_csv(os.path.join(self.name, 'status.csv'),
                                    index = False)
        return status


    def _check_completion(self):
        status = self._fetch_status()
        codes_active = ['RUNNING', 'PENDING']
        n_active = 0
        for code in codes_active:
            n_active += status['State'].count(code)
        completed = True if n_active == 0 else False
        return completed

    def _reduce(self):
        if self._verbose:
            print("Reducing results...", flush = True)
        result = pd.concat([pd.read_csv(x) for x in self.outputs])
        return result

    def _kill_jobs(self):
        for id in self.jobids:
            _ = subprocess.run(['scancel', id])
        return


    def cleanup(self):
        if self._verbose:
            print("Recursively removing directory: {}".format(self.name), 
                  flush = True)
        shutil.rmtree(self.name)
        return


    def submit_jobs(self, wait = True, cleanup = True):

        try:
            if self._verbose:
                print("Submitting jobs...", flush = True)
            for job in self.jobs:
                _ = subprocess.run(['sbatch', job], capture_output = True)
            self.jobids = self._fetch_job_ids()

            if wait: 
                print("Waiting for results...", flush = True)
                while wait:
                    sleep(10)
                    wait = not self._check_completion()
                self.outputs = [os.path.join(self.paths['outputs'], out) 
                                for out in os.listdir(self.paths['outputs'])]
                return self._reduce()

        except KeyboardInterrupt:
            if self._verbose:
                print("Terminating jobs...", flush = True)
            self._kill_jobs()
            self._fetch_status()

        # else:
        #     self.outputs = [os.path.join(self.paths['outputs'], out) 
        #                     for out in os.listdir(self.paths['outputs'])]
        #     return self._reduce()

        finally:
            if cleanup:
                self.cleanup()
                



@contextmanager
def catchtime() -> float:
    t1 = t2 = perf_counter() 
    yield lambda: t2 - t1
    t2 = perf_counter() 


def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ti = time()
        result = f(*args, **kw)
        tf = time()
        print('%r: %2.4f minutes' % (f.__name__, (tf - ti) / 60))
        return result

    return wrap



# Modules ---------------------------------------------------------------------

def execute_local(script, kwargs = None):
    """
    Execute a program locally.

    Parameters
    ---------
    script: str
        Name of the program to execute.
    kwargs: dict, optional
        Key-value pairs containing command line arguments to pass to the script.

    Returns
    -------
    None
    """

    kwargs = {key:val for key, val in kwargs.items() if val is not None}
    for key, val in kwargs.items():
        if type(val) is tuple or type(val) is list:
            kwargs[key] = [str(x) for x in val]
        else:
            kwargs[key] = [str(val)]
    kwargs = [sum([['--' + str(key)], val], []) for key, val in kwargs.items()]
    kwargs = sum(kwargs, [])
    cmd = [script] + kwargs
    subprocess.run(cmd)
    # output = subprocess.run(cmd, capture_output = True)
    # if output.stderr != b'':
    #     raise Exception("Subprocess failed:\n{}"
    #                     .format(output.stderr.decode("UTF-8")))
    return


def mkdir_from_list(strings, outdir = './', sep = '_'):
    """
    Create a directory from a list of strings.
    
    Parameters
    ---------
    strings: list
        List of strings to use in directory name.
    outdir: str, default './'
        Path to the base directory in which to create the new directory.
    sep: str, default '_'
        Separator used to concatenate the strings in `inlist`.
    
    Returns
    -------
    outdir: str
        Path to the new directory.
    """

    strings = sep.join(strings)
    outdir = os.path.join(outdir, strings)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir


def fetch_params_metadata(metadata, **kwargs):
    """
    Fetch the metadata for a set of parameters.

    Parameters
    ----------
    metadata: str
        Path to the file (.csv) containing parameter set metadata.
    kwargs: dict, optional
        The parameter set to look for. If none is provided, all parameter sets
        found in `metadata` are returned.

    Returns
    -------
    df_match: pandas.core.frame.DataFrame or None
        Parameter sets that matched the parameters provided, or None if no
        matches found.
    """

    # Check if metadata file exists
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


def fetch_params_id(metadata, params):
    """
    Fetch the ID for a parameter set.

    Parameters
    ----------
    metadata: str
        Path to the file (.csv) containing parameter set metadata.
    params: dict
        The parameter set to look for.

    Returns
    -------
    str or None
        The ID of the parameter set, or None if no parameter set matches the
        one provided.

    Raises
    ------
    Exception
        If multiple entries are found for the parameters provided.
    """

    df_metadata = fetch_params_metadata(metadata = metadata, **params)
    if df_metadata is not None:
        nmatch = df_metadata.shape[0]
        if nmatch > 1:
            raise Exception("Multiple entries found for the parameters "
                            "provided. Narrow your search.")
        else:
            return df_metadata['id'].values[0]
    else:
        return None


def set_params_id(metadata, params, params_id = None):
    """
    Identify a parameter set.

    Parameters
    ----------
    metadata: str
        Path to the file (.csv) containing parameter set metadata. If the file
        does not yet exist, it will be created.
    params: dict
        The parameter set to identify.
    params_id: str, default None
        The ID to assign to the parameter set. If None, a random numeric ID is
        generated.

    Returns
    -------
    str
        The ID of the parameter set.
    """

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
    """
    Create a unique directory for a set of parameters.

    Parameters
    ----------
    params: dict
        The parameter set for which to create the directory.
    outdir: str
        Path to the directory in which to create the parameter set directory.
    metadata: str, default None
        Path to the metadata file (.csv) in which to store the parameter
        set information. If None, a file called 'metadata.csv' will be created
        in the directory passed to `outdir`.
    params_id: str, default None
        The ID to assign to the parameter set. If None, a random numeric ID is
        generated.

    Returns
    -------
    outdir: str
        The path to the unique directory for parameter set.
    """

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if metadata is None:
        metadata = os.path.join(outdir, 'metadata.csv')

    try:
        params_id_check = fetch_params_id(metadata = metadata,
                                          params = params)
    except OSError:
        params_id_check = None
    if params_id_check is None:
        params_id = set_params_id(metadata = metadata,
                                  params = params,
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
    
    Parameters
    ---------
    src: list
        Paths to the source files.
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


def random_id(n = 3):
    """
    Create a random ID with a specific number of digits.
    
    Parameters
    ---------
    n: int, default 3
        The number of digits in the random ID.
        
    Returns
    -------
    randid: str
        The random ID.
    """

    num = randint(0, (10 ** n) - 1)
    fstring = '{:0' + str(n) + '}'
    randid = fstring.format(num)
    return randid


def gunzip_file(gzfile, keep = True, outdir = None):
    """
    Unzip a compressed file.
    
    Parameters
    ---------
    gzfile: str
        Path to the file to unzip.
    keep: bool, default True
        Option to keep the compressed file after extraction.
    outdir: str, default None
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


def gunzip_files(infiles, keep = True, outdir = None, nproc = 1):
    """
    Unzip a set of compressed files.
    
    Parameters
    ---------
    infiles: list
        Paths to files to unzip.
    keep: bool, default True
        Option to keep the original compressed files after extraction.
    outdir: str, default None
        Path to the output directory. If None, the files are unzipped in their
        native directories.
    nproc: int, default 1
        Number of processors to use. If `nproc` > 1, the program is executed
        in parallel.
        
    Returns
    -------
    outfiles: list
        Paths to the unzipped files.
    """

    gunzip_file_partial = partial(gunzip_file,
                                  keep = keep,
                                  outdir = outdir)
    if nproc > 1:
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
    
    Parameters
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

    subprocess.run(['nii2mnc', '-quiet', '-clobber', infile])
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
    
    Parameters
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
                   keep = True, outdir = None, nproc = None):
    """
    Convert a set of images between NIFTY and MINC formats.
    
    Parameters
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
    nproc: int
        Number of processors to use. Executed in parallel when > 1.
        
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

    if nproc > 1:
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
            raise Exception("Parameters outdir and suffix "
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
    
    Parameters
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
            raise Exception("Parameters outdir and suffix "
                            "cannot both be None.")

    # Autocrop command
    cmd_autocrop = ['autocrop', '-quiet', '-clobber',
                    '-isostep', str(isostep), infile, outfile]
    log = subprocess.run(cmd_autocrop, stdout = subprocess.PIPE)

    return outfile


def resample_images(infiles, isostep, outdir = None,
                    suffix = None, nproc = None):
    """
    Resample a set of MINC images.
     
    Parameters
    ---------
    infiles: list of str
        List of paths to images to resample.
    isostep: float
        Isotropic dimension of voxels in the resampled image (mm).
    outdir: str, default None
        Path to the output directory. If None, the resampled 
        images will be stored in the native directory.
    suffix: str, default None
        Suffix to append to output filename before the file extension.
    nproc: int, default None
        Number of processors to use. If `nproc` > 1, the program is executed
        in parallel.
        
    Returns
    -------
    outfiles: list of str
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

    if nproc > 1:
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


def transform_coordinates(infile, outfile, transforms, orientation = 'RAS'):
    """
    Transform coordinates from one space to another

    Parameters
    ----------
    infile: str
        Path to the file (.csv) containing the coordinates to transform.
    outfile: str
        Path to the file (.csv) in which to save the transformed
        coordinates.
    transforms: tuple of str
        Paths to the transform files specified in the order of application.
    orientation: str, default 'RAS'
        The world coordinate system of the input coordinates.

    Returns
    -------
    outfile: str
        Path to the file (.csv) containing the transformed coordinates.
    """

    # If RAS orientation, convert to LPS
    if orientation == 'RAS':
        coordinates = pd.read_csv(infile)
        coordinates['x'] = -1 * coordinates['x']
        coordinates['y'] = -1 * coordinates['y']
        infile = infile.replace('.csv', '_LPS.csv')
        coordinates.to_csv(infile, index = False)
    elif orientation == 'LPS':
        pass
    else:
        raise ValueError("`orientation` must be one of {'RAS', 'LPS'}")

    # Concatenate transforms in the order following ANTs convention
    transforms = list(reversed(transforms))
    transforms = [['-t', i] for i in transforms]
    transforms = sum(transforms, [])

    # ANTs transform command
    cmd = ['antsApplyTransformsToPoints',
           '-d', str(3),
           '-i', infile,
           '-o', outfile]
    cmd = cmd + transforms
    subprocess.run(cmd)

    # Return to RAS orientation
    if orientation == 'RAS':
        coordinates = pd.read_csv(outfile)
        coordinates['x'] = -1 * coordinates['x']
        coordinates['y'] = -1 * coordinates['y']
        coordinates.to_csv(outfile, index = False)
        os.remove(infile)

    return outfile


def coordinates_to_minc(coordinates, template, outfile, type = 'labels',
                        verbose = False):
    """

    Parameters
    ----------
    coordinates: str
        Path to the file (.csv) containing coordinates.
    template: str
        Path to the template image (.mnc) defining the space for the
        coordinates.
    outfile: str
        Path to the file (.mnc) in which to export the coordinates
        image.
    type: {'labels', 'mask'}
        Option specifying whether to create a label or mask image
        from the coordinates.
    verbose: bool, default False
        Verbosity option.

    Returns
    -------
    outfile: str
        Path to the file (.mnc) containing the coordinates image.
    """
    kwargs = locals().copy()
    kwargs['verbose'] = 'true' if kwargs['verbose'] else 'false'
    script = 'coordinates_to_minc.R'
    execute_local(script = script, kwargs = kwargs)
    return outfile
