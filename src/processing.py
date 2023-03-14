import os
import subprocess
import tempfile
import numpy                as np
import pandas               as pd
import multiprocessing      as mp
from src.utils              import execute_R
from re                     import sub
from glob                   import glob
from functools              import partial
from tqdm                   import tqdm
from datatable              import fread
from pyminc.volumes.factory import volumeFromFile, volumeLikeFile
from warnings               import warn


def normative_growth_norm(infile, demographics, outfile, key = 'file', df = 5, 
                          combat = False, combat_batch = None, 
                          parallel = False, nproc = None):
    
    """
    Calculate human effect sizes using normative growth modelling.
    
    Arguments
    ---------
    infile: str
        Path to the CSV file containing the voxelwise data.
    demographics: str
        Path to the CSV file containing the human demographics data.
    outfile: str
        Path to the CSV file in which to write the effect size data.
    key: str
        Primary key between voxels and demographics data.
    df: int
        Degrees of freedom to use in natural splines.
    combat: bool
        Option to run ComBat normalization on batch variables prior to 
        normative growth modelling.
    combat_batch: list
        Variables to use in ComBat normalization.
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
    
    Returns
    -------
    outfile: str
        Path to the CSV file containing the effect size data.
    """
  
    #Unpack function args into dictionary
    script_args = locals().copy()
  
    #ComBat normalization options
    if combat:
        script_args['combat'] = 'true'
        if combat_batch is None:
            raise Exception("Argument combat_batch must be specified ",
                            "when combat is True")
        else:
            if type(combat_batch) is str:
                 combat_batch = [combat_batch]
            script_args['combat_batch'] = '-'.join(combat_batch)
            script_args['combat-batch'] = script_args.pop('combat_batch')
    else:
        script_args['combat'] = 'false'
        del script_args['combat_batch']

    #Parallel options
    if parallel:
        script_args['parallel'] = 'true'
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the ",
                             "number of processors to use")
    else:
        script_args['parallel'] = 'false'
    
    #Execute script
    script = 'normative_growth_normalization.R'
    execute_R(script = script, args = script_args)
    return outfile
    
    
def propensity_matching_norm(imgdir, demographics, mask, outdir, 
                             ncontrols = 10, parallel = False, nproc = None):
    
    """
    Calculate human effect size images using propensity-matching.
    
    Arguments
    ---------
    imgdir: str
        Path to the directory containing the MINC images to use to 
        compute the effect sizes.
    demographics: str
        Path to the CSV file containing the human demographics data.
    mask: str
        Path to the mask MINC file for the images.
    outdir: str
        Path to the directory in which to save the effect size MINC
        images.
    ncontrols: int
        Number of propensity-matched controls to use when computing 
        the effect sizes.
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
    
    Returns
    -------
    outfiles: list
        List of paths to the effect size images.
    """
    
    #Unpack function args into dictionary
    script_args = locals().copy()
    
    #Parallel options
    if parallel:
        script_args['parallel'] = 'true'
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the ", 
                             "number of processors to use")
    else:
        script_args['parallel'] = 'false'

    script = 'propensity_matching_normalization.R'
    execute_R(script = script, args = script_args)
    outfiles = glob(outdir+'*.mnc')
    return outfiles
    
    
def vector_to_image(x, outfile, maskfile):
    
    """
    Create a MINC image from a vector of voxel values.
    
    Arguments
    ---------
    x: numpy.ndarray
        Vector of voxel values.
    outfile: str
        Path to the MINC file in which to write the image.
    maskfile: str
        Path to the MINC file containing the mask.
      
    Returns
    -------
    None
    """
    
    #Import image mask
    mask_vol = volumeFromFile(maskfile)
    mask = np.array(mask_vol.data)
    mask_vol.closeVolume()

    #Fill voxels in mask with image values
    img = np.zeros_like(mask.flatten())
    img[(mask == 1).flatten()] = x
    img = img.reshape(mask.shape)

    #Write the image to file
    img_vol = volumeLikeFile(likeFilename = maskfile,
                             outputFilename = outfile,
                             labels = False)
    img_vol.data = img
    img_vol.writeFile()
    img_vol.closeVolume()
    

def matrix_to_images(x, outfiles, maskfile):

    """
    Convert a matrix of voxel values to a set of images.
    
    Arguments
    ---------
    x: numpy.ndarray
        Matrix of voxel values.
    outfiles: list
        List of paths to output MINC files.
    maskfile: str
        Path to the MINC file containing the mask.
    
    Returns
    -------
    outfiles: list
        List of paths to the output images.
    """

    #Function to export image
    exporter = lambda i : vector_to_image(x = x[i,], 
                                          outfile = outfiles[i], 
                                          maskfile = maskfile)
    #Error checking
    if type(x) is not np.ndarray:
        raise ValueError("Argument x must have type numpy.ndarray.")
    else:
        if len(x.shape) != 2:
            raise Exception("Argument x must be a 2-dimensional ",
                            "NumPy array.")
    
    if x.shape[0] != len(outfiles):
        raise Exception("Number of rows in x must be equal to the number ",
                        "of entries in outfiles")
    
    #Iterate over number of files
    irange = range(len(outfiles))
    out = list(map(exporter, tqdm(irange)))
        
    return outfiles

    
def calculate_human_effect_sizes(imgdir, demographics, outdir, method, mask, 
                                 parallel = False, nproc = None, **kwargs):
    
    """
    Calculate human effect size images.
    
    Arguments
    ---------
    imgdir: str
        Path to the directory containing the MINC images to use to 
        compute the effect sizes.
    demographics: str
        Path to the CSV file containing the human demographics data.
    outdir: str
        Path to the directory in which to save the effect size 
        MINC images.
    mask: str
        Path to the mask MINC file for the images.
    method: str
        Method to use to compute effect sizes. 
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
        
    Returns
    -------
    outfiles: list
        List of paths to the effect size images.
    """
    
    #Create output dir if needed
    imgdir = os.path.join(imgdir, '')
    outdir = os.path.join(outdir, '')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    #Compute effect sizes using propensity matching
    if method == "propensity-matching":
        
        kwargs.update({'imgdir':imgdir,
                       'demographics':demographics,
                       'mask':mask,
                       'outdir':outdir,
                       'parallel':parallel,
                       'nproc':nproc})
        outfiles = propensity_matching_norm(**kwargs)
        
    #Compute effect sizes using normative growth modelling
    elif method == "normative-growth":
        
        print("Building voxel matrix...")
  
        #Create a voxel matrix from image files
        imgfiles = glob(imgdir+'*.mnc')
        tmpfile = tempfile.mkstemp(dir = outdir, suffix = '.csv')[1]
        df_voxels = build_voxel_matrix(infiles = imgfiles,
                                       mask = mask,
                                       sort = True,
                                       file_col = True,
                                       parallel = parallel,
                                       nproc = nproc)
        if 'key' in kwargs.keys():
            key = kwargs['key']
        else:
            key = 'file'
        df_voxels[key] = [os.path.basename(file) for file in df_voxels['file']]
        
        print("Writing out voxel matrix...")
        df_voxels.to_csv(tmpfile, index = False)
        
        #Run normative growth normalization
        print("Executing normative growth modelling...")
        outfile = tempfile.mkstemp(dir = outdir, suffix = '.csv')[1]
        kwargs.update({'infile':tmpfile,
                       'demographics':demographics,
                       'outfile':outfile,
                       'parallel':parallel,
                       'nproc':nproc})
        outfile = normative_growth_norm(**kwargs)
        
        #Convert normalized matrix to images
        x = fread(outfile, header = True).to_pandas()
        outfiles = x[key].to_list()
        outfiles = [os.path.join(outdir, outfile) for outfile in outfiles]
        x = x.drop(key, axis=1).to_numpy()
        outfiles = matrix_to_images(x = x, 
                                    outfiles = outfiles,
                                    maskfile = mask)
        
        #Remove temporary files
        os.remove(tmpfile)
        os.remove(outfile)
        
    else:
        raise ValueError("Argument method must be one of ",
                         "['propensity-matching', 'normative-growth']: {}"
                         .format(method))
    
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
        If True and a mask is provided, only the voxels in the mask 
        will be returned.
    
    Returns
    -------
    img: numpy.ndarray
        A NumPy array containing the (masked) image.
    """
    
    #Import image
    img_vol = volumeFromFile(img)
    img_dims = img_vol.getSizes()
    img_seps = img_vol.getSeparations()
    img = np.array(img_vol.data)
    img_vol.closeVolume()

    #Flatten if specified
    if flatten:
        img = img.flatten()

    #Apply mask if specified
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


def import_images(infiles, mask = None, output_format = 'list', flatten = True,
                  parallel = False, nproc = None):

    """
    Import a set of MINC images.
    
    Arguments
    ---------
    infiles: list
        List of paths to images to import.
    mask: str
        Optional path to a mask image. 
    output_format: str
        One of {'list', 'numpy', 'pandas'} indicating what format 
        to return.
    flatten: bool
        Option to flatten images into a 1-dimensional array. 
        If True and a mask is provided, only the voxels in the mask 
        will be returned.
        If False and output_format is not 'list', images will be 
        flattened regardless.
    parallel: bool
        Option to run in parallel.
    nproc: int
        Number of processors to use in parallel.
    
    Returns
    -------
    imgs
        A list, NumPy array, or Pandas DataFrame containing 
        the (masked) images.
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
            raise ValueError("Set the nproc argument to specify the "
                             "number of processors to use in parallel.")
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
    
    
def build_voxel_matrix(infiles, mask = None, file_col = False, sort = False, 
                        save = False, outfile = 'voxel_matrix.csv', 
                        parallel = False, nproc = None):
    
    """
    Create a data frame of voxels from a set of images.
    
    Arguments
    ---------
    infiles: list
        List of paths to images.
    mask: str
        Optional path to a mask image. 
    file_col: bool
        Option to store input files in a column. 
        If true, the paths in infiles are stored in a column 'file'.
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


def cluster_human_data(infiles, rownames = None, nk_max = 10, 
                      metric = 'correlation', K = 10, sigma = 0.5, t = 20, 
                      cluster_file = 'clusters.csv', affinity_file = None):

    """
    Cluster matrices of voxel values using similarity network fusion.
    
    Arguments
    ---------
    infiles: list
        List of paths to the CSV files containing voxelwise data.
    rownames: str
        Column in CSV files containing row names.
    nk_max: int
        Maximum number of clusters to identify. 
        Solutions will be obtained for nk = 2 to nk = nk_max.
    metric: str
        Distance metric used to compute the SNF affinity matrices.
    K: int
        Number of nearest-neighbours used to compute the SNF 
        affinity matrices.
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
    
    #Unpack function args into dictionary
    script_args = locals().copy()
    
    #Input files options
    if type(script_args['infiles']) is not list:
        raise TypeError("Argument infiles must be a list.")
    else:
        if len(script_args['infiles']) != 2:
            raise Exception("Argument infiles must have length 2.")

        script_args['file1'] = script_args['infiles'][0]
        script_args['file2'] = script_args['infiles'][1]
        del script_args['infiles']

    #Update dictionary keys to match R command line args
    script_args['nk-max'] = script_args.pop('nk_max')
    script_args['cluster-file'] = script_args.pop('cluster_file')
    script_args['affinity-file'] = script_args.pop('affinity_file')
            
    #Remove command line args if None
    if rownames is None:
        del script_args['rownames']
        
    if affinity_file is None:
        del script_args['affinity-file']

    script = 'cluster_human_data.R'
    execute_R(script = script, args = script_args)
    return cluster_file



def create_cluster_maps(clusters, imgdir, outdir, mask = None, 
                        method = 'mean', verbose = True):

    """
    Create representative voxel-wise maps for clustered images.
    
    Arguments
    ---------
    clusters: str
        Path to the CSV file containing cluster assignments.
    imgdir: str
        Path to the directory containing images to use.
    mask: str
        Path to the mask file to use.
    outdir: str
        Path to the output directory.
    method: str
        Method used to create the representative cluster maps.
    verbose: bool
        Verbosity option.
        
    Returns
    -------
    outfiles: list
        List of paths to the cluster maps.
    """
        
    #Unpack function args into dictionary
    script_args = locals().copy()
    
    #Check existence of cluster file
    if not os.path.exists(script_args['clusters']):
        raise OSError("Cluster file not found: {}".format(script_args['clusters']))
        
    #Convert bools to str
    script_args['verbose'] = 'true' if script_args['verbose'] else 'false'
    
    #Update dictionary keys to match R command line args
    script_args['cluster-file'] = script_args.pop('clusters')

    #Remove command line args if None
    if mask is None:
        del script_args['mask']
        
    script = 'create_cluster_maps.R'
    execute_R(script = script, args = script_args)
    outfiles = glob(outdir+'*.mnc')
    return outfiles



def threshold_intensity(img, threshold = 0.5, symmetric = True,
                        comparison = 'gt'):
    
    """
    Threshold an image using an intensity threshold.

    Arguments
    --------
    img: numpy.ndarray
        Image to threshold.
    threshold: float
        Intensity value at which to threshold.
    symmetric: bool
        Apply threshold symmetrically to negative and
        positive values.
    comparison: str
        One of 'gt' or 'lt' indicating whether the mask should
        return values greater than or lesser than the threshold.

    Returns
    -------
    out: numpy.ndarray
        Thresholded image.
    """

    out = img.copy()
    if symmetric:
        if comparison == 'gt':
            ind = np.abs(img) > np.abs(threshold)
        elif comparison == 'lt':
            ind = np.abs(img) < np.abs(threshold)
        else:
            raise ValueError("Argument comparison must be one of {'gt', 'lt'}: {}".format(comparison))
    else:
        if comparison == 'gt':
            ind = img > threshold
        elif comparison == 'lt':
            ind = img < threshold
        else:
            raise ValueError("Argument comparison must be one of {'gt', 'lt'}: {}".format(comparison))
    if np.sum(ind) == 0:
        raise Exception("No voxels survive the threshold. Select another threshold.")
    out[~ind] = 0
    return out


def threshold_top_n(img, n = 0.2, symmetric = True, tolerance = 1e-5):

    """
    Threshold an image using the top n values.

    Arguments
    ---------
    img: numpy.ndarray
        Image to threshold.
    n: int or float
        If |n| > 1, the top n values will be used.
        If 0 < |n| < 1, the top fraction of values will be used.
        If n < 0, the top negative values will be used.
    symmetric: bool
        Apply the threshold symmetrically.
    tolerance: float
        Absolute values below this threshold will be excluded.

    Returns
    -------
    out: numpy.ndarray
        Thresholded image.
    """

    # Raise error if symmetric is True and n < 0
    if (symmetric) & (n < 0):
        raise ValueError("Setting n < 0 while symmetric = True ",
                         "will return an empty mask.")

    # Flatten image and mask
    values = img.flatten()

    # If symmetric, use absolute values
    if symmetric:
        values = np.abs(values)

    # Sort values and corresponding indices
    sorted_index = values.argsort()
    sorted_values = values[sorted_index]

    # Tolerance filter
    tolerance_filter = np.abs(sorted_values) > tolerance

    # Compute top n values
    if n > 0:
        positive_filter = sorted_values > 0
        sorted_values = sorted_values[positive_filter & tolerance_filter]
        sorted_index = sorted_index[positive_filter & tolerance_filter]
        if n < 1:
            n = int(np.floor(n * len(sorted_values)))
        top_n_index = sorted_index[-n:]
    elif n < 0:
        negative_filter = sorted_values < 0
        sorted_values = sorted_values[negative_filter & tolerance_filter]
        sorted_index = sorted_index[negative_filter & tolerance_filter]
        n = abs(n)
        if n < 1:
            n = int(np.floor(n * len(sorted_values)))
        top_n_index = sorted_index[:n]
    else:
        raise ValueError("Argument n cannot be 0")

    #Threshold the image
    out = img.flatten().copy()
    index = np.arange(len(out))
    out[~np.isin(index, top_n_index)] = 0
    out = out.reshape(img.shape)

    return out


def threshold_image(img, method = 'top_n', threshold = 0.2, symmetric = True,
                    comparison = 'gt'):

    # Thresholding method
    if method == 'intensity':
        img = threshold_intensity(img = img,
                                  threshold = threshold,
                                  symmetric = symmetric,
                                  comparison = comparison)
    elif method == 'top_n':
        img = threshold_top_n(img = img,
                              n = threshold,
                              symmetric = symmetric)
    else:
        raise ValueError("Argument method must be one of "
                         "{'intensity', 'top_n'}")

    return img


def mask_from_image(img, signed = False):

    mask = np.zeros_like(img)
    if signed:
        mask[img > 0] = 1
        mask[img < 0] = -1
    else:
        mask[np.abs(img) > 0] = 1
    mask = np.int8(mask)
    return mask


# def create_image_mask(infile, outfile, mask, method = 'top_n', threshold = 0.2,
#                       symmetric = True, comparison = None, signed = False):
#     """
#     Create a threshold mask based on an input image.

#     Arguments
#     ---------
#     infile: str
#         Path to the image file (.mnc) to mask.
#     outfile: str
#         Path to the image file (.mnc) in which to save the mask.
#     mask: str
#         Path to the mask image (.mnc) to apply before thresholding.
#     method: str
#         One of {'intensity', 'top_n'} indicating which method
#         to use for thresholding.
#     threshold: float
#         Threshold to determine which voxels are in the mask.
#         If method = 'intensity', this is the intensity value at
#         which to threshold.
#         If method = 'top_n', this is the number of fraction of top
#         voxels to use.
#     symmetric: bool
#         Option to apply threshold symmetrically to negative and
#         positive values.
#     comparison: str
#         String indicating how to apply the threshold if
#         method = 'intensity'. This is ignored if method = 'top_n'.
#     signed: bool
#         If True, the output mask will contain -1 where data < 0
#         and 1 where data > 0

#     Returns
#     -------
#     None
#     """

#     # Import image
#     img = import_image(img = infile, mask = mask, flatten = False)

#     # Thresholding method
#     if method == 'intensity':
#         img = threshold_intensity(img = img,
#                                   threshold = threshold,
#                                   symmetric = symmetric,
#                                   comparison = comparison,
#                                   signed = signed)
#     elif method == 'top_n':
#         img = threshold_top_n(img = img,
#                               n = threshold,
#                               symmetric = symmetric,
#                               signed = signed)
#     else:
#         raise ValueError("Argument method must be one of "
#                          "{'intensity', 'top_n'}")

#     # MINC files can't handle negative integer labels
#     labels = False if signed else True
#     img_vol = volumeLikeFile(likeFilename = infile,
#                              outputFilename = outfile,
#                              labels = labels)

#     img_vol.data = img
#     img_vol.writeFile()
#     img_vol.closeVolume()

#     return outfile

# def create_image_masks(imgdir, outdir, mask, method = 'top_n',
#                        threshold = 0.2, symmetric = True, comparison = None,
#                        signed = False):

#     # Cast threshold to integer if number of voxels used
#     if (method == 'top_n') & (abs(threshold) > 1):
#         threshold = int(threshold)

#     # Ensure proper paths
#     imgdir = os.path.join(imgdir, '')
#     outdir = os.path.join(outdir, '')

#     # Create outdir if needed
#     if not os.path.exists(outdir):
#         os.makedirs(outdir)

#     # Get input images
#     infiles = os.listdir(imgdir)
#     infiles = [file for file in infiles if '.mnc' in file]
#     infiles = [os.path.join(imgdir, file) for file in infiles]
#     outfiles = [os.path.basename(file) for file in infiles]
#     outfiles = [os.path.join(outdir, file) for file in outfiles]

#     #Partial function for iteration
#     masker = partial(create_image_mask,
#                      mask = mask,
#                      method = method,
#                      threshold = threshold,
#                      symmetric = symmetric,
#                      comparison = comparison,
#                      signed = signed)

#     #Iterate over images
#     outfiles = list(map(masker, tqdm(infiles), outfiles))

#     return outfiles
