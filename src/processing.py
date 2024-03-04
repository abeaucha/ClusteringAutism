import os
import random
import tempfile
import multiprocessing as mp
import numpy as np
import pandas as pd
from functools import partial
from glob import glob
from pyminc.volumes.factory import volumeFromFile, volumeLikeFile
from shutil import rmtree
from tqdm import tqdm
from utils import execute_R
from warnings import warn


def normative_growth_norm(imgdir, demographics, mask, outdir, key = 'file',
                          group = 'patients', df = 3, batch = None,
                          nbatches = 1, nproc = None):
    """
    Calculate human effect sizes using normative growth modelling.

    Arguments
    ---------
    imgdir: str
        Path to the directory containing the images (.mnc) to use to
        compute the effect sizes.
    demographics: str
        Path to the file (.csv) containing the demographics data.
    mask: str
        Path to the mask file (.mnc).
    outdir: str
        Path to the directory in which to save the effect size images.
    key: str, default 'file'
        Primary key between demographics data and constructed voxel matrix.
    group: {'patients', 'controls', 'all'}
        Group of participants for which to compute effect sizes.
    df: int, default 3
        Degrees of freedom to use in normative model natural splines.
    batch: list of str
        Variables to use in normalization prior to modelling.
    nbatches: int, default 1
        Number of voxel batches to use for computation of voxel-wise
        normative models.
    nproc: int, default None
        Number of processors to use in parallel. Executed sequentially when
        None.

    Returns
    -------
    outfiles: list of str
        List of paths to the effect size images.
    """

    # Unpack function args into dictionary
    script_args = dict(imgdir = imgdir,
                       demographics = demographics,
                       mask = mask,
                       outdir = outdir,
                       key = key,
                       group = group,
                       df = df,
                       batch = batch,
                       nproc = nproc)

    # ComBat normalization options
    if batch is not None:
        if type(batch) is str:
            batch = [batch]
        script_args['batch'] = '-'.join(batch)
    else:
        del script_args['batch']

    # Parallel check
    if nproc is None:
        raise ValueError("Set the nproc argument to specify the number of "
                         "processors to use")

    # R script
    script = 'normative_growth_normalization.R'

    # Option to run in batches
    if nbatches > 1:

        print("Executing script in batches...")

        # Create voxel batches
        mask_array = import_image(img = mask, mask = mask)
        mask_ind = np.arange(len(mask_array))
        batches = np.array_split(mask_ind, nbatches)

        # Batch directories
        batch_dirs = ['batch_{}'.format(batch) for batch in range(nbatches)]
        batch_dirs = [os.path.join(outdir, dir, '') for dir in batch_dirs]

        # Execute script in batches
        for i, batch in enumerate(batches):

            print("Processing batch {}...".format(i))

            # Create batch mask
            batch_mask = np.zeros(len(mask_array))
            batch_mask[batch] = 1

            # Export batch mask
            batch_dir = batch_dirs[i]
            if not os.path.exists(batch_dir):
                os.makedirs(batch_dir)
            batch_maskfile = os.path.join(batch_dir, 'batch_mask.mnc')
            vector_to_image(x = batch_mask, outfile = batch_maskfile,
                            maskfile = mask)

            # Update script args with batch paths
            script_args.update(dict(mask = batch_maskfile, outdir = batch_dir))

            # Execute script
            execute_R(script = script, **script_args)

        # Collate batch images
        print("Collating batched images...")
        outfiles = os.listdir(batch_dirs[0])
        outfiles = [file for file in outfiles if file != 'batch_mask.mnc']
        for outfile in outfiles:

            img = np.zeros_like(mask_array)
            for b, batch in enumerate(batches):
                batch_img = os.path.join(batch_dirs[b], outfile)
                batch_mask = os.path.join(batch_dirs[b], 'batch_mask.mnc')
                img[batch] = import_image(img = batch_img, mask = batch_mask)

            outfile = os.path.join(outdir, outfile)
            vector_to_image(x = img, outfile = outfile, maskfile = mask)

        print("Cleaning up...")
        for batch_dir in batch_dirs:
            rmtree(batch_dir)

    else:

        execute_R(script = script, **script_args)

    outfiles = os.listdir(outdir)
    outfiles = [file for file in outfiles if '.mnc' in file]
    outfiles = [os.path.join(outdir, file) for file in outfiles]
    return outfiles


def propensity_matching_norm(imgdir, demographics, mask, outdir,
                             ncontrols = 10, parallel = False, nproc = None):
    """
    Calculate human effect size images using propensity-matching.
    
    Arguments
    ---------
    imgdir: str
        Path to the directory containing the images (.mnc) to use to
        compute the effect sizes.
    demographics: str
        Path to the file (.csv) containing the demographics data.
    mask: str
        Path to the mask file (.mnc).
    outdir: str
        Path to the directory in which to save the effect size images.
    ncontrols: int, default 10
        Number of propensity-matched controls to use when computing 
        the effect sizes.
    parallel: bool, default False
        Option to run in parallel.
    nproc: int, default None
        Number of processors to use in parallel.
    
    Returns
    -------
    outfiles: list of str
        List of paths to the effect size images.
    """

    # Unpack function args into dictionary
    script_args = locals().copy()

    # Parallel options
    if parallel:
        script_args['parallel'] = 'true'
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the ",
                             "number of processors to use")
    else:
        script_args['parallel'] = 'false'

    script = 'propensity_matching_normalization.R'
    execute_R(script = script, **script_args)
    outfiles = glob(outdir + '*.mnc')
    return outfiles


def vector_to_image(x, outfile, maskfile, version = "v1"):
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

    # Import image mask
    mask_vol = volumeFromFile(maskfile)
    mask = np.array(mask_vol.data)
    mask_vol.closeVolume()

    # Fill voxels in mask with image values
    img = np.zeros_like(mask.flatten())
    if version == "v2":
        img[(mask > 0.5).flatten()] = x
    else:
        img[(mask == 1).flatten()] = x
    img = img.reshape(mask.shape)

    # Write the image to file
    img_vol = volumeLikeFile(likeFilename = maskfile,
                             outputFilename = outfile,
                             labels = False)
    img_vol.data = img
    img_vol.writeFile()
    img_vol.closeVolume()


def matrix_to_images(x, outfiles, maskfile, version = "v1"):
    """
    Convert a matrix of voxel values to a set of images.
    
    Arguments
    ---------
    x: numpy.ndarray
        Matrix of voxel values.
    outfiles: list
        Paths to the output files (.mnc).
    maskfile: str
        Path to the MINC file containing the mask.
    
    Returns
    -------
    outfiles: list
        List of paths to the output images.
    """

    # Function to export image
    exporter = lambda i:vector_to_image_v2(x = x[i,],
                                           outfile = outfiles[i],
                                           maskfile = maskfile,
                                           version = version)

    # Error checking
    if type(x) is not np.ndarray:
        raise ValueError("Argument x must have type numpy.ndarray.")
    else:
        if len(x.shape) != 2:
            raise Exception("Argument x must be a 2-dimensional ",
                            "NumPy array.")

    if x.shape[0] != len(outfiles):
        raise Exception("Number of rows in x must be equal to the number ",
                        "of entries in outfiles")

    # Iterate over number of files
    irange = range(len(outfiles))
    out = list(map(exporter, tqdm(irange)))

    return outfiles


def calculate_human_effect_sizes(imgdir, demographics, mask, outdir,
                                 method = 'normative-growth', parallel = False,
                                 nproc = None, **kwargs):
    """
    Calculate human effect size images.
    
    Arguments
    ---------
    imgdir: str
        Path to the directory containing the images (.mnc) to use to
        compute the effect sizes.
    demographics: str
        Path to the file (.csv) containing the demographics data.
    mask: str
        Path to the mask file (.mnc).
    outdir: str
        Path to the directory in which to save the effect size images.
    method: {'normative-growth', 'propensity-matching'}
        Method to use to compute effect sizes.
    parallel: bool, default False
        Option to run in parallel.
    nproc: int, default None
        Number of processors to use in parallel.
    **kwargs: dict, optional
        Extra arguments passed to `normative_growth_norm` and
        `propensity_matching_norm` sub-routines. Refer to documentation for a
        list of possible arguments.
        
    Returns
    -------
    outfiles: list of str
        List of paths to the effect size images.
    """

    # Create output dir if needed
    imgdir = os.path.join(imgdir, '')
    outdir = os.path.join(outdir, '')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Compute effect sizes using propensity matching
    if method == 'propensity-matching':

        kwargs.update({'imgdir':imgdir,
                       'demographics':demographics,
                       'mask':mask,
                       'outdir':outdir,
                       'parallel':parallel,
                       'nproc':nproc})
        outfiles = propensity_matching_norm(**kwargs)

    # Compute effect sizes using normative growth modelling
    elif method == 'normative-growth':

        kwargs.update(dict(imgdir = imgdir,
                           demographics = demographics,
                           mask = mask,
                           outdir = outdir,
                           nproc = nproc))
        outfiles = normative_growth_norm(**kwargs)

    else:
        raise ValueError("Argument method must be one of ",
                         "('propensity-matching', 'normative-growth'): {}"
                         .format(method))

    return outfiles


def import_image(img, mask = None, flatten = True, version = "v1"):
    """
    Import a MINC image.
    
    Arguments
    ---------
    img: str
        Path to the image (.mnc) to import.
    mask: str, default None
        Path to a mask image (.mnc).
    flatten: bool, default True
        Option to flatten image into a 1-dimensional array. If True and a mask
        is provided, only the voxels in the mask will be returned.
    
    Returns
    -------
    img: numpy.ndarray
        Image voxel values.
    """

    # Import image
    img_vol = volumeFromFile(img)
    img_dims = img_vol.getSizes()
    img_seps = img_vol.getSeparations()
    img = np.array(img_vol.data)
    img_vol.closeVolume()

    # Flatten if specified
    if flatten:
        img = img.flatten()

    # Apply mask if specified
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
            if version == "v2":
                img = img[mask > 0.5]
            else:
                img = img[mask == 1]
        else:
            if version == "v2":
                img[mask < 0.5] = 0
            else :
                img[mask == 0] = 0

    return img


def import_images(imgfiles, mask = None, output_format = 'list', flatten = True,
                  version = "v1", parallel = False, nproc = None):
    """
    Import a set of MINC images.
    
    Arguments
    ---------
    imgfiles: list of str
        List of paths to images.
    mask: str, default None
        Path to a mask image (.mnc).
    output_format: {'list', 'numpy', 'pandas'}
        String indicating what format to return.
    flatten: bool, default True
        Option to flatten images into a 1-dimensional array. If True and a mask
        is provided, only the voxels in the mask will be returned. If False and
        `output_format` is not 'list', images will be  flattened regardless.
    parallel: bool, default False
        Option to run in parallel.
    nproc: int, default None
        Number of processors to use in parallel.
    
    Returns
    -------
    imgs: list or numpy.ndarray or pandas.DataFrame
        Object containing image voxel values.
    """

    format_opts = ['list', 'numpy', 'pandas']
    format_test = sum([output_format == opt for opt in format_opts])
    format_err = ("Argument output_format must be one of {}: {}"
                  .format(format_opts, output_format))
    if format_test != 1:
        raise ValueError(format_err)

    if not flatten:
        if output_format != 'list':
            msg_warn = (
                "flatten = False is only valid when output_format = 'list'. "
                "Proceeding with flattened images.")
            warn(msg_warn)
            flatten = True

    importer = partial(import_image,
                       mask = mask,
                       flatten = flatten,
                       version = version)

    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the "
                             "number of processors to use in parallel.")
        pool = mp.Pool(nproc)
        imgs = []
        for img in tqdm(pool.imap(importer, imgfiles), total = len(imgfiles)):
            imgs.append(img)
        pool.close()
        pool.join()
    else:
        imgs = list(map(importer, tqdm(imgfiles)))

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


def build_voxel_matrix(imgfiles, mask = None, file_col = False, sort = False,
                       save = False, outfile = 'voxel_matrix.csv',
                       version = "v1", parallel = False, nproc = None):
    """
    Create a data frame of voxels from a set of images.
    
    Arguments
    ---------
    imgfiles: list of str
        List of paths to images.
    mask: str, default None
        Path to a mask image (.mnc).
    file_col: bool, default False
        Option to store input files in a column. If true, the paths in
        `imgfiles` are stored in a column 'file'.
    sort: bool, default False
        Option to sort rows based on file names.
    save: bool, default False
        Option to save to CSV.
    outfile: str, default 'voxel_matrix.csv'
        Path to the output CSV file. Ignored if `save` = False.
    parallel: bool, default False
        Option to run in parallel.
    nproc: int, default None
        Number of processors to use in parallel.
    
    Returns
    -------
    df_imgs: pandas.DataFrame
        A data frame containing the image voxel values.
    """

    df_imgs = import_images(imgfiles = imgfiles,
                            mask = mask,
                            output_format = 'pandas',
                            flatten = True,
                            version = version,
                            parallel = parallel,
                            nproc = nproc)
    df_imgs['file'] = imgfiles

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
    infiles: list of str
        List of paths to the files (.csv) containing voxel-wise data.
    rownames: str, default None
        Column in the input files containing row names.
    nk_max: int, default 10
        Maximum number of clusters to identify. Solutions will be obtained for
        nk = 2 to nk = `nk_max`.
    metric: str, default 'correlation'
        Similarity metric used to compute the SNF affinity matrices.
    K: int, default 10
        Number of nearest-neighbours used to compute the SNF 
        affinity matrices.
    sigma: float, default 0.5
        Variance for the local model in the SNF affinity matrices.
    t: int, default 20
        Number of iterations for the diffusion process in SNF.
    cluster_file: str, default 'clusters.csv'
        Path to the file (.csv) in which to write the cluster assignments.
    affinity_file: str, default None
        Path to the file (.csv) in which to write the affinity matrix. If None,
        the affinity matrix is not saved.
        
    Returns
    -------
    cluster_file: str
    Path to the file (.csv) containing the cluster assignments
    """

    # Unpack function args into dictionary
    script_args = locals().copy()

    # Input files options
    if type(script_args['infiles']) is not list:
        raise TypeError("Argument infiles must be a list.")
    else:
        if len(script_args['infiles']) != 2:
            raise Exception("Argument infiles must have length 2.")

        script_args['file1'] = script_args['infiles'][0]
        script_args['file2'] = script_args['infiles'][1]
        del script_args['infiles']

    # Update dictionary keys to match R command line args
    script_args['nk-max'] = script_args.pop('nk_max')
    script_args['cluster-file'] = script_args.pop('cluster_file')
    script_args['affinity-file'] = script_args.pop('affinity_file')

    # Remove command line args if None
    if rownames is None:
        del script_args['rownames']

    if affinity_file is None:
        del script_args['affinity-file']

    script = 'cluster_human_data.R'
    execute_R(script = script, **script_args)
    return cluster_file


def create_cluster_maps(clusters, imgdir, outdir, mask = None,
                        method = 'mean', nproc = 2, verbose = True):
    """
    Create cluster centroid images.
    
    Arguments
    ---------
    clusters: str
        Path to the file (.csv) containing cluster assignments.
    imgdir: str
        Path to the directory containing images.
    outdir: str
        Path to the output directory.
    mask: str, default None
        Path to a mask image (.mnc).
    method: {'mean', 'median'}
        Method used to compute the centroid images.
    nproc: int, default 2
        Number of processors to use in parallel.
    verbose: bool
        Verbosity option.
        
    Returns
    -------
    outfiles: list of str
        Paths to the cluster centroid images.
    """

    # Unpack function args into dictionary
    script_args = locals().copy()

    # Check existence of cluster file
    if not os.path.exists(script_args['clusters']):
        raise OSError(
            "Cluster file not found: {}".format(script_args['clusters']))

    # Convert bools to str
    script_args['verbose'] = 'true' if script_args['verbose'] else 'false'

    # Update dictionary keys to match R command line args
    script_args['cluster-file'] = script_args.pop('clusters')

    # Remove command line args if None
    if mask is None:
        del script_args['mask']

    script = 'create_cluster_maps.R'
    execute_R(script = script, **script_args)
    outfiles = glob(outdir + '*.mnc')
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
            raise ValueError(
                "Argument comparison must be one of ['gt', 'lt']: {}".format(
                    comparison))
    else:
        if comparison == 'gt':
            ind = img > threshold
        elif comparison == 'lt':
            ind = img < threshold
        else:
            raise ValueError(
                "Argument comparison must be one of ['gt', 'lt']: {}".format(
                    comparison))
    if np.sum(ind) == 0:
        raise Exception(
            "No voxels survive the threshold. Select another threshold.")
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
    if symmetric & (n < 0):
        raise ValueError("Setting n < 0 while symmetric = True "
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

    # Threshold the image
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
                         "('intensity', 'top_n')")

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


def permute_cluster_labels(cluster_file, outdir, npermutations = 100, start = 1,
                           min_per_k = None):
    outdir = os.path.join(outdir, '')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    df_clusters = pd.read_csv(cluster_file)
    cols = df_clusters.drop('ID', axis = 1).columns

    min_per_k = 0 if min_per_k is None else min_per_k

    outfiles = []
    for p in range(start, start + npermutations):

        df_permute = df_clusters.copy()
        for col in cols:
            k = df_clusters[col].to_numpy()
            k_vals, k_freq = np.unique(k, return_counts = True)
            k_lt_min = k_vals[k_freq < min_per_k]
            lt_min = np.isin(k, k_lt_min)

            random.seed(p)
            np.random.seed(p)
            k[~lt_min] = np.random.choice(k[~lt_min],
                                          size = len(k[~lt_min]),
                                          replace = False)

            k[lt_min] = -1
            df_permute[col] = k

        outfile = 'clusters_permutation_{}.csv'.format(p)
        outfile = os.path.join(outdir, outfile)
        df_permute.to_csv(outfile, index = False)
        outfiles.append(outfile)

    return outfiles

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
