import contextlib
import os
import sys
import utils
import numpy as np
import pandas as pd
import processing
#from datatable import fread
from dask.distributed import Client, progress
from dask_jobqueue import SLURMCluster
from functools import partial
from io import StringIO
from itertools import product
from pyminc.volumes.factory import volumeFromFile


@contextlib.contextmanager
def silence():
    save_stdout = sys.stdout
    sys.stdout = StringIO()
    yield
    sys.stdout = save_stdout


def fetch_microarray_coordinates(metadata, outfile, labels = True,
                                 verbose = False):
    """
    Fetch the AHBA microarray sample coordinates from the web.

    Parameters
    ----------
    metadata: str
        Path to file (.csv) containing the AHBA microarray sample
        metadata.
    outfile: str
        Path to the file (.csv) in which to save the sample coordinates.
    labels: bool, default True
        Option to attach atlas labels for the individual microarray
        samples.
    verbose: bool, default False
        Verbosity option.

    Returns
    -------
    outfile: str
        Path to the file (.csv) containing the sample coordinates.
    """
    kwargs = locals().copy()
    kwargs['labels'] = 'true' if kwargs['labels'] else 'false'
    kwargs['verbose'] = 'true' if kwargs['verbose'] else 'false'
    script = 'fetch_microarray_coordinates.R'
    utils.execute_local(script = script, kwargs = kwargs)
    return outfile


def prepare_microarray_coordinates(metadata, transforms,
                                   outdir = 'data/human/expression/',
                                   annotations = None):
    """
    Prepare the AHBA microarray coordinates for the imaging study.

    Parameters
    ----------
    metadata: str
        Path to file (.csv) containing the AHBA microarray sample
        metadata.
    transforms: tuple of str
        Transforms from the MNI ICBM 152 NLIN 09c space to the
        registration consensus average space, to be passed to
        antsApplyTransformsToPoints.
    outdir: str
        Path to the directory in which to export the microarray
        sample coordinates.
    annotations: str
        Path to the file (.csv) containing AHBA microarray sample
        annotations indicating whether to keep or discard the sample.

    Returns
    -------
    coords: str
        Path to the file (.csv) containing the AHBA microarray sample
        coordinates in the imaging study space.
    """

    # MNI coordinates output file
    coords_mni = 'AHBA_microarray_coordinates_mni.csv'
    coords_mni = os.path.join(outdir, coords_mni)

    # Fetch microarray coordinates
    coords_mni = fetch_microarray_coordinates(metadata = metadata,
                                              outfile = coords_mni,
                                              labels = True)

    # Microarray coordinates label definitions
    defs = coords_mni.replace('.csv', '_defs.csv')

    if annotations is not None:
        # Microarray sample annotations
        annotations = pd.read_csv(annotations)

        # Filter MNI coordinates
        (pd.read_csv(coords_mni)
         .loc[annotations['keep']]
         .to_csv(coords_mni, index = False))

        # Filter MNI coordinates label definitions
        (pd.read_csv(defs)
         .loc[annotations['keep']]
         .to_csv(defs, index = False))

    # Transform MNI coordinates to study space
    coords_study = coords_mni.replace('mni', 'study')
    coords = utils.transform_coordinates(infile = coords_mni,
                                         outfile = coords_study,
                                         transforms = transforms,
                                         orientation = 'RAS')

    return coords


def get_latent_spaces(expr, ids = None):
    """
    Fetch gene expression latent space files

    Parameters
    ----------
    expr: list of str
        Path to the directory containing the latent space files (.csv).
    ids: list of int
        List containing latent space IDs to return.

    Returns
    -------
    latent_spaces: list of tuple of str
        Paths to the files (.csv) containing the latent space data.
    """
    latent_spaces = []
    for e in expr:
        files = os.listdir(e)
        ids_all = [int(i.split('_')[-1]) for i in
                   [os.path.splitext(file)[0] for file in files]]
        ids_files_all = sorted(zip(ids_all, files))
        if ids is not None:
            ls = [os.path.join(e, file) for i, file in ids_files_all if
                  i in ids]
        else:
            ls = [os.path.join(e, file) for i, file in ids_files_all]
        latent_spaces.append(ls)
    latent_spaces = list(zip(latent_spaces[0], latent_spaces[1]))

    return latent_spaces


def mouse_signature(img: str, expr: str, mask: str, signed: bool = True,
                    threshold: str = 'top_n', threshold_value: float = 0.2,
                    threshold_symmetric: bool = True):
    """
    Compute the transcriptomic signature of a mouse image.

    Description
    -----------


    Arguments
    ---------
    img: str
        Path to the image file (.mnc).
    expr: str
        Path to the expression matrix file (.csv).
    mask: str
        Path to the mask image file used to construct the expression matrix
        (.mnc).
    signed: bool, default True
    threshold: {'top_n', 'intensity'}
        Method used to threshold the image prior to compute the transcriptomic
        signature. If None, all non-zero voxels that fall within the mask will
        be used.
    threshold_value: float, default 0.2
        Threshold for the thresholding method. Ignored if `threshold` is None.
    threshold_symmetric: bool, default True
        Option to apply the thresholding method symmetrically to positive and
        negative voxels.

    Returns
    -------
    signature:
    """

    # Import image
    img = processing.import_image(img = img, mask = mask, flatten = False)

    # Import brain mask
    mask = processing.import_image(img = mask, flatten = True)

    # Import transcriptomic data
    with silence():
        # expr = fread(expr, header = True).to_pandas()
        expr = pd.read_csv(expr)
    if np.sum(mask) != expr.shape[0]:
        raise Exception("The number of voxels in `mask` does not match the "
                        "number of rows in `expr`.")

    # Threshold the image if desired
    if threshold is not None:
        img = processing.threshold_image(img = img,
                                         method = threshold,
                                         threshold = threshold_value,
                                         symmetric = threshold_symmetric)
    else:
        img = processing.threshold_image(img = img,
                                         method = 'intensity',
                                         threshold = 0,
                                         symmetric = True)

    # Create a mask from the (thresholded) image
    img_mask = processing.mask_from_image(img = img, signed = signed)
    img_mask = img_mask.flatten()
    img_mask = img_mask[mask == 1]

    # If signed, return signatures for positive and negative image values
    if signed:
        signature_pos = expr.loc[img_mask == 1].mean().to_numpy()
        signature_neg = expr.loc[img_mask == -1].mean().to_numpy()
        signature = (signature_pos, signature_neg)
    else:
        signature = expr.loc[img_mask == 1].mean().to_numpy()

    return signature


def human_signature(img, expr, mask, coords, signed = True, threshold = 'top_n',
                    threshold_value = 0.2, threshold_symmetric = True):
    """
    Compute the transcriptomic signature for a human image
    """

    # Import image
    img = processing.import_image(img = img, mask = mask, flatten = False)

    # Import transcriptomic data
    with silence():
        # expr = fread(expr, header = True).to_pandas()
        expr = pd.read_csv(expr)

    # Import microarray coordinates
    coords = pd.read_csv(coords)

    # Threshold the image if desired
    if threshold is not None:
        img = processing.threshold_image(img = img,
                                         method = threshold,
                                         threshold = threshold_value,
                                         symmetric = threshold_symmetric)
    else:
        img = processing.threshold_image(img = img,
                                         method = 'intensity',
                                         threshold = 0,
                                         symmetric = True)

    # Create a mask from the (thresholded) image
    img_mask = processing.mask_from_image(img = img, signed = signed)

    # Create a mask for voxels with microarray data
    vol = volumeFromFile(mask)
    microarray_mask = np.zeros(coords.shape[0], dtype = int)
    for i, row in coords.iterrows():
        coords_world = np.array([row['x'], row['y'], row['z']])
        coords_voxel = vol.convertWorldToVoxel(coords_world)
        coords_voxel = np.round(coords_voxel)
        x = int(coords_voxel[0]) - 1
        y = int(coords_voxel[1]) - 1
        z = int(coords_voxel[2]) - 1
        microarray_mask[i] = img_mask[x, y, z]
    vol.closeVolume()

    if signed:
        signature_pos = expr.loc[microarray_mask == 1].mean().to_numpy()
        signature_neg = expr.loc[microarray_mask == -1].mean().to_numpy()
        signature = (signature_pos, signature_neg)
    else:
        signature = expr.loc[microarray_mask == 1].mean().to_numpy()

    return signature


def similarity(x, y, metric = 'correlation'):
    """
    Compute the similarity between two vectors
    """

    if (type(x) != np.ndarray) | (type(y) != np.ndarray):
        raise TypeError("x and y must have type numpy.ndarray")

    if (len(x.shape) != 1) | (len(y.shape) != 1):
        raise Exception("x and y must be 1-dimensional arrays")

    if len(x) != len(y):
        raise Exception("x and y must have the same length")

    if metric == 'correlation':
        d = np.corrcoef(x = x, y = y)
        d = d[0, 1]
    elif metric == 'cosine':
        d = np.dot(x, y) / np.sqrt((np.dot(x, x) * np.dot(y, y)))
    elif metric == 'euclidean':
        d = np.linalg.norm(x - y, ord = 2)
    elif metric == 'manhattan':
        d = np.sqrt(np.linalg.norm(x - y, ord = 1))
    elif metric == 'KL-divergence':
        d = np.sum(np.where(x != 0, x * np.log(x / y), 0))
    else:
        raise ValueError("Invalid metric: {}".format(metric))

    return d


def compute_transcriptomic_similarity(imgs, species, expr, masks,
                                      microarray_coords = None,
                                      metric = 'correlation', signed = True,
                                      threshold = 'top_n',
                                      threshold_value = 0.2,
                                      threshold_symmetric = True,
                                      return_signed = False):
    """
    Compute the similarity between a human image and a mouse image.

    Parameters
    ----------
    imgs: tuple of str
        A tuple of length 2 containing the paths to the human and mouse images
        (.mnc) to compare.
    expr: tuple of str
        A tuple of length 2 containing the paths to the human and mouse
        expression matrices (.csv).
    masks: tuple of str
        A tuple of length 2 containing the paths to the human and mouse mask
        images (.mnc).
    microarray_coords: str
        The path to the human microarray sample coordinates file (.csv).
    metric: str, default 'correlation'
        The metric used to compute the similarity of mouse and human images.
    signed: bool, default True
        Option to evaluate similarity based on positive and negative image
        values before combining into a single similarity value. If False,
        similarity is calculated based on absolute voxel values.
    threshold: {'top_n', 'intensity', None}
        Method used to threshold mouse and human images before evaluating
        similarity.
    threshold_value: float, default 0.2
        Threshold value to use with threshold method.
    threshold_symmetric: bool, default True
        Option to apply threshold symmetrically to positive and negative voxel
        values.

    Returns
    -------
    sim: float
        The similarity between the images.
    """

    # Unpack the input tuples
    img_human, img_mouse = imgs
    expr_human, expr_mouse = expr
    mask_human, mask_mouse = masks

    # Evaluate the expression signatures
    signatures = []
    for inputs in zip(species, imgs, expr, masks):
        kwargs = dict(
            img = inputs[1],
            expr = inputs[2],
            mask = inputs[3],
            signed = signed,
            threshold = threshold,
            threshold_value = threshold_value,
            threshold_symmetric = threshold_symmetric
        )
        if inputs[0] == 'human':
            signature = human_signature
            kwargs['coords'] = microarray_coords
        elif inputs[0] == 'mouse':
            signature = mouse_signature
        else:
            raise ValueError
        signatures.append(signature(**kwargs))

    # Evaluate the similarity of the expression signatures
    if signed:
        sim = []
        for sig in zip(signatures[0], signatures[1]):
            sim.append(similarity(x = sig[0], y = sig[1],
                                  metric = metric))

        is_nan = np.isnan(sim)
        if all(is_nan):
            raise Exception("No surviving voxels.")
        if not return_signed:
            sim = np.array(sim)
            sim = sim[~is_nan]
            sim = np.mean(sim)

    else:
        sim = similarity(x = signatures[0],
                         y = signatures[1],
                         metric = metric)

    return sim


def _compute_transcriptomic_similarity(inputs, **kwargs):
    return compute_transcriptomic_similarity(imgs = inputs[0],
                                             expr = inputs[1],
                                             **kwargs)


def transcriptomic_similarity(imgs, species, expr, masks,
                              microarray_coords = None,
                              gene_space = 'avg-mlp-latent-space',
                              n_latent_spaces = 100, latent_space_id = 1,
                              metric = 'correlation', signed = True,
                              threshold = 'top_n', threshold_value = 0.2,
                              threshold_symmetric = True, return_signed = False,
                              execution = 'local', nproc = 1, verbose = True,
                              slurm_kwargs = None):
    """
    Compute the transcriptomic similarity for a set of image pairs.

    Parameters
    ----------
    imgs: tuple of str or list of tuple of str
        A tuple of length 2 containing the paths to the images (.mnc)
        to compare. Multiple pairs of images can be passed as a list
        of tuples.
    species: tuple of str
        A tuple of length 2 indicating which species are being compared.
    expr: tuple of str
        A tuple of length 2 containing the paths to the expression
        directories.
    masks: tuple of str
        A tuple of length 2 containing the paths to the images (.mnc).
    microarray_coords: str, default None
        The path to the human microarray sample coordinates file (.csv).
        Ignored if only mouse images are compared.
    gene_space: {'avg-mlp-latent-space', 'mlp-latent-space', 'vae-latent-space', 'homologous-genes'}
        The gene expression common space to use for comparison.
    n_latent_spaces: int, default 100
        The number of latent spaces to aggregate when `gene_space` =
        'avg-mlp-latent-space'. Ignored otherwise.
    latent_space_id: int, default 1
        The ID of the latent space to use when 'gene_space' = 'mlp-latent-space'.
        Ignored otherwise.
    metric: str, default 'correlation'
        The metric used to compute the similarity of the images.
    signed: bool, default True
        Option to evaluate similarity based on positive and negative image
        values before combining into a single similarity value. If False,
        similarity is calculated based on absolute voxel values.
    threshold: {'top_n', 'intensity', None}
        Method used to threshold the images before evaluating similarity.
    threshold_value: float, default 0.2
        Threshold value to use with threshold method.
    threshold_symmetric: bool, default True
        Option to apply threshold symmetrically to positive and negative voxel
        values.
    nproc: int, default 1
        Number of processors to use in parallel.

    Returns
    -------
    results: pandas.DataFrame
        A data frame containing the similarity values for all input images.
    """

    if imgs is None:
        raise ValueError("Argument `imgs` is required.")
    if species is None:
        raise ValueError("Argument `species` is required.")
    if expr is None:
        raise ValueError("Argument `expr` is required.")
    if masks is None:
        raise ValueError("Argument `masks` is required.")
    if 'human' in species:
        if microarray_coords is None:
            raise ValueError("Argument `microarray_coords` is required when "
                             "one of the species is 'human'.")

    # If imgs is tuple, convert to list of tuple
    if type(imgs) is tuple:
        imgs = [imgs]

    # Check imgs lengths
    imgs_test = [len(i) for i in imgs]
    if set(imgs_test) != {2}:
        raise Exception("`imgs` must either be a tuple of length 2 or a list "
                        "of tuples of length 2.")

    # Check imgs types
    imgs_test = [type(i) for i in imgs]
    if set(imgs_test) != {tuple}:
        raise Exception("`imgs` must either be a tuple of length 2 or a list "
                        "of tuples of length 2.")

    if gene_space == 'homologous-genes':

        # In case of homologous gene space, return a list of length 1
        # containing a tuple with the expression files
        expr = list(expr)
        expr_files = dict(
            human = 'HumanExpressionMatrix_samples_pipeline_abagen_homologs_scaled.csv',
            mouse = 'MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv'
        )
        for i, path in enumerate(expr):
            expr[i] = os.path.join(path, 'input_space', expr_files[species[i]])
        expr = [tuple(expr)]

    elif gene_space == 'mlp-latent-space':

        expr = [os.path.join(path, 'mlp_latent_space') for path in expr]
        expr = get_latent_spaces(expr = expr, ids = [latent_space_id])

    elif gene_space == 'avg-mlp-latent-space':

        expr = [os.path.join(path, 'mlp_latent_space') for path in expr]
        expr = get_latent_spaces(expr = expr,
                                 ids = range(1, n_latent_spaces + 1))

    elif gene_space == 'vae-latent-space':
        expr = list(expr)
        expr_files = dict(
            human = 'human_vae_embedding.csv',
            mouse = 'mouse_vae_embedding.csv'
        )
        for i, path in enumerate(expr):
            expr[i] = os.path.join(path, 'vae_latent_space', expr_files[species[i]])
        expr = [tuple(expr)]

    else:
        raise ValueError("`gene_space` must be one of "
                         "{'homologous-genes', "
                         "'mlp-latent-space', "
                         "'avg-mlp-latent-space', "
                         "'vae-latent-space'}")

    # Expand all combinations of images and expression spaces
    inputs = list(product(imgs, expr))

    iterator = partial(_compute_transcriptomic_similarity,
                       species = species,
                       masks = masks,
                       microarray_coords = microarray_coords,
                       signed = signed,
                       metric = metric,
                       threshold = threshold,
                       threshold_value = threshold_value,
                       threshold_symmetric = threshold_symmetric,
                       return_signed = return_signed)

    if execution == 'local':
        client = Client(processes = True,
                        n_workers = nproc,
                        threads_per_worker = 1)
    elif execution == 'slurm':
        cluster = SLURMCluster(
            cores = slurm_kwargs.get('cores', 1),
            memory = slurm_kwargs.get('memory', '16GB'),
            walltime = slurm_kwargs.get('walltime', '08:00:00')
        )
        cluster.scale(jobs = slurm_kwargs.get("jobs", nproc))
        client = Client(cluster)
    else:
        raise ValueError("Argument `execution` must be one of {'local', 'slurm'}")

    futures = client.map(iterator, inputs)

    if verbose:
        progress(futures)

    sim = client.gather(futures)

    if return_signed:
        results = pd.DataFrame(
            dict(img1 = [x[0][0] for x in inputs],
                 img2 = [x[0][1] for x in inputs],
                 similarity_pos = [x[0] for x in sim],
                 similarity_neg = [x[1] for x in sim])
        )
    else:
        results = pd.DataFrame(
            dict(img1 = [x[0][0] for x in inputs],
                 img2 = [x[0][1] for x in inputs],
                 similarity = sim)
        )

    if gene_space == 'avg-mlp-latent-space':
        if return_signed:
            raise Exception
        else:
            results = (results
                       .groupby(by = ['img1', 'img2'],
                                as_index = False)
                       .mean()
                       .copy())

    return results
