import contextlib
import os
import sys
import utils
import multiprocessing as mp
import numpy as np
import pandas as pd
import processing
from datatable import fread
from functools import partial
from io import StringIO
from itertools import product
from pyminc.volumes.factory import volumeFromFile
from tqdm import tqdm


@contextlib.contextmanager
def silence():
    save_stdout = sys.stdout
    sys.stdout = StringIO()
    yield
    sys.stdout = save_stdout


def fetch_microarray_coordinates(metadata, outfile, labels = True,
                                 verbose = False):
    script_args = locals().copy()
    script_args['labels'] = 'true' if script_args['labels'] else 'false'
    script_args['verbose'] = 'true' if script_args['verbose'] else 'false'
    script = 'fetch_microarray_coordinates.R'
    utils.execute_R(script = script, **script_args)
    return outfile


def coordinates_to_minc(coordinates, template, outfile, type = 'labels',
                        verbose = False):
    script_args = locals().copy()
    script_args['verbose'] = 'true' if script_args['verbose'] else 'false'
    script = 'coordinates_to_minc.R'
    utils.execute_R(script = script, **script_args)
    return outfile


def prepare_microarray_coordinates(metadata, transforms,
                                   outdir = 'data/human/expression/',
                                   annotations = None):
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
        expr = fread(expr, header = True).to_pandas()
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
        expr = fread(expr, header = True).to_pandas()

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

    # Create a mask for voxels with microarray date
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


def compute_transcriptomic_similarity(imgs, expr, masks, microarray_coords,
                                      metric = 'correlation', signed = True,
                                      threshold = 'top_n',
                                      threshold_value = 0.2,
                                      threshold_symmetric = True):
    """
    Compute the similarity between a human image and a mouse image. 
    
    """

    img_human, img_mouse = imgs
    expr_human, expr_mouse = expr
    mask_human, mask_mouse = masks

    human = human_signature(img = img_human,
                            expr = expr_human,
                            mask = mask_human,
                            coords = microarray_coords,
                            signed = signed,
                            threshold = threshold,
                            threshold_value = threshold_value,
                            threshold_symmetric = threshold_symmetric)
    
    mouse = mouse_signature(img = img_mouse,
                            expr = expr_mouse,
                            mask = mask_mouse,
                            signed = signed,
                            threshold = threshold,
                            threshold_value = threshold_value,
                            threshold_symmetric = threshold_symmetric)

    if signed:
        sim = []
        for signatures in zip(human, mouse):
            sim.append(similarity(x = signatures[0],
                                  y = signatures[1],
                                  metric = metric))

        is_nan = np.isnan(sim)
        if all(is_nan):
            raise Exception("No surviving voxels.")
        sim = np.array(sim)
        sim = sim[~is_nan]
        sim = np.mean(sim)

    else:
        sim = similarity(x = human,
                         y = mouse,
                         metric = metric)

    return sim


def get_latent_spaces(expr, ids = None):
    latent_spaces = []
    for e in expr:
        e = os.path.join(e, 'latent_space')
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


def tempfunc(inputs, **kwargs):
    return compute_transcriptomic_similarity(imgs = inputs[0],
                                             expr = inputs[1],
                                             **kwargs)


def transcriptomic_similarity(imgs, expr, masks, microarray_coords,
                              gene_space = 'average-latent-space',
                              n_latent_spaces = 100, latent_space_id = 1,
                              metric = 'correlation', signed = True,
                              threshold = 'top_n', threshold_value = 0.2,
                              threshold_symmetric = True, parallel = False,
                              nproc = None):

    """
    Compute the transcriptomic similarity for a set of mouse and human images.

    Parameters
    ----------
    imgs: tuple of str or list of tuple of str
        A tuple of length 2 containing the paths to the human and mouse images
        (.mnc) to compare. Multiple pairs of human and mouse images can be
        passed as a list of tuples.
    expr: tuple of str
        A tuple of length 2 containing the paths to the human and mouse
        expression directories.
    masks: tuple of str
        A tuple of length 2 containing the paths to the human and mouse mask
        images (.mnc).
    microarray_coords: str
        The path to the human microarray sample coordinates file (.csv).
    gene_space: {'average-latent-space', 'latent-space', 'homologous-genes'}
        The transcriptomic common space to use for comparison.
    n_latent_spaces: int, default 100
        The number of latent spaces to aggregate when `gene_space` =
        'average-latent-space'. Ignored otherwise.
    latent_space_id: int, default 1
        The ID of the latent space to use when 'gene_space' = 'latent-space'.
        Ignored otherwise.
    metric: str, default 'correlation'
        The metric used to compute the similarity of mouse and human images.
    signed: bool, default True
    threshold: {'top_n', 'intensity', None}
    threshold_value: float, default 0.2
    threshold_symmetric: bool, default True
    parallel: bool, default True
    nproc: int, default None

    Returns
    -------
    """

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

        expr = list(expr)
        expr_files = (
            'HumanExpressionMatrix_samples_pipeline_abagen_homologs_scaled.csv'
            'MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv',
        )
        for i, path in enumerate(expr):
            expr[i] = os.path.join(path, 'input_space', expr_files[i])
        expr = [tuple(expr)]

    elif gene_space == 'latent-space':

        expr = get_latent_spaces(expr = expr, ids = [latent_space_id])

    elif gene_space == 'average-latent-space':

        expr = get_latent_spaces(expr = expr,
                                 ids = range(1, n_latent_spaces + 1))

    else:
        raise ValueError("`gene_space` must be one of "
                         "{'homologous-genes', " 
                         "'latent-space', " 
                         "'average-latent-space'}")
        
    inputs = list(product(imgs, expr))
    
    tempfunc_partial = partial(tempfunc,
                               masks = masks,
                               microarray_coords = microarray_coords,
                               signed = signed,
                               metric = metric,
                               threshold = threshold,
                               threshold_value = threshold_value,
                               threshold_symmetric = threshold_symmetric)

    if parallel:
        if nproc is None:
            raise ValueError("Set the nproc argument to specify the "
                             "number of processors to use in parallel.")
        pool = mp.Pool(nproc)
        sim = []
        for s in tqdm(pool.imap(tempfunc_partial, inputs), total = len(inputs)):
            sim.append(s)
        pool.close()
        pool.join()
    else:
        sim = list(map(tempfunc_partial, tqdm(inputs)))

    out = pd.DataFrame(
        dict(human_img = [x[0][0] for x in inputs],
             mouse_img = [x[0][1] for x in inputs],
             similarity = sim)
    )

    if gene_space == 'average-latent-space':
        out = out.groupby(by = ['human_img','mouse_img'], as_index = False).mean().copy()

    return out