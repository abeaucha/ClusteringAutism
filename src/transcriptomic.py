import contextlib
import os
import sys
import multiprocessing as mp
import numpy as np
import pandas as pd
from datatable import fread
from functools import partial
from io                     import StringIO
from pyminc.volumes.factory import volumeFromFile
from src import processing
from tqdm import tqdm


@contextlib.contextmanager
def silence():
    save_stdout = sys.stdout
    sys.stdout = StringIO()
    yield
    sys.stdout = save_stdout


def mouse_signature(img, expr, mask, signed = True, threshold = 'top_n',
                    threshold_value = 0.2, threshold_symmetric = True):

    """
    Compute the transcriptomic signature for a mouse image
    """

    #Import image
    img = processing.import_image(img = img, mask = mask, flatten = False)

    #Import brain mask
    mask = processing.import_image(img = mask, flatten = True)

    #Import transcriptomic data
    with silence():
        expr = fread(expr, header = True).to_pandas()
    if np.sum(mask) != expr.shape[0]:
        raise Exception("The number of voxels in mask does not match the "
                        "number of rows in expr.")

    #Threshold the image if desired
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

    #Create a mask from the (thresholded) image
    img_mask = processing.mask_from_image(img = img, signed = signed)
    img_mask = img_mask.flatten()
    img_mask = img_mask[mask == 1]

    #If signed, return signatures for positive and negative image values
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

    #Import image
    img = processing.import_image(img = img, mask = mask, flatten = False)

    #Import transcriptomic data
    with silence():
        expr = fread(expr, header = True).to_pandas()

    #Import microarray coordinates
    coords = pd.read_csv(coords)

    #Threshold the image if desired
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

    #Create a mask from the (thresholded) image
    img_mask = processing.mask_from_image(img = img, signed = signed)

    #Create a mask for voxels with microarray date
    vol = volumeFromFile(mask)
    microarray_mask = np.zeros(coords.shape[0], dtype = int)
    for i, row in coords.iterrows():
        coords_world = np.array([row['x'], row['y'], row['z']])
        coords_voxel = vol.convertWorldToVoxel(coords_world)
        coords_voxel = np.round(coords_voxel)
        x = int(coords_voxel[0])-1
        y = int(coords_voxel[1])-1
        z = int(coords_voxel[2])-1
        microarray_mask[i] = img_mask[x,y,z]
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

    if (len(x) != len(y)):
        raise Exception("x and y must have the same length")

    if metric == 'correlation':
        d = np.corrcoef(x = x, y = y)
        d = d[0,1]
    elif metric == 'cosine':
        d = np.dot(x, y)/np.sqrt((np.dot(x, x)*np.dot(y,y)))
    elif metric == 'euclidean':
        d = np.linalg.norm(x-y, ord = 2)
    elif metric == 'manhattan':
        d = np.sqrt(np.linalg.norm(x-y, ord = 1))
    elif metric == 'KL-divergence':
        d = np.sum(np.where(x != 0, x * np.log(x / y), 0))
    else:
        raise ValueError("Invalid metric: {}".format(metric))

    return d


def compute_transcriptomic_similarity(expr, imgs, masks, microarray_coords,
                                      metric = 'correlation', signed = True,
                                      threshold = 'top_n', threshold_value = 0.2,
                                      threshold_symmetric = True):

    """
    Compute the similarity between a mouse image and a human image. 
    
    """

    img_mouse, img_human = imgs
    expr_mouse, expr_human = expr
    mask_mouse, mask_human = masks

    mouse = mouse_signature(img = img_mouse,
                            expr = expr_mouse,
                            mask = mask_mouse,
                            signed = signed,
                            threshold = threshold,
                            threshold_value = threshold_value,
                            threshold_symmetric = threshold_symmetric)

    human = human_signature(img = img_human,
                            expr = expr_human,
                            mask = mask_human,
                            coords = microarray_coords,
                            signed = signed,
                            threshold = threshold,
                            threshold_value = threshold_value,
                            threshold_symmetric = threshold_symmetric)

    if signed:
        sim = []
        for signatures in zip(mouse, human):
            sim.append(similarity(x = signatures[0],
                                  y = signatures[1],
                                  metric = metric))
        sim = np.mean(np.array(sim))

    else:
        sim = similarity(x = mouse,
                         y = human,
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
            ls = [os.path.join(e, file) for i, file in ids_files_all if i in ids]
        else:
            ls = [os.path.join(e, file) for i, file in ids_files_all]
        latent_spaces.append(ls)
    latent_spaces = list(zip(latent_spaces[0], latent_spaces[1]))

    return latent_spaces


def transcriptomic_similarity(imgs, expr, masks, microarray_coords,
                              gene_space = 'average-latent-space',
                              n_latent_spaces = 100, latent_space_id = 1,
                              metric = 'correlation', signed = True,
                              threshold = 'top_n', threshold_value = 0.2,
                              threshold_symmetric = True, show_progress = True):

    if gene_space == 'homologous-genes':

        expr = [i for i in expr]
        expr_files = ('MouseExpressionMatrix_voxel_coronal_log2_grouped_imputed_homologs_scaled.csv',
                      'HumanExpressionMatrix_samples_pipeline_abagen_homologs_scaled.csv')
        for i, path in enumerate(expr):
            expr[i] = os.path.join(path, 'input_space', expr_files[i])
        expr = (expr[0], expr[1])

        sim = compute_transcriptomic_similarity(expr = expr,
                                                imgs = imgs,
                                                masks = masks,
                                                microarray_coords = microarray_coords,
                                                signed = signed,
                                                metric = metric,
                                                threshold = threshold,
                                                threshold_value = threshold_value,
                                                threshold_symmetric = threshold_symmetric)

    elif gene_space == 'latent-space':

        expr = get_latent_spaces(expr = expr, ids = [latent_space_id])

        sim = compute_transcriptomic_similarity(expr = expr[0],
                                                imgs = imgs,
                                                masks = masks,
                                                microarray_coords = microarray_coords,
                                                signed = signed,
                                                metric = metric,
                                                threshold = threshold,
                                                threshold_value = threshold_value,
                                                threshold_symmetric = threshold_symmetric)

    elif gene_space == 'average-latent-space':

        expr = get_latent_spaces(expr = expr,
                                 ids = range(1, n_latent_spaces+1))
        if show_progress: expr = tqdm(expr)

        sim = []
        for e in expr:
            s = compute_transcriptomic_similarity(imgs = imgs,
                                                    expr = e,
                                                    masks = masks,
                                                    microarray_coords = microarray_coords,
                                                    signed = signed,
                                                    metric = metric,
                                                    threshold = threshold,
                                                    threshold_value = threshold_value,
                                                    threshold_symmetric = threshold_symmetric)
            sim.append(s)
        sim = np.mean(sim)

    else:
        raise ValueError("Argument gene_space must be one of "
                         "{'homologous-genes', 'latent-space', 'average-latent-space'}")

    return sim