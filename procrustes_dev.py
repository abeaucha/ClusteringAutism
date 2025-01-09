import os
from processing import import_images
import numpy as np

cv_dir = 'data/human/derivatives/v3/916/cross_validation/'

target_dir = 'data/human/derivatives/v3/916/cross_validation/sample_1/centroids/resolution_0.8/'

jacobians = ('absolute', 'relative')

# for j in jacobians:
jacobian = 'absolute'

# Import centroids 1 and 2 for given jacobians
target_dir_j = os.path.join(target_dir, jacobian, '')
target_imgs = ['centroid_nk_2_k_1.mnc', 'centroid_nk_2_k_2.mnc']
target_imgs = [os.path.join(target_dir_j, file) for file in target_imgs]

mask = 'data/human/registration/v3/reference_files/mask_0.8mm.mnc'

target_matrix = import_images(imgfiles = target_imgs, mask = mask,
                              output_format = 'numpy')

input_dir = 'data/human/derivatives/v3/916/cross_validation/'
nsamples = len(os.listdir(input_dir))

# for i in range(2, nsamples+1):
for i in range(2, 10):

    # Get input dir for sample i and Jacobians j
    input_dir_i = os.path.join(input_dir, 'sample_{}'.format(i))
    input_dir_i = os.path.join(input_dir_i, 'centroids', 'resolution_0.8', jacobian)

    # Get centroids 1 and 2
    input_imgs = ['centroid_nk_2_k_1.mnc', 'centroid_nk_2_k_2.mnc']
    input_imgs = [os.path.join(input_dir_i, file) for file in input_imgs]

    input_matrix = import_images(imgfiles = input_imgs, mask = mask,
                                 output_format = 'numpy')

    # from scipy.linalg import orthogonal_procrustes
    #
    # R, scale = orthogonal_procrustes(A = input_matrix,
    #                                  B = target_matrix)
    #
    # from scipy.spatial import procrustes
    #
    # mtx1, mtx2, disparity = procrustes(data1 = target_matrix,
    #                                    data2 = input_matrix)

    cor_matrix = np.corrcoef(target_matrix, input_matrix)[:2, 2:]

    for j in range(len(input_matrix)):
        idx_max = np.where(cor_matrix[:, j] == np.max(cor_matrix[:, j]))[0][0]
        print("Sample {} image {} best match to target {}. Correlation: {:.3f}"
              .format(i, j+1, idx_max+1, np.max(cor_matrix[:, j])))



