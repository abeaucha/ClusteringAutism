#!.venv/bin/python3

# Packages
import os
import transcriptomic
from pyminc.volumes.factory import volumeFromFile

# Global variables
# Directories
expr_dir = 'data/human/expression'
registration_dir = 'data/human/registration/'
version = 'v2'

# Expression inputs
metadata = 'SampleInformation_pipeline_abagen.csv'
annotations = 'AHBA_microarray_sample_annotations.csv'

# Registration inputs
template = 'model_0.8mm.mnc'
transforms = ['to_target_0GenericAffine.mat', 'to_target_1Warp.nii']

if __name__ == '__main__':

    # Paths to expression inputs
    metadata = os.path.join(expr_dir, metadata)
    annotations = os.path.join(expr_dir, annotations)

    # Paths to registration inputs
    registration_dir = os.path.join(registration_dir, version)
    template = os.path.join(registration_dir, 'reference_files', template)
    transforms = [os.path.join(registration_dir, 'average_to_MNI', i) for i in
                  transforms]

    # Fetch microarray coordinates and transform to study space
    print("Fetching microarray coordinates...")
    coords = transcriptomic.prepare_microarray_coordinates(metadata = metadata,
                                                           annotations = annotations,
                                                           transforms = tuple(
                                                               transforms))

    # Rename coordinates file
    coords_new = coords.replace('.csv', '_{}.csv'.format(version))
    os.rename(coords, coords_new)
    coords = coords_new

    # Get template resolution
    vol = volumeFromFile(template)
    resolution = vol.getSeparations()
    if len(set(resolution)) == 1:
        resolution = resolution[0]
    else:
        raise Exception
    vol.closeVolume()

    # Create label image from microarray samples
    print("Creating label image...")
    labels = 'AHBA_microarray_labels_study_{}_{}mm.mnc'.format(version,
                                                               resolution)
    labels = os.path.join(expr_dir, labels)
    labels = transcriptomic.coordinates_to_minc(coordinates = coords,
                                                template = template,
                                                outfile = labels,
                                                type = 'labels')

    # Create mask image from microarray samples
    print("Creating mask image...")
    mask = 'AHBA_microarray_mask_study_{}_{}mm.mnc'.format(version, resolution)
    mask = os.path.join(expr_dir, mask)
    mask = transcriptomic.coordinates_to_minc(coordinates = coords,
                                              template = template,
                                              outfile = mask,
                                              type = 'mask')