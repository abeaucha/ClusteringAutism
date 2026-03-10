#!/bin/bash

preprocess_human_images.py \
--imgdir data/human/registration/v3/jacobians/ \
--jacobians absolute \
--nproc 8
