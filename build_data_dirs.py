#!/usr/bin/env python3
# ----------------------------------------------------------------------------
# build_data_dirs.R
# Author: Antoine Beauchamp
# Created: March 9th, 2026

# Packages --------------------------------------------------------------------

import os


# Environment Variables ------------------------------------------------------

PROJECTPATH = os.getenv("PROJECTPATH")


# Main ------------------------------------------------------------------------

# Define the directories to create
directories = [
        "data/mouse",
        "data/mouse/derivatives",
        "data/mouse/expression",
        "data/mouse/registration",
        "data/mouse/atlas",
        "data/mouse/registration/jacobians",
        "data/mouse/registration/reference_files",
        "data/mouse/registration/transforms",
        "data/mouse/registration/resources",
        "data/human",
        "data/human/derivatives",
        "data/human/expression",
        "data/human/registration",
        "data/human/registration/jacobians",
        "data/human/registration/reference_files",
        "data/human/registration/subject_info",
        "data/cross_species",
        "data/enrichment"
]

# Create each directory
for directory in directories:
    os.makedirs(os.path.join(PROJECTPATH, directory), exist_ok=True)
    print(f"Created directory: {directory}")