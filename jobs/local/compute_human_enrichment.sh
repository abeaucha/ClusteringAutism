#!/bin/bash

human_cluster_enrichment.R \
  --pipeline-dir data/human/derivatives/v3/ \
  --params-id 700 \
  --gene-score 950 \
  --stringdb-version 12.0 \
  --bader-version 2025