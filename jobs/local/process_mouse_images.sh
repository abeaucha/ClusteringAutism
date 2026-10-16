#!/bin/bash

set -Eeuo pipefail

trap 'status=$?; printf "ERROR: %s:%s: %s (exit %s)\n" \
  "${BASH_SOURCE[0]}" "$LINENO" "$BASH_COMMAND" "$status" >&2; \
  exit "$status"' ERR

process_mouse_images.py \
  --pipeline-dir data/mouse/derivatives/v3 \
  --params-id 107 \
  --nproc 16

mouse_cluster_enrichment.R \
  --pipeline-dir data/mouse/derivatives/v3/ \
  --params-id 107 \
  --gene-score 950 \
  --stringdb-version 12.0 \
  --bader-version 2025
