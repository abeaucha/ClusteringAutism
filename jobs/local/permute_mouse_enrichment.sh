#!/bin/bash

set -Eeuo pipefail

trap 'status=$?; printf "ERROR: %s:%s: %s (exit %s)\n" \
  "${BASH_SOURCE[0]}" "$LINENO" "$BASH_COMMAND" "$status" >&2; \
  exit "$status"' ERR

permute_mouse_cluster_enrichment.r \
  --pipeline-dir $PROJECTPATH/data/mouse/derivatives/v3/ \
  --params-id 107 \
  --gene-score 950 \
  --n-permutations 16 \
  --n-workers 4