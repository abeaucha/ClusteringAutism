#!/bin/bash

set -Eeuo pipefail

trap 'status=$?; printf "ERROR: %s:%s: %s (exit %s)\n" \
  "${BASH_SOURCE[0]}" "$LINENO" "$BASH_COMMAND" "$status" >&2; \
  exit "$status"' ERR

human_cluster_enrichment.R \
  --pipeline-dir data/human/derivatives/v3/ \
  --params-id 700 \
  --gene-score 950 \
  --stringdb-version 12.0 \
  --bader-version 2025