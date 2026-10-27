#!/bin/bash

set -Eeuo pipefail

trap 'status=$?; printf "ERROR: %s:%s: %s (exit %s)\n" \
  "${BASH_SOURCE[0]}" "$LINENO" "$BASH_COMMAND" "$status" >&2; \
  exit "$status"' ERR

batch_size=1000
B=10000

for seed_start in $(seq 1 "$batch_size" "$((B - batch_size + 1))"); do

    echo "Starting permutation batch at seed $seed_start"

   permute_mouse_cluster_enrichment.r \
       --pipeline-dir "$PROJECTPATH/data/mouse/derivatives/v3/" \
       --params-id 107 \
       --gene-score 950 \
       --seed-start "$seed_start" \
       --n-permutations 1000 \
       --n-workers 32

done

