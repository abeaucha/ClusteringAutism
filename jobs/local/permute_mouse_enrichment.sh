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
       --equal-clusters true \
       --output-dir "$PROJECTPATH/data/mouse/derivatives/v3/107/enrichment/StringDB_12.0_Bader_2025/950/permutations_equal_clusters/" \
       --n-permutations "$batch_size" \
       --n-workers 32

done

