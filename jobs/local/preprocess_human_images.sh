#!/bin/bash

set -Eeuo pipefail

trap 'status=$?; printf "ERROR: %s:%s: %s (exit %s)\n" \
  "${BASH_SOURCE[0]}" "$LINENO" "$BASH_COMMAND" "$status" >&2; \
  exit "$status"' ERR

preprocess_human_images.py \
--imgdir data/human/registration/v3/jacobians/ \
--jacobians absolute \
--nproc 8
