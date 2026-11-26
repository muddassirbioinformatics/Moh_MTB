#!/usr/bin/env bash
# Load every pinned container image from MTB_DATA/containers/*.tar.gz into the local Docker.
# Run after unpacking the data zip. If a tarball is missing, `docker pull` the image named in conf/containers.tsv
# (mtb-report:1.1.0 is a local image: docker build -t mtb-report:1.1.0 containers/mtb-report).
set -euo pipefail
C="${1:-$HOME/mtbdata/MTB_DATA/containers}"
for f in "$C"/*.tar.gz; do echo "loading $f"; docker load -i "$f"; done
echo "--- images now present:"; docker images --format '{{.Repository}}:{{.Tag}}' | sort
echo "--- compare with env/docker_images_pinned.tsv (digests) if needed"
