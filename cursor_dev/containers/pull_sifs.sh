#!/usr/bin/env bash
# Pull tool SIFs into cursor_dev/containers/sif/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
SIFDIR="$ROOT/sif"
mkdir -p "$SIFDIR"

# Assemble mainline: vsearch + abpoa
# Reassess optional: blast + mmseqs2
declare -A IMAGES=(
  [vsearch]="docker://quay.io/biocontainers/vsearch:2.31.0--hd2be7a0_0"
  [blast]="docker://quay.io/biocontainers/blast:2.16.0--h66d330f_5"
  [mmseqs2]="docker://quay.io/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_3"
  [abpoa]="docker://quay.io/biocontainers/abpoa:1.5.6--hb7acf71_0"
)

for name in "${!IMAGES[@]}"; do
  echo "==> pulling ${name}: ${IMAGES[$name]}"
  apptainer pull -F "$SIFDIR/${name}.sif" "${IMAGES[$name]}"
done
ls -lh "$SIFDIR"/*.sif
