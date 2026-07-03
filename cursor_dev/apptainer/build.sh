#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

IMAGE_DIR="${IMAGE_DIR:-$ROOT/images}"
mkdir -p "$IMAGE_DIR"

DEF="$ROOT/miaoseq-refless.def"
SIF="$IMAGE_DIR/miaoseq-refless.sif"
SANDBOX="$IMAGE_DIR/miaoseq-refless.sandbox"

if command -v apptainer >/dev/null 2>&1; then
    RUNNER=apptainer
elif command -v singularity >/dev/null 2>&1; then
    RUNNER=singularity
else
    echo "apptainer or singularity is required" >&2
    exit 1
fi

echo "Building sandbox: $SANDBOX"
rm -rf "$SANDBOX"
"$RUNNER" build --sandbox "$SANDBOX" "$DEF"

echo "Converting sandbox to SIF: $SIF"
rm -f "$SIF"
if "$RUNNER" build "$SIF" "$SANDBOX"; then
    echo "SIF build complete: $SIF"
else
    echo "SIF conversion failed; sandbox is available at $SANDBOX" >&2
fi

echo "Done."
echo "export MIAOSEQ_REFLESS_SIF=${SIF:-$SANDBOX}"
