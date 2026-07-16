# Apptainer tools for miaoseq external backends

Provides `vsearch`, `abpoa`, `blastn`/`makeblastdb`, and `mmseqs` via Apptainer biocontainers.

- **Assemble (mainline):** `vsearch` (clustering) + `abpoa` (consensus)
- **Reassess (optional):** `blastn`, `mmseqs` (pairwise comparison)

## Layout

- `sif/` — pulled SIF images
- `bin/` — thin wrappers that `apptainer exec` the matching SIF

## Setup

```bash
# pull (or re-pull) SIFs
./cursor_dev/containers/pull_sifs.sh

# put wrappers first on PATH for R / shell sessions
export PATH="/home/ftom/01_wd/softDevel/miaoseq/cursor_dev/containers/bin:$PATH"
```

## Images

| Tool | Image |
|------|--------|
| VSEARCH | `quay.io/biocontainers/vsearch:2.31.0--hd2be7a0_0` |
| abPOA | `quay.io/biocontainers/abpoa:1.5.6--hb7acf71_0` |
| BLAST+ | `quay.io/biocontainers/blast:2.16.0--h66d330f_5` |
| MMseqs2 | `quay.io/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_3` |
