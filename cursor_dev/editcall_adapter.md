Status: **Phase J implemented** — `doEditcall()` is **per-read only**.
There is no consensus-based editcall path. Guide-level Plan A and dual-cut Plan A′ are supported.

Policy: [pipeline_revise.md](pipeline_revise.md) §7.3–7.4, §17, §18.

## Two paths after gene assignment

| Path | Purpose | Function |
|------|---------|----------|
| **B (primary)** | Edit frequency near PAM | `doEditcall(gene_assign, …)` + FASTQ / `by_sample` |
| **A (optional)** | Full-length amplicon consensus / QC | `doAssembleAmplicons` → `amplicon/` |

Assemble and Editcall do not share clustering settings (§15.2).

## API

```r
edit <- doEditcall(
  gene_assign = ga,              # doAssignGenes() result or amplicon_assign/
  pam_list = pam_list,
  genome_fn = genome_fn,
  amplicon_fn = amplicon_fn,     # prepAmpliconDB() output
  primer_list = primer_list,
  editcall_dir = editcall_dir,
  sample_list = sample_list,     # optional plate layout
  check_window = 10L,
  anchor_bp = 5L,                # §15.7 adaptive window anchors
  max_expand = 50L,              # discard if anchors not found
  min_span_bp = 30L,             # Plan A′: min expected cut–cut span
  excision_tol_bp = 20L,         # Plan A′: |del_span - expected| tolerance
  min_count = 5L,
  fastq = fastq,                 # or sample_fastq_dir
  n_core = 1L
)
```

Internal steps (current):

1. Resolve `gene_assignments` from `gene_assign` (or fallback `amplicon_out$gene_assignments`).
2. Load read sequences from `fastq` or `sample_fastq_dir`.
3. Map PAM/cut from genome → amplicon → `ref_insert` (`matchPattern`; optional strand offset on genome).
4. Primer-trim once; align insert↔ref_insert once; extract §15.7 window **per PAM / `target_gene`**.
5. If ≥2 PAMs on the gene: classify adjacent pairs (Plan A′); on `both_cut_excision`, write junction alleles into Plan A rows (`allele_source = "excision"`, `intact = FALSE`).
6. Exact-match aggregate `sample × target_gene × read_seq`; apply `min_count` and legacy top-fraction rule.
7. Write `editcall_joint.csv` / `editcall_joint_summary.csv` alongside guide-level tables.
8. Safety-net: drop local `read_seq` longer than `2 * check_window + 2 * max_expand + 20` (excision alleles exempt).

Assemble purity (`doAssembleAmplicons(min_cluster_purity = 0.8)`) applies only to Path A outputs.

## `pam_list`

CSV, no header:

| Col | Meaning |
|-----|---------|
| 1 | gene ID — must match `primer_list` / `amplicon.fa` gene names |
| 2 | chromosome / seqname — **exact match** to `genome_fn` FASTA headers (no `chr%02d` rewriting) |
| 3 | PAM start (1-based) |
| 4 | optional guide ID (required when multiple rows share a gene) |
| 5 | optional strand `+`/`-` (`+` → cut = pam−3, `-` → cut = pam+3 on genome; missing → cut = pam start) |

Mismatch of col 2 vs genome → **stop**. See §7.3 / §18.

## Multi-guide on one amplicon

| Case | How to declare | Tool behavior |
|------|----------------|---------------|
| **1** Separate primer pairs / amplicons | Distinct gene IDs (`SD1_1`, `SD1_2`) | Normal per-gene editcall |
| **2** Same primers / same amplicon, multiple PAMs | Same gene ID + multiple `pam_list` rows with distinct guides | **Plan A** per-guide windows + **Plan A′** joint both-cut / inter-cut excision |

- Plan A: one global insert↔ref alignment per read; extract §15.7 window once per PAM; aggregate by `target_gene = gene_guide`.
- Plan A′: for adjacent guide pairs (sorted by insert cut, then `guide_id`), classify `event_class` (`wt` / `g_i_only` / `g_j_only` / `both_local` / `both_cut_excision`); write `editcall_joint.csv` / `editcall_joint_summary.csv`.
- `excision_rate` denominator = reads with ≥1 Plan A row for that gene.

## Reporting

- `evalMiao(out_dir)` prefers `amplicon_assign/gene_assignments.tsv`, falls back to `amplicon/`.
- Editcall-only runs (no `amplicon/`) work with `demultiplex/` + `amplicon_assign/` + `editcall/`.
- `editViewer(out_dir)` reads `editcall/editcall_summary.csv` (guide-level).
- Joint / excision plate display: deferred.

## Migration

Named `amplicon_out = amp` still resolves gene assignments if present; prefer `gene_assign = ga`.
