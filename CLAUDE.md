# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

This project validates SARS-CoV-2 lineage detections from wastewater sequencing data. The workflow confirms that SRA samples initially flagged as containing a variant of interest (e.g., BA.2.75.X, BF.7, BQ.1, XBB.1.5, XBB.1.9.X) actually carry lineage-defining mutations, using co-variant (covar) haplotype data.

## Running the validation

```bash
# Run the main validation script (reads samples/, covar/, lineage_def/)
python check_detections.py

# Outputs to validated/<lineage>.csv
```

## Full pipeline (upstream of check_detections.py)

The scripts require the `freyja-sc2` micromamba environment (sra-tools, minimap2, samtools, covar) and the `freyja-global` repo on disk for the reference FASTA (`Assets/NC_045512_Hu-1.fasta`) and GFF annotation.

```bash
# Step 1: Download SRA FASTQs and align to SARS-CoV-2 reference
bash scripts/fetch_and_align_sra.sh          # all CSVs in samples/
bash scripts/fetch_and_align_sra.sh --dry-run  # preview without running

# Step 2: Run covar on aligned BAMs
bash scripts/run_covar.sh                    # all BAMs in bam/
bash scripts/run_covar.sh --force            # re-run even if output exists
```

Both scripts auto-detect `WORKDIR` relative to the script location, or accept `WORKDIR=<path>` as an environment variable. Per-sample failures are logged and processing continues; check `failed_accessions.log` / `failed_covar.log`.

## Data layout

| Path | Contents |
|------|----------|
| `samples/initial_detections_<LINEAGE>.csv` | Input: `collection_date_ww`, `state`, `accession` (SRR IDs) |
| `covar/<accession>.covar.tsv` | Per-sample haplotype clusters: `nt_mutations`, `cluster_depth`, etc. |
| `lineage_def/<LINEAGE>_unique.txt` | Lineage-defining nucleotide mutations (one per line, e.g. `C2790T`) |
| `<lineage>.txt` (root-level) | All mutations for a lineage (superset of `lineage_def/`) |
| `validated/validated_detections_<LINEAGE>.csv` | Output with cluster counts per sample |

## Validation logic (`check_detections.py`)

For each sample, the script counts covar clusters that contain lineage-defining mutations (LDMs):

- **1 LDM** → tallied in `n_clusters_1_ldm` / `cluster_depth_1_ldm`
- **2+ LDMs** → tallied in `n_clusters_2plus_ldm` / `cluster_depth_2plus_ldm`

Cluster depth is the absolute read count supporting that haplotype. Samples without a `.covar.tsv` file receive zero counts.

Lineage name resolution: `BA.2.75.X` looks for `lineage_def/BA.2.75.X_unique.txt` first, then falls back to `lineage_def/BA.2.75_unique.txt` (strips trailing `.X`).

## Notebook (`validate_detects.ipynb`)

Exploratory analysis using `outbreak_data`, `outbreak_tools`, and `yaml`. It loads `lineages.yml`, queries the outbreak.info API for lineage prevalence, and computes Bayesian posterior probabilities (P(lineage | observed mutations)) for covariant clusters. Requires authentication via `authenticate_user.authenticate_new_user()`.
