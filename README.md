# Nextstrain build for Orthohantavirus

This repository contains the [Nextstrain](https://nextstrain.org) build for
**Orthohantavirus** (hantavirus). It fetches genomic data from NCBI, classifies each
sequence by species and genome segment, curates the metadata, builds per-species/segment
phylogenies, and exports [Auspice](https://docs.nextstrain.org/projects/auspice/) datasets
for interactive visualization.

Hantaviruses have a tri-segmented negative-sense RNA genome — the **L** (large), **M**
(medium), and **S** (small) segments — so builds are produced independently for each
species/segment combination.

## Repository organization

The build is split into two independent [Snakemake](https://snakemake.readthedocs.io/)
workflows, following the standard Nextstrain pathogen-repo layout:

| Directory    | Purpose                                                                               |
| ------------ | ------------------------------------------------------------------------------------- |
| `ingest/`    | Fetch sequences from NCBI, classify them by species/segment, and curate metadata.     |
| `phylogeny/` | Align sequences, build and refine trees, annotate, and export Auspice JSONs.          |
| `shared/`    | Per-species/segment RefSeq GenBank files used as alignment and annotation references. |

The two workflows are connected only through the filesystem: `phylogeny/` reads the curated
output written to `ingest/results/`. There is no top-level Snakefile — each stage is run
separately.

## Installation

The workflows depend on [Augur](https://docs.nextstrain.org/projects/augur/) and a number of
bioinformatics tools. The Nextstrain [`nextstrain-base`](https://docs.nextstrain.org/en/latest/reference/glossary.html)
conda package bundles most of them (`augur`, `snakemake`, the NCBI Datasets CLI, `csvtk`,
`iqtree`, plus `biopython`/`pandas`/`numpy`), but **not `diamond`**, which the ingest workflow
needs for sequence classification. Create an environment with both in one command:

```bash
mamba create -n nextstrain -c conda-forge -c bioconda -c nextstrain nextstrain-base diamond
```

Then activate it before running either workflow:

```bash
conda activate nextstrain
```

Alternatively, use the [Nextstrain CLI](https://docs.nextstrain.org/projects/cli/) with a
managed runtime, ensuring `diamond` is additionally available on the `PATH`.

## Usage

Each workflow uses `workdir: workflow.current_basedir`, so run it from within its own
directory (or point Snakemake at its `snakefile`); all paths inside the rules are relative to
that directory.

### 1. Ingest

```bash
cd ingest
snakemake --cores 4
```

This downloads the NCBI virus dataset for the configured taxon, assigns each sequence a
species and segment by aligning it against RefSeq proteins with DIAMOND, filters by percent
identity, and runs the `augur curate` pipeline. Curated output is written to:

```
ingest/results/{species}/{segment}/metadata_curated.tsv
ingest/results/{species}/{segment}/sequences_curated.fasta
```

### 2. Phylogeny

```bash
cd phylogeny
snakemake --cores 4
```

This filters and aligns the curated sequences, builds a tree with IQ-TREE, refines it
(midpoint rooting by default), reconstructs ancestral nucleotide and amino-acid mutations,
assigns colors, and exports one Auspice dataset per species/segment:

```
phylogeny/auspice/{species}_{segment}.json
```

Use `snakemake -n` (dry run) in either directory to preview the DAG before running.

### 3. Visualize

```bash
auspice view --datasetDir phylogeny/auspice
```

## Configuration

Most behavior is controlled through the `defaults/config.yaml` file in each workflow, without
editing rule logic:

- **`ingest/defaults/config.yaml`** — NCBI taxon ID, metadata field mapping, and the `augur
  curate` settings (strain-name regex, date formats, titlecasing, geolocation rules).
- **`phylogeny/defaults/config.yaml`** — per-species/segment sequence length bounds, rerooting
  strategy, and the metadata columns exported to Auspice.

The set of **species and segments** to build is defined as Python lists at the top of each
workflow's `snakefile` (`ingest/snakefile`, `phylogeny/snakefile`) — edit those to add or
remove a build. Other tunable inputs include `ingest/defaults/rename_species.csv` (species
name normalization), `ingest/defaults/special_cases.csv` (manual classification overrides),
and `phylogeny/defaults/exclude.txt` (strains to drop).

## References

Corresponding RefSeq references for each species/segment live in `shared/` as
`{species}_{segment}_refseq.gb`. They serve as both the alignment reference and the
gene-annotation source for amino-acid translation.
