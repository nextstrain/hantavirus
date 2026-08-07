# Phylogeny

This workflow builds a phylogenetic tree for each species/segment from the curated output of
the [ingest](../ingest) workflow and exports [Auspice](https://docs.nextstrain.org/projects/auspice/)
datasets for visualization.

## Prerequisites

Run the [ingest](../ingest) workflow first — this workflow reads its curated sequences and
metadata from `../ingest/results/`. It also uses the per-species/segment reference GenBank
files in [`../shared/`](../shared).

## Running

From within this directory:

```bash
snakemake --cores 4
```

See the [top-level README](../README.md) for environment setup.

## Steps

The workflow (`snakefile`, including `rules/*.smk`) runs:

1. **Prepare** (`rules/prepare_sequences.smk`) — `augur filter` (per-species/segment length
   bounds and exclusion list) then `augur align` against the segment's RefSeq.
2. **Construct** (`rules/construct_phylogeny.smk`) — `augur tree` (IQ-TREE) then
   `augur refine` (midpoint rooting by default).
3. **Annotate** (`rules/annotate_phylogeny.smk`) — `augur ancestral` (nucleotide mutations)
   and `augur translate` (amino-acid mutations, using the RefSeq gene annotations).
4. **Export** (`rules/export.smk`) — assigns trait colors (`scripts/assign-colors.py`) and
   runs `augur export v2` to produce the final Auspice dataset.

## Outputs

One Auspice dataset per species and segment:

```
auspice/{species}_{segment}.json
```

View them with:

```bash
auspice view --datasetDir auspice
```

## Configuration

- `defaults/config.yaml` — filter length bounds (nested per species then segment), rerooting
  strategy, and the metadata columns exported to Auspice.
- `defaults/auspice_config.json` — Auspice display settings (colorings, filters, panels,
  maintainers).
- `defaults/description.md` — footer text shown in Auspice.
- `defaults/color_schemes.tsv` / `defaults/color_orderings.tsv` — color palette and trait
  value ordering.
- `defaults/exclude.txt` — strains to drop during filtering.

The species and segments built are set at the top of `snakefile`.
