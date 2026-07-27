# Ingest

This workflow downloads Orthohantavirus sequences and metadata from NCBI, classifies each
sequence by species and genome segment, and curates the metadata into a standardized form
for the [phylogeny](../phylogeny) workflow.

## Running

From within this directory (paths in the rules are relative to it):

```bash
snakemake --cores 4
```

See the [top-level README](../README.md) for environment setup.

## Tests

Unit tests for the helper scripts live in `tests/` (stdlib `unittest`):

```bash
python -m unittest discover -s tests
```

## Steps

The workflow (`snakefile`, including `rules/*.smk`) runs in three stages:

1. **Fetch** (`rules/fetch_from_ncbi.smk`) — downloads the full NCBI Virus dataset package
   for the configured taxon (`ncbi_taxon_id` in `defaults/config.yaml`) with the NCBI
   Datasets CLI, extracts sequences and a metadata report, and pulls the per-species RefSeq
   GenBank records.
2. **Classify** (`rules/classify_sequences.smk`) — builds a DIAMOND protein database from the
   RefSeq proteins, aligns every downloaded sequence against it with `diamond blastx`, and
   assigns each sequence a species and segment from its best hit.
3. **Curate** (`rules/curate.smk`) — filters sequences by percent identity to the reference
   (`PIDENT` in `snakefile`), then runs the `augur curate` chain (field renaming, string
   normalization, strain-name/date standardization, GenBank location parsing, titlecasing,
   author abbreviation, and geolocation rules) driven by `defaults/config.yaml`.

## Outputs

Per species and segment, consumed by the phylogeny workflow:

```
results/{species}/{segment}/metadata_curated.tsv
results/{species}/{segment}/sequences_curated.fasta
```

## Configuration

- `defaults/config.yaml` — NCBI taxon ID, metadata field mapping, and all `augur curate`
  settings.
- `defaults/rename_species.csv` — normalizes NCBI species names to the canonical names.
- `defaults/special_cases.csv` — manual overrides for protein extraction / classification.
- `defaults/geolocation_rules.tsv` — geographic name corrections applied during curation.

The species built and the percent-identity threshold (`PIDENT`) are set at the top of
`snakefile`.
