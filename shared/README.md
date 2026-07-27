# `shared/` — reference GenBank records

Per-species, per-segment RefSeq GenBank files, named
`{species}_{segment}_refseq.gb` (e.g. `Orthohantavirus_puumalaense_l_refseq.gb`).

These are the **canonical reference sequences** for the phylogeny workflow.

They are **hand-curated and committed to git**.

## Manually corrected records

- **`Orthohantavirus_puumalaense_l_refseq.gb` (NC_005225.1, 2018)** — the source RefSeq ships
  **without any CDS/translation**. The L polymerase ORF (largest ORF, start set to its
  first Met → `37..6504`, `/locus_tag="PUUVsLgp1"`) was added **in-house** so the segment
  can be translated and displayed like the others. The added ORF was
  sanity-checked by alignment against the L polymerases of the other Orthohantavirus
  references (full-length, ~2153 aa, >85% identity to the nearest).
