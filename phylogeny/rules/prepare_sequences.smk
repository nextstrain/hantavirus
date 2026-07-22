"""
This part of the workflow prepares sequences for constructing the phylogenetic tree.

REQUIRED INPUTS:

    metadata    = data/curated_metadata/{species}_{segment}.tsv
    sequences   = data/curated_sequences/{species}_{segment}.fasta
    reference   = data/shared/{species}_{segment}_refseq.fasta

OUTPUTS:

    prepared_sequences = results/prepared_sequences.fasta

This part of the workflow usually includes the following steps:

    - augur index
    - augur filter
    - augur align
    - augur mask

See Augur's usage docs for these commands for more details.
"""

rule filter:
    """
    Filtering to
      - sequences per {params.species} and {params.segment}
      - excluding strains in {input.exclude}
    """
    input:
        sequences = "../ingest/results/{species}/{segment}/sequences_curated.fasta",
        metadata = "../ingest/results/{species}/{segment}/metadata_curated.tsv",
        exclude = config['filter']['exclude']
    output:
        sequences = "results/{species}/{segment}/filtered.fasta"
    params:
        strain_id_field = config["strain_id_field"],
        min_length = lambda w: config['filter']['min_length'][w.species][w.segment],
        max_length = lambda w: config['filter']['max_length'][w.species][w.segment],
        exclude = config['filter']['exclude'],
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-sequences {output.sequences} \
            --min-length {params.min_length} \
            --max-len {params.max_length} \
            --exclude {input.exclude} 
        """

from Bio import SeqIO

def find_reference_name(species, segment):
    record = SeqIO.read(f"../shared/{species}_{segment}_refseq.gb", "genbank")
    reference, version = record.id.split('.')

    return reference

rule align:
    """
    Aligning sequences to {params.reference}
    - filling gaps with -
    """
    input:
        sequences = "results/{species}/{segment}/filtered.fasta",
    params:
        reference = lambda wildcards: find_reference_name(wildcards.species, wildcards.segment)
    output:
        alignment = "results/{species}/{segment}/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-name {params.reference} \
            --output {output.alignment} \
            --remove-reference
        """
