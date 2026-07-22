"""
This part of the workflow creates additonal annotations for the phylogenetic tree.

REQUIRED INPUTS:

    metadata            = ../data/curated_metadata/{species}_{segment}.tsv
    prepared_sequences  = results/{species}/{segment}/aligned.fasta
    tree                = results/trees/{species}_{segment}_tree.nwk

OUTPUTS:

    node_data = results/{species}/{segment}/*.json

    There are no required outputs for this part of the workflow as it depends
    on which annotations are created. All outputs are expected to be node data
    JSON files that can be fed into `augur export`.

    See Nextstrain's data format docs for more details on node data JSONs:
    https://docs.nextstrain.org/page/reference/data-formats.html

This part of the workflow usually includes the following steps:

    - augur traits
    - augur ancestral
    - augur translate
    - augur clades

See Augur's usage docs for these commands for more details.

Custom node data files can also be produced by build-specific scripts in addition
to the ones produced by Augur commands.
"""
rule ancestral:
    """Reconstructing ancestral sequences and mutations"""
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment,
    output:
        node_data = "results/{species}/{segment}/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    """Translating amino acid sequences"""
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = "../shared/{species}_{segment}_refseq.gb"
    output:
        node_data = "results/{species}/{segment}/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """
