"""
This part of the workflow constructs the phylogenetic tree.

REQUIRED INPUTS:

    metadata            = ../data/curated_metadata/{species}_{segment}.tsv
    prepared_sequences  = results/{species}/{segment}/aligned.fasta

OUTPUTS:

    tree            = results/trees/{species}_{segment}_tree.nwk
    branch_lengths  = results/{species}/{species}_{segment}_branch_lengths.json

This part of the workflow usually includes the following steps:

    - augur tree
    - augur refine

See Augur's usage docs for these commands for more details.
"""
rule tree:
    """Building tree"""
    input:
        alignment = rules.align.output
    output:
        tree = "results/{species}/{segment}/prelim_tree.nwk"
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method}
        """

rule refine:
    """
    Rerooting tree at midpoint
    """
    input:
        tree = rules.tree.output,
        alignment = rules.align.output,
        metadata = "data/{species}/{segment}/metadata.tsv",
    output:
        tree = "results/{species}/{segment}/final_tree.nwk",
        node_data = "results/{species}/{segment}/branch_lengths.json"
    params:
        strain_id_field = config["strain_id_field"],
        root = config['refine']['root'],
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root}
        """
