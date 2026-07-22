
rule add_missing_translations:
    input:
        gb_records = "data/refseqs.gb"
    output:
        all_gb_records = "data/complete_refseqs.gb"
    shell:
        """
        python scripts/add_missing_translations.py \
        --gb_file {input.gb_records} \
        --output {output.all_gb_records}
        """

rule extract_proteins:
    input: 
        all_gb_records = "data/complete_refseqs.gb",
        special_cases = "defaults/special_cases.csv"
    output:
        proteins = "data/refseq_proteins.fasta",
    shell:
        """
        python scripts/extract_proteins.py \
        --input {input.all_gb_records} \
        --special_cases {input.special_cases} \
        --output {output.proteins}
        """

rule diamond_align:
    input:
        fasta_seq="data/ncbi_dataset_sequences.fasta",
        diamond_db= "data/refseq_proteins.fasta",
    output:
        diamond_tsv="results/diamond_alignments.tsv",
    threads: 8
    shell:
        """
        diamond blastx \
            --query {input.fasta_seq} \
            --threads {threads} \
            --max-target-seqs 1 \
            --db {input.diamond_db} \
            --out {output.diamond_tsv} \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        """

rule add_to_metadata:
    input:
        diamond_tsv="results/diamond_alignments.tsv",
        metadata="data/ncbi_dataset_report.tsv",
        rename_species=config["species_rename_file"],
    output:
        "data/ncbi_dataset_seq_assigned.tsv",
    shell:
        """
        python scripts/add_to_metadata.py \
        --diamond_tsv {input.diamond_tsv} \
        --rename_species {input.rename_species} \
        --metadata {input.metadata} \
        --output {output}
        """
