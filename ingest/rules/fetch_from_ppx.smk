rule download_ppx_seqs:
    output:
        sequences="data/{species}/{segment}/sequences.fasta",
    params:
        sequences_url=lambda w: config["ppx_fetch"]["seqs"][w.segment],
    shell:
        """
        curl -fsSL {params.sequences_url:q} -o {output.sequences}
        """


rule download_ppx_meta:
    output:
        metadata="data/{species}/{segment}/full_ppx_metadata.csv",
    params:
        metadata_url=lambda w: config["ppx_fetch"]["meta"][w.segment],
    shell:
        """
        curl -fsSL '{params.metadata_url}' \
        | csvtk mutate2 -n is_reference -e '""' > {output.metadata}
        """


rule handle_segmented_metadata_fields:
    input:
        metadata="data/{species}/{segment}/full_ppx_metadata.csv",
    output:
        subset_metadata="data/{species}/{segment}/ppx_metadata.csv",
    run:
        # 1. Throw out all columns that have the wrong segment, for `L` throw out `_M` and `_S` columns, etc.
        # 2. Then strip trailing `_L` etc

        import pandas as pd

        df = pd.read_csv(input.metadata)
        # Subset to columns that are either not segment-specific or specific to the current segment
        segment_specific_cols = [
            col for col in df.columns if col.endswith(f"_{wildcards.segment}")
        ]
        cols_to_keep = [
            col
            for col in df.columns
            if not any(col.endswith(suffix) for suffix in ["_L", "_M", "_S"])
        ] + segment_specific_cols
        df = df[cols_to_keep]
        # Rename columns to strip segment suffixes
        new_column_names = {
            col: (
                col.rsplit("_", 1)[0]
                if col.endswith(f"_{wildcards.segment}")
                else col
            )
            for col in df.columns
        }
        df = df.rename(columns=new_column_names)
        df.to_csv(output.subset_metadata, index=False)


rule strip_segment_from_id:
    input:
        sequences="data/{species}/{segment}/sequences.fasta",
    output:
        sequences="data/{species}/{segment}/sequences_stripped.fasta",
    shell:
        """
        seqkit replace -p '\|.*' -r '' {input.sequences} > {output.sequences}
        """


rule format_ppx_ndjson:
    input:
        sequences="data/{species}/{segment}/sequences_stripped.fasta",
        metadata="data/{species}/{segment}/ppx_metadata.csv",
    output:
        ndjson="data/{species}/{segment}/ppx.ndjson",
    shell:
        """
        augur curate passthru \
            --metadata {input.metadata} \
            --fasta {input.sequences} \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            --seq-id-column accessionVersion \
            --seq-field sequence \
            > {output.ndjson}
        """
