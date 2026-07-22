
rule filter_by_pident:
    input:
        metadata="data/ncbi_dataset_seq_assigned.tsv",
        sequences="data/ncbi_dataset_sequences.fasta"
    params:
        species="{species}",
        segment="{segment}",
        pident=PIDENT,
    output:
        metadata = "results/{species}/{segment}/metadata_filtered.csv",
        sequences = "results/{species}/{segment}/sequences_filtered.fasta",
    shell:
        """
        python scripts/filter_by_pident.py \
        --metadata {input.metadata} \
        --sequences {input.sequences} \
        --species {params.species} \
        --segment {params.segment} \
        --pident {params.pident} \
        --output_metadata {output.metadata} \
        --output_sequences {output.sequences}
        """

rule format_ncbi_datasets_ndjson:
    input:
        sequences="results/{species}/{segment}/sequences_filtered.fasta",
        metadata="results/{species}/{segment}/metadata_filtered.csv",
    output:
        ndjson="results/{species}/{segment}/ncbi.ndjson",
    shell:
        """
        augur curate passthru \
            --metadata {input.metadata} \
            --fasta {input.sequences} \
            --seq-id-column accession_version \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            > {output.ndjson}
        """

def format_field_map(field_map: dict[str, str]) -> list[str]:
    """
    Format entries to the format expected by `augur curate --field-map`.

    When used in a Snakemake shell block, the list is automatically expanded and
    spaces are handled by quoted interpolation.
    """
    return  [f'{key}={value}' for key, value in field_map.items()]

rule curate:
    input:
        sequences_ndjson="results/{species}/{segment}/ncbi.ndjson",
        local_geolocation_rules=config["curate"]["local_geolocation_rules"],
    output:
        metadata="results/{species}/{segment}/metadata_curated.tsv",
        sequences="results/{species}/{segment}/sequences_curated.fasta",
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        strain_regex=config["curate"]["strain_regex"],
        strain_backup_fields=config["curate"]["strain_backup_fields"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        genbank_location_field=config["curate"]["genbank_location_field"],
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],
        authors_field=config["curate"]["authors_field"],
        authors_default_value=config["curate"]["authors_default_value"],
        # abbr_authors_field=config["curate"]["abbr_authors_field"],
        # annotations_id=config["curate"]["annotations_id"],
        id_field=config["curate"]["output_id_field"],
        sequence_field=config["curate"]["output_sequence_field"],
    shell:
        """
        cat {input.sequences_ndjson} \
            | augur curate rename \
                --field-map {params.field_map:q} \
            | augur curate normalize-strings \
            | augur curate transform-strain-name \
                --strain-regex {params.strain_regex} \
                --backup-fields {params.strain_backup_fields} \
            | augur curate format-dates \
                --date-fields {params.date_fields} \
                --expected-date-formats {params.expected_date_formats} \
            | augur curate parse-genbank-location \
                --location-field {params.genbank_location_field} \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields} \
                --articles {params.articles} \
                --abbreviations {params.abbreviations} \
            | augur curate abbreviate-authors \
                --authors-field {params.authors_field} \
                --default-value {params.authors_default_value} \
            | augur curate apply-geolocation-rules \
                --geolocation-rules {input.local_geolocation_rules} \
                --output-metadata {output.metadata} \
                --output-fasta {output.sequences} \
                --output-id-field {params.id_field} \
                --output-seq-field {params.sequence_field}
        """