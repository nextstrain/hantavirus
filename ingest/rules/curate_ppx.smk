"""
This part of the workflow handles the curation of data from PPX

"""


def format_field_map(field_map: dict[str, str]) -> str:
    """
    Format dict to `"key1"="value1" "key2"="value2"...` for use in shell commands.
    """
    return " ".join([f'"{key}"="{value}"' for key, value in field_map.items()])


# This curate pipeline is based on existing pipelines for pathogen repos using NCBI data.
# You may want to add and/or remove steps from the pipeline for custom metadata
# curation for your pathogen. Note that the curate pipeline is streaming NDJSON
# records between scripts, so any custom scripts added to the pipeline should expect
# the input as NDJSON records from stdin and output NDJSON records to stdout.
# The final step of the pipeline should convert the NDJSON records to two
# separate files: a metadata TSV and a sequences FASTA.
rule curate_ppx:
    input:
        sequences_ndjson="data/{species}/{segment}/ppx.ndjson",
        geolocation_rules=config["curate"]["local_geolocation_rules"]
    output:
        metadata="results/{species}/{segment}/metadata.tsv",
        sequences="results/{species}/{segment}/sequences.fasta",
    params:
        field_map=format_field_map(config["curate"]["ppx_field_map"]),
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
        abbr_authors_field=config["curate"]["abbr_authors_field"],
        id_field=config["curate"]["output_id_field"],
        sequence_field=config["curate"]["output_sequence_field"],
    wildcard_constraints:
        species="andv",
    shell:
        r"""
        cat {input.sequences_ndjson:q} \
            | augur curate rename \
                --field-map {params.field_map} \
            | augur curate normalize-strings \
            | augur curate transform-strain-name \
                --strain-regex {params.strain_regex:q} \
                --backup-fields {params.strain_backup_fields:q} \
            | augur curate format-dates \
                --date-fields {params.date_fields:q} \
                --expected-date-formats {params.expected_date_formats:q} \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields:q} \
                --articles {params.articles:q} \
                --abbreviations {params.abbreviations:q} \
            | augur curate abbreviate-authors \
                --authors-field {params.authors_field:q} \
                --default-value {params.authors_default_value:q} \
                --abbr-authors-field {params.abbr_authors_field:q} \
            | python {workflow.basedir}/scripts/curate-urls.py \
            | augur curate apply-geolocation-rules \
                --geolocation-rules {input.geolocation_rules:q} \
                --output-metadata {output.metadata:q} \
                --output-fasta {output.sequences:q} \
                --output-id-field {params.id_field:q} \
                --output-seq-field {params.sequence_field:q}
        """


rule extract_ppx_data:
    input:
        metadata="results/{species}/{segment}/metadata.tsv",
        sequences="results/{species}/{segment}/sequences.fasta",
    output:
        metadata="results/{species}/{segment}/metadata_{ppx_dut}.tsv",
        sequences="results/{species}/{segment}/sequences_{ppx_dut}.fasta",
    wildcard_constraints:
        ppx_dut="open|restricted",
    params:
        ppx_dut=lambda w: w.ppx_dut.upper(),
        # Warn on empty output for restricted since it's feasible that
        # none of the data is restricted
        empty_output=lambda w: "warn" if w.ppx_dut == "restricted" else "error",
    shell:
        """
        augur filter --metadata {input.metadata} \
                     --sequences {input.sequences} \
                     --metadata-id-columns accession \
                     --exclude-all \
                     --include-where "dataUseTerms={params.ppx_dut:q}" \
                     --output-metadata {output.metadata} \
                     --output-sequences {output.sequences} \
                     --empty-output-reporting {params.empty_output:q}
        """
