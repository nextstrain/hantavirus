
import pandas as pd
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add assigned segments and species to metadata"
    )

    parser.add_argument(
        "--diamond_tsv",
        help="Path to CSV file with assigned segments & species"
        )
    parser.add_argument("--metadata", help="Path to metadata file")
    parser.add_argument(
        "--rename_species", help="Dictionary of old species name -> new species name"
        )
    parser.add_argument("--output", help="Path to output file")

    return parser.parse_args()


def main():
    args = parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    diamond_tsv = pd.read_csv(args.diamond_tsv, sep='\t', header=None)
    diamond_tsv.columns = [
        "accession", "sseqid", "pident", "length", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore"
        ]
    split_sseqid = diamond_tsv["sseqid"].str.split('|', expand=True)
    diamond_tsv["assigned_species_id"] = split_sseqid[0]

    diamond_tsv["assigned_species"] = split_sseqid[1]
    species_renames_df = pd.read_csv(args.rename_species)
    species_renames = {
        row["old_name"]: row["new_name"] for idx, row in species_renames_df.iterrows()
        }
    diamond_tsv["assigned_species"] = diamond_tsv["assigned_species"].replace(species_renames)

    diamond_tsv["assigned_segment"] = split_sseqid[2]
    split_accession = diamond_tsv["accession"].str.split('.', expand=True)
    diamond_tsv["accession"] = split_accession[0]
    diamond_tsv["accession_version"] = split_accession[1]

    metadata = pd.read_csv(args.metadata, sep='\t')

    diamond_subset = diamond_tsv[
        ["accession", "assigned_species", "assigned_segment", "pident"]
        ]
    new_df = pd.merge(metadata, diamond_subset, on="accession")

    new_df.to_csv(args.output, sep='\t', index=False)


if __name__ == "__main__":
    main()
