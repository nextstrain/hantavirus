
import argparse
import pandas as pd
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="extract protein sequences out of genbank files"
    )
    parser.add_argument("--metadata", help="Path to input metadata file")
    parser.add_argument("--sequences", help="path to input sequences file")
    parser.add_argument("--species", help="Species")
    parser.add_argument("--segment", help="Segment")
    parser.add_argument("--pident", type=float, help="cutoff pident value")
    parser.add_argument("--output_metadata", help="Path to output metadata file")
    parser.add_argument("--output_sequences", help="Path to output sequences file")

    return parser.parse_args()


def main():
    args = parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')

    new_df = metadata[
        (metadata["assigned_species"] == args.species)
        & (metadata["pident"] > args.pident)
        & (metadata["assigned_segment"] == args.segment.upper())
    ]
    new_df.to_csv(args.output_metadata, sep='\t', index=False)

    wanted_ids = set(new_df["accession_version"])

    records_to_write = []
    for record in SeqIO.parse(args.sequences, "fasta"):
        if record.id in wanted_ids:
            records_to_write.append(record)
    SeqIO.write(records_to_write, args.output_sequences, "fasta")


if __name__ == "__main__":
    main()
