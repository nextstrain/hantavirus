
from Bio import Entrez, SeqIO
import argparse
import pandas as pd
from pathlib import Path

Entrez.email = "hello@nextstrain.org"


def parse_args():
    parser = argparse.ArgumentParser(
        description="download genbank records for list of accession ids"
    )
    parser.add_argument("--metadata", help="tsv file containing accession ids")
    parser.add_argument("--output", help="Path to output file")

    return parser.parse_args()


def main():
    args = parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.metadata, sep='\t')
    refseq_df = df[df["sourcedb"] == "RefSeq"]
    accessions = list(refseq_df["accession"])

    records = []
    for accession in accessions:
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="gb", retmode="text"
            )
        record = SeqIO.read(handle, "genbank")
        records.append(record)
    handle.close()

    SeqIO.write(records, args.output, "genbank")


if __name__ == "__main__":
    main()
