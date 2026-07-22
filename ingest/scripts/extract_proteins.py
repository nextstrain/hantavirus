
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import numpy as np
import pandas as pd
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="extract protein sequences out of genbank files"
    )
    parser.add_argument("--input", help="Path to input file")
    parser.add_argument("--special_cases", help="Path to CSV file with special cases")
    parser.add_argument("--output", help="Path to output fasta file")

    return parser.parse_args()


def main():
    args = parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    records = SeqIO.parse(args.input, "genbank")
    records_to_write = []

    # This extracts the AA sequence out of the genbank records
    for record in records:
        largest_record = None
        for feature in record.features:
            if feature.type == "CDS":
                translation = feature.qualifiers.get("translation")
                if translation:
                    product = feature.qualifiers.get("product", ["unknown"])[0]
                    accession = record.annotations.get(
                        "accessions", ["unknown"]
                        )[0]
                    source = record.annotations.get("source", "unknown")
                    source = source.replace(" ", "_")
                    protein_record = SeqRecord(
                            seq=Seq(translation[0]),
                            id=accession + '|' + source,
                            description=product,
                        )

                    if largest_record is None or len(protein_record.seq) > len(largest_record.seq):
                        largest_record = protein_record

        if largest_record is None:
            print("Error")
        else:
            records_to_write.append(largest_record)

    # This classifies the records as S, M or L segment
    special_cases = pd.read_csv(args.special_cases)
    special_cases_dict = dict(zip(special_cases["seq_id"], special_cases["classification"]))
    s_records = []
    m_records = []
    l_records = []

    for record in records_to_write:
        record_split = record.id.split("|")
        if record_split[0] in special_cases_dict.keys():
            classification = special_cases_dict[record_split[0]]

            if classification == "exclude":
                pass
            elif classification == "L":
                l_records.append(record)
            elif classification == "M":
                m_records.append(record)
            elif classification == "S":
                s_records.append(record)

        elif np.abs(len(record.seq)-450) < 100:
            s_records.append(record)
        elif np.abs(len(record.seq)-1200) < 100:
            m_records.append(record)
        elif np.abs(len(record.seq) - 2000) < 300:
            l_records.append(record)

    clean_records = []
    for segment_name, record_list in [("S", s_records), ("M", m_records), ("L", l_records)]:
        for record in record_list:
            record.id = record.id + "|" + segment_name
            record.description = ""
            clean_records.append(record)

    SeqIO.write(clean_records, args.output, "fasta")


if __name__ == "__main__":
    main()
