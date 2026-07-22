
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="add genbank protein records"
    )
    parser.add_argument("--gb_file", help="Path to file containing genbank records")
    parser.add_argument("--output", help="Path to output file")

    return parser.parse_args()


def main():
    args = parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    records = SeqIO.parse(args.gb_file, "genbank")
    records_to_write = []
    missing_records = []
    for record in records:
        found_translation = False

        for feature in record.features:
            if feature.type == "CDS":
                translation = feature.qualifiers.get("translation")
                if translation:
                    found_translation = True
                    break

        if found_translation:
            records_to_write.append(record)
        else:
            missing_records.append(record.id)

            source = record.annotations.get("source", "unknown")
            source = source.replace(" ", "_")

            longest_translation = None
            for frame in range(3):
                length = 3 * ((len(record.seq) - frame) // 3)  # Multiple of three

                aa_start = 0
                translation = str(record.seq[frame: frame + length].translate())

                while aa_start < len(translation):
                    aa_end = translation.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = len(translation)

                    current_protein = translation[aa_start:aa_end]
                    m_index = current_protein.find("M")
                    if m_index != -1:
                        current_protein = current_protein[m_index:]

                    if (
                        longest_translation is None
                        or len(current_protein) > len(longest_translation)
                    ):
                        longest_translation = current_protein
                        nt_start = frame + aa_start*3
                        nt_end = frame + aa_end*3

                    aa_start = aa_end + 1

                feature = SeqFeature(
                    FeatureLocation(nt_start, nt_end),
                    type="CDS",
                    qualifiers={"translation": [str(longest_translation)]}
                )

            record.features.append(feature)
            records_to_write.append(record)

    SeqIO.write(records_to_write, args.output, "genbank")


if __name__ == "__main__":
    main()
