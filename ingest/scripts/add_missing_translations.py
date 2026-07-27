
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="add a synthetic CDS to genbank records that have none"
    )
    parser.add_argument("--gb_file", help="Path to file containing genbank records")
    parser.add_argument("--output", help="Path to output file")

    return parser.parse_args()


def find_largest_orf(seq):
    """Return ``(start, end, protein)`` for the largest ORF across the three
    forward reading frames, where an ORF runs from its first Met to the next stop.

    Coordinates are 0-based, end-exclusive, and the stop codon is excluded, so they
    satisfy ``str(seq[start:end].translate()) == protein``. If a stop-delimited
    stretch contains no Met, the whole stretch is used. Returns ``(0, 0, "")`` for a
    sequence with no codons.
    """
    best = None  # (protein_length, start, end, protein)
    for frame in range(3):
        length = 3 * ((len(seq) - frame) // 3)  # Multiple of three
        translation = str(seq[frame:frame + length].translate())

        aa_start = 0
        while aa_start < len(translation):
            aa_end = translation.find("*", aa_start)
            if aa_end == -1:
                aa_end = len(translation)

            # Start the ORF at the first Met of this stop-delimited stretch. The
            # coordinate must be advanced by the same offset as the protein,
            # otherwise the location and the translation disagree.
            m_index = translation.find("M", aa_start, aa_end)
            orf_aa_start = m_index if m_index != -1 else aa_start
            protein = translation[orf_aa_start:aa_end]

            if best is None or len(protein) > best[0]:
                best = (
                    len(protein),
                    frame + orf_aa_start * 3,
                    frame + aa_end * 3,
                    protein,
                )

            aa_start = aa_end + 1

    if best is None:
        return 0, 0, ""
    _, start, end, protein = best
    return start, end, protein


def has_cds_translation(record):
    """True if the record already has a CDS feature carrying a translation."""
    return any(
        feature.type == "CDS" and feature.qualifiers.get("translation")
        for feature in record.features
    )


def main():
    args = parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(args.gb_file, "genbank"))
    added = 0
    for record in records:
        if has_cds_translation(record):
            continue
        start, end, protein = find_largest_orf(record.seq)
        record.features.append(
            SeqFeature(
                FeatureLocation(start, end),
                type="CDS",
                qualifiers={
                    "translation": [protein],
                    "codon_start": ["1"],
                    "note": [
                        "CDS added in-house: largest ORF (first Met to stop); "
                        "the source record carried no CDS/translation"
                    ],
                },
            )
        )
        added += 1

    SeqIO.write(records, args.output, "genbank")
    print(f"Added a synthetic largest-ORF CDS to {added} of {len(records)} records.")


if __name__ == "__main__":
    main()
