from Bio import SeqIO
from pathlib import Path
import csv


def gc_content(seq):
    seq = seq.upper()
    g = seq.count("G")
    c = seq.count("C")
    return round((g + c) / len(seq) * 100, 2)


def load_sequence(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    return str(record.seq)


if __name__ == "__main__":

    fasta_folder = Path("data/raw/fasta files")
    output_file = Path("results/gc_content.csv")

    results = []

    for fasta_file in fasta_folder.glob("*.fasta"):

        sequence = load_sequence(fasta_file)
        gc = gc_content(sequence)

        results.append({
            "gene": fasta_file.stem,
            "sequence_length": len(sequence),
            "gc_content": gc
        })

        print(f"{fasta_file.name}:")
        print("  Sequence length:", len(sequence))
        print("  GC content:", gc, "%")

    with open(output_file, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=["gene", "sequence_length", "gc_content"]
        )

        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved results to {output_file}")