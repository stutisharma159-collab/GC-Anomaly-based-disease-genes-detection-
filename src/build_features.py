import os
import csv

from gc_content import load_sequence
from sliding_window import sliding_gc


RAW_DIR = "data/raw/fasta files"
OUTPUT_FILE = "results/features.csv"


def compute_features(gene_name, seq):

    profile = sliding_gc(seq)

    mean_gc = sum(profile) / len(profile)
    max_gc = max(profile)
    min_gc = min(profile)

    std_gc = (
        sum((x - mean_gc) ** 2 for x in profile) / len(profile)
    ) ** 0.5

    # How much the highest GC window rises above
    # the overall average GC of the sliding windows
    spikiness = max_gc - mean_gc

    return {
        "gene": gene_name,
        "mean_gc": round(mean_gc, 2),
        "max_gc": round(max_gc, 2),
        "min_gc": round(min_gc, 2),
        "std_gc": round(std_gc, 2),
        "spikiness": round(spikiness, 2),
    }


if __name__ == "__main__":

    rows = []

    for filename in os.listdir(RAW_DIR):

        if filename.endswith(".fasta"):

            gene_name = filename.replace(".fasta", "")

            fasta_path = os.path.join(RAW_DIR, filename)

            seq = load_sequence(fasta_path)

            features = compute_features(gene_name, seq)

            rows.append(features)

    if len(rows) == 0:
        print("No FASTA files found.")
    else:

        with open(OUTPUT_FILE, "w", newline="") as f:

            writer = csv.DictWriter(
                f,
                fieldnames=rows[0].keys()
            )

            writer.writeheader()
            writer.writerows(rows)

        print(f"Saved features for {len(rows)} genes to {OUTPUT_FILE}")