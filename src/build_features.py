import os
import csv

from gc_content import load_sequence
from sliding_window import sliding_window_gc, calculate_spikiness


RAW_DIR = "data/raw/fasta files"
OUTPUT_FILE = "results/features.csv"


def compute_features(gene_name, seq):

    window_results = sliding_window_gc(seq)
    window_results, mean_gc, std_gc = calculate_spikiness(window_results)

    profile = [row["gc_content"] for row in window_results]
    z_scores = [row["z_score"] for row in window_results]

    max_gc = max(profile)
    min_gc = min(profile)

    spikiness = max_gc - mean_gc
    max_z_score = max(z_scores)

    return {
        "gene": gene_name,
        "mean_gc": round(mean_gc, 2),
        "max_gc": round(max_gc, 2),
        "min_gc": round(min_gc, 2),
        "std_gc": round(std_gc, 2),
        "spikiness": round(spikiness, 2),
        "max_z_score": round(max_z_score, 2),
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