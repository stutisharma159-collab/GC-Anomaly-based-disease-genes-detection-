from Bio import SeqIO
from pathlib import Path
import csv
import statistics


def gc_content(seq):
    seq = seq.upper()
    g = seq.count("G")
    c = seq.count("C")
    return (g + c) / len(seq) * 100


def load_sequence(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    return str(record.seq)


def sliding_window_gc(sequence, window_size=500, step_size=100):

    results = []

    for start in range(0, len(sequence) - window_size + 1, step_size):

        window = sequence[start:start + window_size]

        gc = gc_content(window)

        results.append({
            "start": start + 1,
            "end": start + window_size,
            "gc_content": round(gc, 2)
        })

    return results


def calculate_spikiness(window_results):

    gc_values = [row["gc_content"] for row in window_results]

    mean_gc = statistics.mean(gc_values)
    std_gc = statistics.stdev(gc_values)

    for row in window_results:

        if std_gc == 0:
            row["z_score"] = 0
        else:
            row["z_score"] = round(
                (row["gc_content"] - mean_gc) / std_gc,
                2
            )

    return window_results, mean_gc, std_gc


if __name__ == "__main__":

    fasta_folder = Path("data/raw/fasta files")
    output_folder = Path("results/sliding_windows")

    output_folder.mkdir(parents=True, exist_ok=True)

    for fasta_file in fasta_folder.glob("*.fasta"):

        sequence = load_sequence(fasta_file)

        window_results = sliding_window_gc(
            sequence,
            window_size=500,
            step_size=100
        )

        window_results, mean_gc, std_gc = calculate_spikiness(
            window_results
        )

        output_file = (
            output_folder /
            f"{fasta_file.stem}_sliding_gc.csv"
        )

        with open(output_file, "w", newline="") as f:

            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "start",
                    "end",
                    "gc_content",
                    "z_score"
                ]
            )

            writer.writeheader()
            writer.writerows(window_results)

        max_spike = max(
            window_results,
            key=lambda row: row["z_score"]
        )

        print(f"\n{fasta_file.name}")
        print(f"  Mean window GC: {mean_gc:.2f}%")
        print(f"  Window GC SD: {std_gc:.2f}")
        print(
            f"  Maximum GC spike: "
            f"{max_spike['gc_content']:.2f}% "
            f"(z = {max_spike['z_score']:.2f})"
        )
        print(f"  Saved: {output_file}")