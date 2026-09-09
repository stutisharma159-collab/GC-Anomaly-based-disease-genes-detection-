from Bio import Entrez, SeqIO
from pathlib import Path
import time

Entrez.email = "stutisharma159@gmail.com"


OUTPUT_FOLDER = Path("data/raw/fasta files")
OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)


def fetch_gene_fasta(gene_name, organism="Homo sapiens", output_folder=OUTPUT_FOLDER):
    """
    Search NCBI for a human RefSeq mRNA for the given gene
    and save the sequence as a FASTA file.
    Returns the output path on success, None on failure.
    """
    gene_name = gene_name.strip().upper()
    if not gene_name:
        return None

    output_file = output_folder / f"{gene_name}.fasta"
    if output_file.exists():
        overwrite = input(f"{gene_name}.fasta already exists — overwrite? (y/n): ").strip().lower()
        if overwrite != "y":
            print(f"Skipped {gene_name}")
            return None

    print(f"\nSearching NCBI for {gene_name}...")
    search_term = (
        f"{gene_name}[Gene Name] "
        f"AND {organism}[Organism] "
        f"AND RefSeq[Filter] "
        f"AND biomol_mrna[PROP]"
    )

    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        print(f"No sequence found for {gene_name} — check spelling, or it may not be a valid gene symbol")
        return None

    sequence_id = record["IdList"][0]
    handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
    sequence_record = SeqIO.read(handle, "fasta")
    handle.close()

    sequence_record.id = gene_name
    sequence_record.name = gene_name
    sequence_record.description = f"{gene_name} downloaded from NCBI"

    SeqIO.write(sequence_record, output_file, "fasta")
    print(f"Saved {gene_name}: {len(sequence_record.seq)} bp -> {output_file}")

    time.sleep(0.4)  # stay polite to NCBI servers
    return output_file


def fetch_multiple_genes(gene_list, organism="Homo sapiens"):
    results = {"success": [], "failed": []}
    for gene in gene_list:
        path = fetch_gene_fasta(gene, organism)
        if path:
            results["success"].append(gene)
        else:
            results["failed"].append(gene)
    return results


if __name__ == "__main__":
    raw_input = input("Enter gene names, separated by commas (e.g. ACTB, FMR1, HTT): ")
    gene_list = [g.strip() for g in raw_input.split(",") if g.strip()]

    print(f"\nFetching {len(gene_list)} genes: {gene_list}")
    results = fetch_multiple_genes(gene_list)

    print(f"\nDone. Success: {len(results['success'])}, Failed: {len(results['failed'])}")
    if results["failed"]:
        print(f"Failed genes (check spelling/symbol): {results['failed']}")