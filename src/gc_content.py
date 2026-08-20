from Bio import SeqIO


def gc_content(seq):
    seq = seq.upper()
    g = seq.count("G")
    c = seq.count("C")
    return round((g + c) / len(seq) * 100, 2)


def load_sequence(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    return str(record.seq)