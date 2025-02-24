from Bio import Entrez, SeqIO
import logging

# NCBI configuration
Entrez.email = "btc.cchrys@gmail.com"  # Replace with your actual email

def read_sequence_ids(file_path):
    """Read sequence IDs from a file."""
    with open(file_path, "r") as f:
        return [line.strip() for line in f if line.strip()]

def fetch_sequences_from_ncbi(sequence_ids):
    """Fetch sequences from NCBI given a list of sequence IDs."""
    sequences = []
    for seq_id in sequence_ids:
        try:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            sequences.append(SeqIO.read(handle, "fasta"))
            handle.close()
        except Exception as e:
            logging.error(f"Error fetching {seq_id}: {e}")
    return sequences