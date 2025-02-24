from Bio.Seq import Seq
from collections import Counter
import hashlib
import os
from utils.plot import plot_nucleotide_distribution, plot_amino_acid_frequencies

def compute_gc_content(sequence):
    """Compute the GC content percentage of a sequence."""
    sequence = sequence.upper()
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    return gc_content

def compute_nucleotide_percentage(sequence):
    """Compute the percentage of each nucleotide in the sequence."""
    length = len(sequence)
    return {nt: sequence.count(nt) / length * 100 for nt in 'ATGC'}

def count_amino_acids(protein_sequence):
    """Count the occurrences of each amino acid in a protein sequence."""
    return Counter(protein_sequence)

def extract_top_proteins(protein_sequence):
    """Extract the top 5 longest protein sequences and save them as FASTA files."""
    # Split the protein sequence by '*' to get individual proteins
    proteins = protein_sequence.split('*')
    
    # Filter out short sequences and sort by length in descending order
    long_proteins = sorted([protein for protein in proteins if len(protein) > 20], key=len, reverse=True)[:5]

    fasta_files = []
    for idx, protein in enumerate(long_proteins, start=1):
        # Generate a unique filename using a hash of the protein sequence
        protein_hash = hashlib.md5(protein.encode()).hexdigest()
        filename = os.path.join('top5proteins', f"protein_{idx}_{protein_hash}.fasta")
        with open(filename, "w") as f:
            f.write(f">Protein {idx}\n{protein}\n")
        fasta_files.append(filename)

    return fasta_files

def process_sequences(sequences):
    results = []

    for seq_record in sequences:
        seq_str = str(seq_record.seq)
        trimmed_sequence = seq_str[:len(seq_str) - (len(seq_str) % 3)]
        dna_seq = Seq(trimmed_sequence)
        rna_sequence = dna_seq.transcribe()
        protein_sequence = rna_sequence.translate()
        gc_content = compute_gc_content(trimmed_sequence)
        nucleotide_counts = compute_nucleotide_percentage(trimmed_sequence)
        amino_acid_counts = count_amino_acids(str(protein_sequence))
        total_amino_acids = sum(amino_acid_counts.values())

        # Extract top 5 longest protein sequences
        top_5_proteins = sorted(
            [protein for protein in str(protein_sequence).split('*') if len(protein) > 20],
            key=len, reverse=True
        )[:5]

        top_5_proteins_files = extract_top_proteins(str(protein_sequence))

        results.append({
            'id': seq_record.id,
            'genome_name': seq_record.description,
            'dna': seq_str,
            'rna': str(rna_sequence),
            'protein': str(protein_sequence),
            'gc_content': gc_content,
            'nucleotide_distribution': plot_nucleotide_distribution(nucleotide_counts),
            'amino_acid_counts': amino_acid_counts.most_common(20),
            'amino_acid_distribution': plot_amino_acid_frequencies(amino_acid_counts),
            'total_amino_acids': total_amino_acids,
            'top_long_proteins': top_5_proteins,
            'top_long_proteins_files': top_5_proteins_files
        })

    return results