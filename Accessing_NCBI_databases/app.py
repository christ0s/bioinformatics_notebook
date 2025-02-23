from flask import Flask, render_template, request, jsonify
from markupsafe import Markup
import time
import os
import io
import base64
import logging
import pickle
from collections import Counter
from Bio import Entrez, SeqIO, SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
import matplotlib.pyplot as plt
import requests

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['TOP_PROTEINS_DIR'] = 'top5proteins'
app.config['CACHE_FOLDER'] = 'cache'

# Ensure necessary directories exist
os.makedirs(app.config['CACHE_FOLDER'], exist_ok=True)
os.makedirs(app.config['TOP_PROTEINS_DIR'], exist_ok=True)

# Set your email (NCBI requires this for access)
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
            print(f"Fetching {seq_id} from NCBI...")
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            genome_name = record.description 
            sequences.append(record)
            handle.close()
        except Exception as e:
            print(f"Error fetching {seq_id}: {e}")
    return sequences

def compute_gc_content(sequence):
    """Compute the GC content percentage of a sequence."""
    g_count = sequence.count('G') + sequence.count('g')
    c_count = sequence.count('C') + sequence.count('c')
    gc_content = (g_count + c_count) / len(sequence) * 100
    return gc_content

def compute_nucleotide_percentage(sequence):
    """Compute the percentage of each nucleotide in the sequence."""
    length = len(sequence)
    return {
        'A': sequence.count('A') / length * 100,
        'T': sequence.count('T') / length * 100,
        'G': sequence.count('G') / length * 100,
        'C': sequence.count('C') / length * 100
    }

def count_amino_acids(protein_sequence):
    """Count the occurrences of each amino acid in a protein sequence."""
    return Counter(protein_sequence)

def plot_nucleotide_distribution(nucleotide_counts, title="Nucleotide Composition"):
    """Create a bar chart and return it as an HTML-safe image."""
    plt.figure(figsize=(4, 4))
    plt.bar(nucleotide_counts.keys(), nucleotide_counts.values())
    plt.ylabel("Percentage")
    plt.title(title)
    # Convert plot to PNG image
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close()
    return Markup(f'<img src="data:image/png;base64,{image_base64}" />')
def plot_amino_acid_frequencies(amino_acid_counts):
    """Generate a bar plot for the top 20 amino acid frequencies and return it as an HTML-safe image."""
    # Extract amino acids and counts, sort by count in descending order
    sorted_aa_counts = sorted(amino_acid_counts.items(), key=lambda x: x[1], reverse=True)[:20]
    amino_acids = [aa for aa, count in sorted_aa_counts]
    counts = [count for aa, count in sorted_aa_counts]
    
    # Create the bar plot
    plt.figure(figsize=(8, 4))
    plt.bar(amino_acids, counts)
    plt.xlabel('Amino Acids')
    plt.ylabel('Frequency')
    plt.title('Top 20 Amino Acid Frequencies')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    # Convert plot to PNG image
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close()
    return Markup(f'<img src="data:image/png;base64,{image_base64}" />')
@app.route("/", methods=["GET", "POST"])
def home():
    if request.method == "POST":
        # Check if a file was submitted
        if "file" not in request.files:
            return "No file part"
        file = request.files["file"]
        if file.filename == "":
            return "No selected file"
        if file:
            # Ensure upload folder exists
            os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)

            # Process the uploaded file
            sequence_ids = read_sequence_ids(filepath)
            sequences = fetch_sequences_from_ncbi(sequence_ids)

            results = []
            for seq_record in sequences:
                seq_str = str(seq_record.seq)
                # Adjust sequence length to be a multiple of three
                trimmed_sequence = seq_str[:len(seq_str) - (len(seq_str) % 3)]
                dna_seq = Seq(trimmed_sequence)
                rna_sequence = dna_seq.transcribe()
                protein_sequence = rna_sequence.translate()
                gc_content = compute_gc_content(trimmed_sequence)
                nucleotide_counts = compute_nucleotide_percentage(trimmed_sequence)
                amino_acid_counts = count_amino_acids(str(protein_sequence))
                total_amino_acids = sum(amino_acid_counts.values())
                protein_fragments = (protein_sequence).split('*')
                long_proteins = [protein for protein in protein_fragments if len(protein) > 20]

                long_proteins.sort(key=len, reverse=True)
                top_5_proteins = long_proteins[:5]


                # Generate nucleotide distribution plot
                nucleotide_distribution_plot = plot_nucleotide_distribution(
                    nucleotide_counts, title=f"Nucleotide Composition for {seq_record.id}"
                )
                # Generate amino acid frequencies plot
                amino_acid_frequencies_plot = plot_amino_acid_frequencies(amino_acid_counts)

                results.append({
                    'id': seq_record.id,
                    'genome_name': seq_record.description,
                    'dna': seq_str,
                    'rna': str(rna_sequence),
                    'protein': str(protein_sequence),
                    'gc_content': gc_content,
                    'amino_acid_counts': amino_acid_counts.most_common(20),
                    'nucleotide_distribution_plot': nucleotide_distribution_plot,
                    'amino_acid_frequencies_plot': amino_acid_frequencies_plot,
                    'total_amino_acids': total_amino_acids,
                    'top_5_proteins': top_5_proteins
                })

            return render_template("results.html", results=results)

    return render_template("index.html")

if __name__ == "__main__":
    app.run(debug=True)
