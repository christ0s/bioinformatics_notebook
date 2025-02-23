from flask import Flask, render_template, request, jsonify
from markupsafe import Markup
import os
import io
import base64
import logging
import sys
import pickle
from collections import Counter
from Bio import Entrez, SeqIO, SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
import matplotlib.pyplot as plt
import hashlib

# Flask application setup
app = Flask(__name__)
app.config.update(
    UPLOAD_FOLDER='uploads',
    TOP_PROTEINS_DIR='top5proteins',
    CACHE_FOLDER='cache'
)

# Ensure necessary directories exist
for folder in [app.config['CACHE_FOLDER'], app.config['TOP_PROTEINS_DIR'], app.config['UPLOAD_FOLDER']]:
    os.makedirs(folder, exist_ok=True)

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
            print(f"Error fetching {seq_id}: {e}")
    return sequences

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

def plot_distribution(counts, title, x_label='Nucleotide', y_label='Percentage'):
    """Create a bar chart and return it as an HTML-safe image."""
    plt.figure(figsize=(5, 5))
    plt.bar(counts.keys(), counts.values())
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close()
    return Markup(f'<img src="data:image/png;base64,{image_base64}" />')

plot_nucleotide_distribution = lambda counts: plot_distribution(counts, "Nucleotide Composition")
plot_amino_acid_frequencies = lambda counts: plot_distribution(
    dict(sorted(counts.items(), key=lambda x: x[1], reverse=True)[:20]),
    "Top 20 Amino Acid Frequencies", x_label='Amino Acids', y_label='Frequency')

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
        filename = os.path.join(app.config['TOP_PROTEINS_DIR'], f"protein_{idx}_{protein_hash}.fasta")
        with open(filename, "w") as f:
            f.write(f">Protein {idx}\n{protein}\n")
        fasta_files.append(filename)

    return fasta_files

def run_blast_search(fasta_file):
    try:
        
        protein_seq = SeqIO.read(fasta_file, "fasta")
        print(f"BLAST search start with protein seq  :", {protein_seq.seq})
        result_handle = NCBIWWW.qblast("blastp", "pdb", protein_seq.seq)
        
        # Ensure to parse results correctly.
        blast_records = list(SearchIO.parse(result_handle, 'blast-xml'))  # Use parse instead of read

        results = []
        for blast_record in blast_records:
            for hit in blast_record.hits:
                for hsp in hit.hsps:
                    results.append({
                        "Sequence ID": hit.id,
                        "Description": hit.description,
                        "E-value": hsp.evalue,
                        "Bit Score": hsp.bitscore,
                        "Alignment": str(hsp.aln)
                    })

        cache_file = os.path.join(app.config['CACHE_FOLDER'], f"{os.path.basename(fasta_file)}.pkl")  # Define cache_file
        with open(cache_file, "wb") as f:
            pickle.dump(results, f)

        print(f"BLAST search completed successfully. Found {len(results)} hits.")
        return results

    except Exception as e:
        print(f"Error during BLAST search: {e}")
        return []
    
@app.route("/blast_search", methods=["POST"])
def blast_search():
    print("BLAST search endpoint was called.")
    data = request.get_json()

    if not data or "protein_sequence" not in data:
        print("No protein sequence provided for BLAST search.")
        return jsonify({"error": "No protein sequence provided"}), 400

    protein_sequence = data.get("protein_sequence")
    print(f"Processing BLAST search for sequence: {protein_sequence[:50]}...")

    # Generate a unique filename using a hash of the protein sequence
    protein_hash = hashlib.md5(protein_sequence.encode()).hexdigest()
    fasta_file = os.path.join(app.config['TOP_PROTEINS_DIR'], f"protein_{protein_hash}.fasta")
    
    # Save the protein sequence to the file
    with open(fasta_file, "w") as f:
        f.write(f">Protein\n{protein_sequence}\n")

    print(f"Protein sequence saved to file: {fasta_file}")

    blast_results = run_blast_search(fasta_file)  # Pass the correct file path
    print(f"BLAST search returned {len(blast_results)} results.")

    return jsonify(blast_results)


@app.route("/", methods=["GET", "POST"])
def home():
    if request.method == "POST":
        print("Processing uploaded file.")
        file = request.files.get("file")

        if not file:
            print("No file part in the request.")
            return "No file part", 400

        if file.filename == "":
            print("No selected file.")
            return "No selected file", 400

        filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
        file.save(filepath)
        print(f"File saved at: {filepath}")

        sequence_ids = read_sequence_ids(filepath)
        sequences = fetch_sequences_from_ncbi(sequence_ids)
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

        print("File processing completed successfully.")
        return render_template("results.html", results=results)

    return render_template("index.html")

if __name__ == "__main__":
    app.run(debug=True)
