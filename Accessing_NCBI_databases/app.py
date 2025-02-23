from flask import Flask, render_template, request, jsonify
import os
import hashlib
from utils.fetch import read_sequence_ids, fetch_sequences_from_ncbi
from utils.process import process_sequences, extract_top_proteins
from utils.blast import run_blast_search

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
        results = process_sequences(sequences)

        print("File processing completed successfully.")
        return render_template("results.html", results=results)

    return render_template("index.html")

if __name__ == "__main__":
    app.run(debug=True)
