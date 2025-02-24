from flask import Flask, render_template, request, jsonify
from flask_caching import Cache
import os
import subprocess
import tempfile
from utils.fetch import read_sequence_ids, fetch_sequences_from_ncbi
from utils.process import process_sequences, extract_top_proteins
from utils.blast import run_blast_search  # Use Biopython qblast here

# Flask application setup
app = Flask(__name__)
app.config.update(
    UPLOAD_FOLDER='uploads',
    TOP_PROTEINS_DIR='top5proteins',
    CACHE_FOLDER='cache',
    CACHE_TYPE='simple',
    CACHE_DEFAULT_TIMEOUT=3600  # Cache for 1 hour
)
cache = Cache(app)

# Ensure necessary directories exist
for folder in [app.config['CACHE_FOLDER'], app.config['TOP_PROTEINS_DIR'], app.config['UPLOAD_FOLDER']]:
    os.makedirs(folder, exist_ok=True)

def perform_blast_search_qblast(protein_sequence):
    """
    Write the protein sequence to a temporary FASTA file and perform BLAST search
    using Biopython's NCBIWWW.qblast via the run_blast_search function.
    """
    try:
        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as tmp_file:
            tmp_file.write(f">query\n{protein_sequence}")
            tmp_file_path = tmp_file.name

        # Run the BLAST search using our utils/blast.py function.
        results = run_blast_search(tmp_file_path)
        return results
    except Exception as e:
        raise Exception(f"BLAST search failed: {e}")
    finally:
        # Clean up the temporary file.
        if os.path.exists(tmp_file_path):
            os.remove(tmp_file_path)

@app.route("/blast_search", methods=["POST"])
def blast_search():
    data = request.get_json()
    protein_sequence = data.get('protein_sequence')
    if not protein_sequence:
        return jsonify({'error': 'No protein sequence provided'}), 400

    # Use protein sequence as the cache key
    cached_result = cache.get(protein_sequence)
    if cached_result:
        return jsonify(cached_result)

    try:
        # Use Biopython qblast to search (no local blast installation needed)
        result = perform_blast_search_qblast(protein_sequence)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

    cache.set(protein_sequence, result)
    return jsonify(result)

@app.route("/", methods=["GET", "POST"])
def home():
    if request.method == "POST":
        file = request.files.get("file")
        if not file:
            return "No file part", 400
        if file.filename == "":
            return "No selected file", 400
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
        file.save(filepath)

        sequence_ids = read_sequence_ids(filepath)
        sequences = fetch_sequences_from_ncbi(sequence_ids)
        results = process_sequences(sequences)
        return render_template("results.html", results=results)

    return render_template("index.html")

if __name__ == "__main__":
    app.run(debug=True)
