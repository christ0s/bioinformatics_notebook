import os
from flask import Flask, render_template, request, jsonify, send_file
import requests
from flask_caching import Cache
import subprocess
import tempfile
from utils.fetch import read_sequence_ids, fetch_sequences_from_ncbi
from utils.process import process_sequences, extract_top_proteins
from utils.blast import run_blast_search  # Use Biopython qblast here
from utils.structure import ProteinStructureHandler

app = Flask(__name__)
app.config.update(
    UPLOAD_FOLDER='uploads',
    TOP_PROTEINS_DIR='top5proteins',
    CACHE_FOLDER='cache',
    CACHE_TYPE='filesystem',             # Use file‑based caching
    CACHE_DEFAULT_TIMEOUT=3600,
    CACHE_DIR=os.path.join(os.getcwd(), 'cache'),  # Directory for cache files
    CACHE_THRESHOLD=1000
)
cache = Cache(app)

# Add cache directory to store PDB files
PDB_CACHE_DIR = os.path.join(app.config['CACHE_DIR'], 'pdb_structures')
os.makedirs(PDB_CACHE_DIR, exist_ok=True)

# Initialize structure handler (it handles pdb download)
structure_handler = ProteinStructureHandler(app.config['CACHE_FOLDER'])

# Ensure necessary directories exist
for folder in [app.config['CACHE_FOLDER'], app.config['TOP_PROTEINS_DIR'], app.config['UPLOAD_FOLDER']]:
    os.makedirs(folder, exist_ok=True)

def perform_blast_search_qblast(protein_sequence):
    """Write the protein sequence to a temporary FASTA file and perform BLAST search."""
    try:
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as tmp_file:
            tmp_file.write(f">query\n{protein_sequence}")
            tmp_file_path = tmp_file.name
        results = run_blast_search(tmp_file_path)
        return results
    except Exception as e:
        raise Exception(f"BLAST search failed: {e}")
    finally:
        if os.path.exists(tmp_file_path):
            os.remove(tmp_file_path)

@app.route("/blast_search", methods=["POST"])
def blast_search():
    data = request.get_json()
    protein_sequence = data.get('protein_sequence')
    if not protein_sequence:
        return jsonify({'error': 'No protein sequence provided'}), 400
    cached_result = cache.get(protein_sequence)
    if cached_result:
        return jsonify(cached_result)
    try:
        result = perform_blast_search_qblast(protein_sequence)
    except Exception as e:
        return jsonify({'error': str(e)}), 500
    cache.set(protein_sequence, result)
    return jsonify(result)

@app.route("/get_protein_structure/<pdb_id>", methods=["GET"])
def get_protein_structure(pdb_id):
    print(f"[DEBUG] Fetching structure for PDB ID: {pdb_id}")
    try:
        pdb_file = os.path.join(PDB_CACHE_DIR, f"{pdb_id.upper()}.pdb")
        
        # Download PDB file if not in cache
        if not os.path.exists(pdb_file):
            print(f"[DEBUG] Downloading PDB file from RCSB")
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            response = requests.get(url)
            
            if response.status_code != 200:
                raise Exception(f"Failed to download PDB file: {response.status_code}")
            
            # Save PDB file
            with open(pdb_file, 'wb') as f:
                f.write(response.content)
            print(f"[DEBUG] PDB file saved to {pdb_file}")
        
        structure_info = {
            'pdb_id': pdb_id,
            'file_path': pdb_file
        }
        return jsonify(structure_info)
        
    except Exception as e:
        print(f"[ERROR] {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route("/pdb/<pdb_id>", methods=["GET"])
def serve_pdb_file(pdb_id):
    print(f"[DEBUG] Serving PDB file for ID: {pdb_id}")
    try:
        pdb_file = os.path.join(PDB_CACHE_DIR, f"{pdb_id.upper()}.pdb")
        if not os.path.exists(pdb_file):
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            response = requests.get(url)
            if response.status_code != 200:
                raise Exception(f"Failed to download PDB file: {response.status_code}")
            with open(pdb_file, 'wb') as f:
                f.write(response.content)
            print(f"[DEBUG] PDB file saved successfully")
        return send_file(
            pdb_file,
            mimetype='chemical/x-pdb',
            as_attachment=True,
            download_name=f"{pdb_id.upper()}.pdb"
        )
    except Exception as e:
        print(f"[ERROR] {str(e)}")
        return jsonify({'error': str(e)}), 500

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
