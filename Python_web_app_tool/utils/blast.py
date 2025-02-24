from Bio import Entrez,SeqIO, SearchIO
from Bio.Blast import NCBIWWW
import os
import pickle
Entrez.email = "your_email@example.com"  # Optional by NCBI
Entrez.api_key = "ceaa76cc578dc1673192930863df8bd6cc08"  # Set API key globally

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

        cache_file = os.path.join('cache', f"{os.path.basename(fasta_file)}.pkl")  # Define cache_file
        with open(cache_file, "wb") as f:
            pickle.dump(results, f)

        print(f"BLAST search completed successfully. Found {len(results)} hits.")
        return results

    except Exception as e:
        print(f"Error during BLAST search: {e}")
        return []