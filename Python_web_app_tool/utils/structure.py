from Bio import PDB
from Bio.PDB import *
import os
import logging
import requests

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class ProteinStructureHandler:
    def __init__(self, cache_dir):
        self.cache_dir = cache_dir
        self.pdb_list = PDB.PDBList()
        self.parser = PDB.PDBParser(QUIET=True)
        logger.info(f"Initialized ProteinStructureHandler with cache dir: {cache_dir}")
        
    def fetch_structure(self, pdb_id):
        """Fetch and parse a PDB structure."""
        try:
            logger.info(f"[DEBUG] Starting structure fetch for PDB ID: {pdb_id}")
            
            # Create PDB cache directory
            pdb_dir = os.path.join(self.cache_dir, 'pdb_structures')
            os.makedirs(pdb_dir, exist_ok=True)
            logger.info(f"[DEBUG] Using PDB directory: {pdb_dir}")
            
            # Check if file already exists in cache with correct name
            pdb_file = os.path.join(pdb_dir, f"{pdb_id.upper()}.pdb")
            
            if os.path.exists(pdb_file):
                logger.info(f"[DEBUG] Found cached PDB file: {pdb_file}")
            else:
                logger.info(f"[DEBUG] Downloading PDB file for {pdb_id}")
                # Download directly from RCSB
                url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
                response = requests.get(url)
                if response.status_code != 200:
                    raise Exception(f"Failed to download PDB file: {response.status_code}")
                    
                with open(pdb_file, 'wb') as f:
                    f.write(response.content)
                logger.info(f"[DEBUG] Downloaded and saved PDB file to: {pdb_file}")
            
            # Parse structure for metadata
            structure = self.parser.get_structure(pdb_id, pdb_file)
            chains = list(structure.get_chains())
            
            structure_info = {
                'pdb_id': pdb_id,
                'file_path': pdb_file,
                'chains': len(chains),
                'residues': sum(1 for _ in structure.get_residues())
            }
            
            logger.info(f"[DEBUG] Successfully processed structure: {structure_info}")
            return structure_info
            
        except Exception as e:
            logger.error(f"[ERROR] Failed to fetch structure {pdb_id}: {str(e)}", exc_info=True)
            raise
            
    def _get_resolution(self, structure):
        """Get structure resolution if available."""
        if hasattr(structure, 'header'):
            return structure.header.get('resolution')
        return None
        
    def _get_chain_info(self, chains):
        """Get information about each chain."""
        info = []
        for chain in chains:
            residues = list(chain.get_residues())
            if residues:
                info.append({
                    'chain_id': chain.id,
                    'length': len(residues),
                    'residue_range': f"{residues[0].id[1]}-{residues[-1].id[1]}"
                })
        return info