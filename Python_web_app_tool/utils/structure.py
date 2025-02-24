from Bio import PDB
from Bio.PDB import *
import os
import logging

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
            logger.info(f"Starting structure fetch for PDB ID: {pdb_id}")
            
            # Create PDB cache directory
            pdb_dir = os.path.join(self.cache_dir, 'pdb_structures')
            os.makedirs(pdb_dir, exist_ok=True)
            logger.info(f"Using PDB directory: {pdb_dir}")
            
            # Check if file already exists in cache
            expected_file = os.path.join(pdb_dir, f"{pdb_id.lower()}.pdb")
            if os.path.exists(expected_file):
                logger.info(f"Found cached PDB file: {expected_file}")
            else:
                logger.info(f"Downloading PDB file for {pdb_id}")
            
            # Fetch PDB file
            pdb_file = self.pdb_list.retrieve_pdb_file(
                pdb_id,
                pdir=pdb_dir,
                file_format='pdb'
            )
            logger.info(f"PDB file path: {pdb_file}")
            
            if not os.path.exists(pdb_file):
                raise FileNotFoundError(f"PDB file not found at {pdb_file}")
            
            # Parse structure
            logger.info(f"Parsing structure from file: {pdb_file}")
            structure = self.parser.get_structure(pdb_id, pdb_file)
            
            # Extract structure information
            chains = list(structure.get_chains())
            logger.info(f"Found {len(chains)} chains in structure")
            
            structure_info = {
                'pdb_id': pdb_id,
                'chains': len(chains),
                'residues': sum(1 for _ in structure.get_residues()),
                'atoms': sum(1 for _ in structure.get_atoms()),
                'resolution': self._get_resolution(structure),
                'chain_info': self._get_chain_info(chains),
                'file_path': pdb_file  # Add file path to response
            }
            
            logger.info(f"Successfully processed structure: {structure_info}")
            return structure_info
            
        except Exception as e:
            logger.error(f"Error fetching structure {pdb_id}: {str(e)}", exc_info=True)
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