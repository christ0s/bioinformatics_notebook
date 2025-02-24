# Genome Analysis Web Tool

A web-based bioinformatics application for analyzing genomic sequences, performing BLAST searches, and visualizing sequence data.

## Features

- 🧬 Load and analyze DNA/RNA sequences from NCBI
- 🔍 Perform BLAST searches on protein sequences
- 📊 Visualize nucleotide and amino acid distributions
- 📋 Copy full DNA/RNA sequences
- 💾 Client-side caching of BLAST results
- 🎯 Identify top 5 proteins (>20 amino acids)

## Prerequisites

- [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Python 3.9+

## Installation

1. **Clone the repository**
```bash
git clone <repository-url>
cd Python_web_app_tool
```

2. **Create and activate conda environment**
```bash
# Create environment from yml file
conda env create -f environment.yml

# Activate the environment
conda activate genome_analyzer
```

## Configuration

1. **Set up NCBI email and API key**
Update `utils/blast.py`:
```python
Entrez.email = "your_email@example.com"  # Required by NCBI
Entrez.api_key = "your_api_key"          # Optional, but recommended
```

## Usage

1. **Activate conda environment (if not already activated)**
```bash
conda activate genome_analyzer
```

2. **Start the application**
```bash
python app.py
```

3. **Access the web interface**
- Open your browser and navigate to `http://localhost:5000`

4. **Upload sequence IDs**
- Create a text file with NCBI sequence IDs (one per line)
- Click "Choose File" and select your file
- Click "Upload" to process sequences

5. **Analyze Results**
- View sequence details (DNA, RNA, protein)
- Click on DNA/RNA sequences to view full sequence
- Use "Copy" buttons to copy sequences
- View amino acid counts
- Select proteins for BLAST search
- View distribution plots

## Directory Structure

```
/Python_web_app_tool/
├── templates/          # HTML templates
├── static/            # CSS and static files
├── utils/             # Python utility modules
├── uploads/           # Uploaded sequence files
├── cache/             # BLAST results cache
└── app.py            # Flask application
```

## Example Input File

```text
MN908947    # Example sequence ID
NC_001722.1 # Example sequence ID
```

## Notes

- BLAST searches are cached to improve performance
- Timeout set to 60 seconds for BLAST searches
- Plots are generated for nucleotide and amino acid distributions
- Top 5 proteins are filtered to include only those >20 amino acids

## Troubleshooting

- If BLAST search fails, verify your internet connection
- For "$ is not defined" errors, check if jQuery is properly loaded
- Ensure all conda dependencies are installed correctly
- Check NCBI email and API key configuration

## Deactivating the Environment

When you're done working with the tool:
```bash
conda deactivate
```

## License

[Your License Here]

## Contributors

[Your Name/Team]

## Contact

For issues and support, please [create an issue](repository-issues-url) or contact [your-email].