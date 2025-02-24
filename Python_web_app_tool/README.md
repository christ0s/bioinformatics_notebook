# Genome Analysis Web Tool

An interactive web application for genomic sequence analysis, BLAST searching, and visualization built with Flask and Biopython.

## 🧬 Features

- Load and analyze DNA/RNA sequences from NCBI databases
- Perform BLAST searches on protein sequences 
- Visualize nucleotide and amino acid distributions
- Copy full DNA/RNA sequences with one click
- Cache BLAST search results for better performance
- Identify top 5 proteins with length >20 amino acids
- Dark mode UI with interactive plots

## 🚀 Quick Start

### Prerequisites

```bash
# Install Python 3.8+ and pip
sudo apt-get update
sudo apt-get install python3 python3-pip
```

### Installation

```bash
# Clone the repository
git clone <repository-url>

# Navigate to the project directory
cd bioinformatics_notebook/Python_web_app_tool

# Install required Python packages
pip install -r requirements.txt
```

### Running the Application

```bash
# Run the Flask application
export FLASK_APP=app.py
flask run
```

### Accessing the Application

Open your web browser and go to `http://127.0.0.1:5000` to access the Genome Analysis Web Tool.

## 📂 Project Structure

```
bioinformatics_notebook/
│
├── Python_web_app_tool/
│   ├── app.py
│   ├── requirements.txt
│   ├── static/
│   ├── templates/
│   └── README.md
│
└── Accessing_NCBI_databases/
    └── All useful python or R scripts for bioinformatics analysis
```

## 🛠️ Technologies Used

- Flask
- Biopython
- HTML/CSS/JavaScript
- Bootstrap
- jQuery
## 🔧 Troubleshooting
Common issues and solutions:

BLAST Search Fails: Check internet connection and NCBI credentials
Loading Issues: Ensure all dependencies are installed
Display Problems: Clear browser cache if plots don't show

## 👥 Contributing
Contributions are welcome! Please read our contributing guidelines and submit pull requests.

## 📄 License
[MIT] - See LICENSE file for details

📮 Contact
For support or queries: 

Create an issue on GitHub
Email: [chrissanthis@gmai.com]
Project Link: [repository-url]