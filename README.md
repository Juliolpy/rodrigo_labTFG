# COLUMBO Design - Genetic Parts Designer

Program for designing genetic parts for the COLUMBO diagnostic system. This tool analyzes DNA sequences to find NGG PAM sites, designs molecular beacons, and creates primers for CRISPR-based diagnostics.

## 📚 Documentation

📖 **Full Documentation**: [https://juliolpy.github.io/rodrigo_labTFG/](https://juliolpy.github.io/rodrigo_labTFG/)

## 🚀 Quick Start with Docker

The easiest way to use COLUMBO Design is with Docker:

```bash
# Pull the Docker image
docker pull ghcr.io/juliolpy/rodrigo_labtfg:latest

# Run the container interactively
docker run -it ghcr.io/juliolpy/rodrigo_labtfg:latest

# Run with your own FASTA file
docker run -v $(pwd):/data ghcr.io/juliolpy/rodrigo_labtfg:latest columbo-cli /data/your_sequence.fasta
```

## 📋 Usage

### Inside the Docker Container

Once you're inside the container, you can use the following commands:

```bash
# Show help and available options
columbo-cli --help

# Analyze a FASTA file (JSON output by default)
columbo-cli <filename>.fasta

# Analyze a FASTA file with pickle output
columbo-cli <filename>.fasta --output pickle

# Alternative ways to run the CLI
python -m columbo_design.cli --help
poetry run columbo-cli --help
```

### Command Line Arguments

- `fasta`: Path to the input FASTA file (required)
- `--output`: Output format - "json" or "pickle" (default: "json")

### Example Output Files

The tool generates two output files:
- `output_NGG.json/pkl`: Contains NGG PAM sites and ColumboParts
- `output_PRIMERS.json/pkl`: Contains designed primers and scores

## 🔬 What COLUMBO Design Does

1. **Sequence Analysis**: Reads FASTA files and finds NGG PAM sites
2. **ColumboParts Creation**: Builds genetic parts around each PAM site
3. **Molecular Beacon Design**: Designs beacons for detection
4. **Primer Design**: Creates primers using Primer3
5. **Scoring**: Evaluates and ranks all components

## 🏗️ Local Development

If you prefer to run the project locally:

### Prerequisites

- Python 3.12+
- Poetry (for dependency management)

### Installation

```bash
# Clone the repository
git clone https://github.com/juliolpy/rodrigo_labTFG.git
cd rodrigo_labTFG

# Install dependencies
poetry install

# Run the CLI
poetry run columbo-cli --help
```

## 📁 Project Structure

- `src/columbo_design/`: Main source code
  - `cli.py`: Command-line interface
  - `core.py`: Core analysis functions
  - `primer_design.py`: Primer design functionality
  - `beacon.py`: Molecular beacon design
- `data/`: Sample data files and outputs
- `tests/`: Test files
- `Dockerfile`: Docker configuration

## 🧪 Dependencies

- **bio**: Bioinformatics library for sequence analysis
- **viennarna**: RNA folding and structure prediction
- **primer3-py**: Primer design algorithms
- **pytest**: Testing framework
- **sphinx**: Documentation generation

## 📊 Output Format

The tool generates detailed analysis including:
- NGG PAM site positions
- Protospacer sequences
- Molecular beacon designs with folding structures
- Primer pairs with scores
- Melting temperatures and hybridization energies

## 🤝 Contributing

This project is part of Julio's TFG (Final Degree Project). For contributions or questions, please contact the authors.

## 📄 License

Unlicense

## 👥 Authors

- Julio <jcochenta@gmail.com>
- luksgrin <lucas.goiriz@csic.es> _(as advisor)_
