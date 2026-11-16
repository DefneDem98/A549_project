# A549 RNA-Seq Analysis Project

## 1. Research Question
Does propionate treatment produce a transcriptional pattern consistent with a more epithelial state in A549 cells?
## 2. Working Directory & Folder Structure
All commands and paths are relative to the project root directory:
A549_project/
├── data/
│   ├── raw/                  # Unprocessed raw data files (FASTQ, raw counts, metadata)
│   └── processed/            # Cleaned, normalized data and AnnData objects
├── notebooks/
│   └── A549_analysis.ipynb   # Main interactive Jupyter notebook for exploratory analysis
├── src/                      # Custom Python scripts and modular functions
├── results/
│   ├── tables/               # Exported CSV/TSV tables (differential expression, GSEA results)
│   └── figures/              # Generated plots, heatmaps, and volcano plots
├── reports/                  # Final compiled analytical reports
├── .gitignore                # Git exclusion rules (virtual environments, caches, secrets)
├── requirements.txt          # Pinned Python package dependencies
└── README.md                 # Project documentation

## 3. Data Sources
* **Raw Data:** Deposited locally in `data/raw/` (excluded from version control). Includes raw count matrices and sample clinical/experimental metadata.
* **Processed Data:** Cleaned expression matrices and AnnData objects stored in `data/processed/`.

## 4. Environment Setup
The project runs on **Python 3.15.0** within an isolated virtual environment (`.venv`).

### Installation Steps:
1. Create the virtual environment:
   ```bash
   python3 -m venv .venv
