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
Workflow A: Plots, Diagnostics, and PyDESeq2 Synthetic Dataset
Compute size factors and load synthetic data (cell f934279e):
loads information and raw counts for the dataset "synthetic" using pydeseq2.datasets.load_example_data.
DeseqDataSet(..., design="~condition", refit_cooks=True) and DefaultInference(n_cpus=8) are initialized.
examines dds.obs ("size_factors") and executes dds.fit_size_factors().
Set up the output directory, filter, and execute DESeq2 (cells 09d8c468 through 10532e9a):
constructs it using os.makedirs(OUTPUT_PATH, exist_ok=True) after setting OUTPUT_PATH = "../output_files/synthetic_example."
filters to genes with total counts >= 10 across samples (counts_df.sum(axis=0) >= 10) and eliminates samples with missing (NaN) condition values.
executes dds.deseq2().
Save AnnData and examine statistics, dispersions, and LFC shrinkage (cells ef396d6e through ffcf42db):
dds.to_picklable_anndata() is serialized to ../output_files/synthetic_example/result_adata.pkl.
prints dds.var ("dispersions").
executes ds.summary() and DeseqStats(dds, contrast=="condition", "B", "A", inference=inference).
examines dds.obsm."design_matrix"].columns and uses ds.lfc_shrink(coeff="condition[T.B]") to apply log2 fold-change shrinkage, saving the results in res_shrunk = ds.results_df.
Example of AnnData compilation and gene occurrence filtering (cells a2a88e67 and 1f833a39):
Genes with >= 10 are filtered
counts in a minimum of three samples (counts_df >= 10).sum(axis=0) >= 3) and creates an anndata.The synthetic dataset's AnnData object (adata.shape -> (100, 10)).
VST PCA plot, MA plot, sample diagnostics, and VST transformation (cells 07f527f0 and c0ea0d9a):
uses dds.vst (use_design=False, fit_type="parametric") to calculate variance-stabilized counts (dds.layers\vst_counts").
creates ma_plot.png (300 dpi) that classifies genes on a log-scale baseMean axis according to padj < 0.05 and log2FoldChange direction (Increased, Decreased, Not significant).
calculates the mean pairwise Pearson correlation between each sample and other samples (Mean_Sample_Corr), displays the sample diagnostics table, saves pca_vst.png (300 dpi), and performs PCA (2 components) on up to the top 500 highest-variance VST genes.
Workflow B: MyGene.info Gene Annotations (Ensembl BioMart Fallback)
Human gene annotations are received in cell 9a22c821 after Ensembl BioMart fails its TLS handshake in the notebook environment.
from MyGene.info ([https://mygene.info/v3/query](https://mygene.info/v3/query)).
List the Ensembl gene IDs for humans, such as ['ENSG00000125285', 'ENSG00000182968'].
[https://mygene.info/v3/query](https://mygene.info/v3/query) with q=f'ensembl.gene:{gene_id}, fields='symbol,entrezgene,go', and species='human'.
Gather distinct GO IDs from all hit.get('go', {}) categories, connect them using ;, and create a DataFrame with ensembl_gene_id, entrezgene_id, external_gene_name, and go_id.
Workflow C: Top-15 VST Gene Heatmap GSE224740
Within the gse224740_heatmap cell:
First, create upstream tables by running src/analyze_gse224740_propionate.py so that:
The project root (or a parent directory) contains results/tables/GSE224740_3day_vst_expression.csv, results/tables/GSE224740_propionate_vs_control_deseq2.csv, and results/tables/GSE224740_3day_library_sizes.csv.
Choose the top 15 DE genes by removing any missing padj values from de_results, sorting by padj in ascending order, and choosing the first 15 Ensembl_ID entries that are present as columns in vst_df.
Sample order and row-wise z-score VST expression:
uses sample_metadata ("condition") to arrange samples as "Control" and "Propionate" respectively.
subtracts the row mean of each gene, transposes VST values to genes × samples,
and divides by the population standard deviation of each gene (ddof=0, substituting NaN for 0).
Z-scores are saved to results/tables/GSE224740_top15_vst_gene_heatmap_zscores.csv.
Plot and save heatmap: Creates a coolwarm heatmap (vmin=-2, vmax=2) and saves results/figures/GSE224740_top15_vst_gene_heatmap.png at 300 dpi with a vertical divider at x=2.5 dividing Control and Propionate samples.
4. The Notebook's Notes and Cautions
Variable/state overwrites across cells: c0ea0d9a uses ds.results_df and dds from previous PyDESeq2 cells, whereas cell 1f833a39 reloads synthetic raw counts and metadata and adds a new adata variable. If you want to replicate the entire synthetic operation, run the cells sequentially.
In order to load the synthetic dataset and query MyGene, network connectivity is necessary.Internet access is necessary for installing or utilizing external enrichment sets with gseapy.
Heatmap sample boundary: ax.axvline(2.5,...) presupposes three control samples, then propionate samples (0, 1, 2 vs. 3, 4,5. Modify the divider position if each group in your GSE224740_3day_library_sizes.csv has a different count.
