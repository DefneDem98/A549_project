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

## Installation & Execution Steps

### Workflow A: Plots, Diagnostics, and PyDESeq2 Synthetic Dataset

1. **Compute size factors and load synthetic data**
   * Loads metadata and raw counts for the `"synthetic"` dataset using `pydeseq2.datasets.load_example_data`.
   * Initializes `DeseqDataSet(..., design="~condition", refit_cooks=True)` and `DefaultInference(n_cpus=8)`.
   * Examines `dds.obs["size_factors"]` and executes `dds.fit_size_factors()`.

2. **Set up the output directory, filter, and execute DESeq2** *(cells `09d8c468` – `10532e9a`)*
   * Sets `OUTPUT_PATH = "../output_files/synthetic_example/"` and creates it using `os.makedirs(OUTPUT_PATH, exist_ok=True)`[cite: 8].
   * Filters to genes with total counts `>= 10` across samples (`counts_df.sum(axis=0) >= 10`) and eliminates samples with missing (`NaN`) condition values[cite: 8].
   * Executes `dds.deseq2()`[cite: 8].

3. **Save `AnnData` and examine statistics, dispersions, and LFC shrinkage** *(cells `ef396d6e` – `ffcf42db`)*
   * Serializes `dds.to_picklable_anndata()` to `../output_files/synthetic_example/result_adata.pkl`[cite: 8].
   * Prints `dds.var["dispersions"]` and executes `ds.summary()` using `DeseqStats(dds, contrast=["condition", "B", "A"], inference=inference)`[cite: 8].
   * Examines `dds.obsm["design_matrix"].columns` and applies log2 fold-change shrinkage via `ds.lfc_shrink(coeff="condition[T.B]")`, storing results in `res_shrunk = ds.results_df`[cite: 8].

4. **Example of `AnnData` compilation and gene occurrence filtering** *(cells `a2a88e67` and `1f833a39`)*
   * Filters genes with `>= 10` counts in a minimum of three samples (`(counts_df >= 10).sum(axis=0) >= 3`)[cite: 8].
   * Creates an `anndata.AnnData` object from the synthetic dataset (`adata.shape -> (100, 10)`)[cite: 8].

5. **VST PCA plot, MA plot, sample diagnostics, and VST transformation** *(cells `07f527f0` and `c0ea0d9a`)*
   * Computes variance-stabilized counts (`dds.layers["vst_counts"]`) via `dds.vst(use_design=False, fit_type="parametric")`[cite: 8].
   * Creates `ma_plot.png` (300 dpi) classifying genes on a log-scale `baseMean` axis according to `padj < 0.05` and `log2FoldChange` direction (`Increased`, `Decreased`, `Not significant`)[cite: 8].
   * Calculates each sample's mean pairwise Pearson correlation with other samples (`Mean_Sample_Corr`), displays the sample diagnostics table, performs PCA (2 components) on up to the top 500 highest-variance VST genes, and saves `pca_vst.png` (300 dpi)[cite: 8].

---

### Workflow B: MyGene.info Gene Annotations (Ensembl BioMart Fallback)

* In cell `9a22c821`, human gene annotations are retrieved from **MyGene.info** (`https://mygene.info/v3/query`) after Ensembl BioMart fails its TLS handshake in the notebook environment[cite: 8].
* Queries `https://mygene.info/v3/query` for human Ensembl gene IDs (e.g., `['ENSG00000125285', 'ENSG00000182968']`) using:
  * `q=f'ensembl.gene:{gene_id}'`[cite: 8]
  * `fields='symbol,entrezgene,go'`[cite: 8]
  * `species='human'`[cite: 8]
* Gathers distinct GO IDs from all `hit.get('go', {})` categories, joins them using `;`, and builds a `DataFrame` with columns `ensembl_gene_id`, `entrezgene_id`, `external_gene_name`, and `go_id`[cite: 8].

---

### Workflow C: Top-15 VST Gene Heatmap (`GSE224740`)

Within the `gse224740_heatmap` cell:

1. **Generate upstream tables first** by running `src/analyze_gse224740_propionate.py` so the required CSV inputs exist (`GSE224740_3day_vst_expression.csv`, `GSE224740_propionate_vs_control_deseq2.csv`, and `GSE224740_3day_library_sizes.csv`)[cite: 8].
2. **Select top 15 DE genes**, order samples (`Control` followed by `Propionate`), and compute row-wise z-scores of VST expression[cite: 8].
3. **Save outputs** to `results/tables/GSE224740_top15_vst_gene_heatmap_zscores.csv` and `results/figures/GSE224740_top15_vst_gene_heatmap.png`[cite: 8].
4. ## Data Setup & Downloading Raw Counts

This analysis uses the six A549 3-day samples from **GEO series `GSE224740`**:
- **Control (`GSM7030651`–`GSM7030653`)**
- **Propionate (`GSM7030654`–`GSM7030656`)**

### 1. Download `GSE224740_Raw_Count_Table_1.txt.gz`
1. Open the **GSE224740** record on NCBI GEO (`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224740`) and scroll to the **Supplementary file** section at the bottom.
2. Download **`GSE224740_Raw_Count_Table_1.txt.gz`** (keep it compressed as `.txt.gz`; do not extract it).
   - *(Optional)* If you have the direct GEO FTP/HTTPS supplementary file URL, you can also download it from the command line:
     ```bash
     curl -L "<DIRECT_GEO_SUPPLEMENTARY_URL_FOR_GSE224740_Raw_Count_Table_1.txt.gz>" \
       -o A549_project/data/raw/GSE224740_Raw_Count_Table_1.txt.gz
     ```

### 2. Place the File in `data/raw/`
Move or copy `GSE224740_Raw_Count_Table_1.txt.gz` into your raw data directory:

```text
A549_project/
├── data/
│   ├── metadata/
│   │   └── GSE224740_A549_3day_samples.csv
│   └── raw/
│       └── GSE224740_Raw_Count_Table_1.txt.gz
├── results/
│   ├── figures/
│   └── tables/
└── src/
    └── analyze_gse224740_propionate.py
