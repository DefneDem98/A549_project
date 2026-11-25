"""Analyze the GSE224740 A549 3-day propionate-versus-control RNA-seq data."""

from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

PROJECT_ROOT = Path(__file__).resolve().parents[1]
COUNT_TABLE = PROJECT_ROOT / "data" / "raw" / "GSE224740_Raw_Count_Table_1.txt.gz"
SAMPLE_METADATA = (
    PROJECT_ROOT / "data" / "metadata" / "GSE224740_A549_3day_samples.csv"
)
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"
ANNOTATION_URL = "https://mygene.info/v3/query"
ANNOTATION_BATCH_SIZE = 100
SAMPLE_IDS = [f"GSM70306{i}" for i in range(51, 57)]
ALPHA = 0.05


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and validate raw counts and metadata for the six A549 3-day samples."""
    for path in (COUNT_TABLE, SAMPLE_METADATA):
        if not path.is_file():
            raise FileNotFoundError(f"Missing required file: {path}")

    sample_table = pd.read_csv(SAMPLE_METADATA).set_index("geo_sample_id")
    if sample_table.index.tolist() != SAMPLE_IDS:
        raise ValueError("Sample metadata must contain GSM7030651-GSM7030656 in order.")
    if sample_table["condition"].value_counts().to_dict() != {"Control": 3, "Propionate": 3}:
        raise ValueError("Expected three Control and three Propionate samples.")

    raw_table = pd.read_csv(COUNT_TABLE, sep="\t", compression="gzip")
    if "Ensembl_ID" not in raw_table.columns:
        raise KeyError("Expected an 'Ensembl_ID' column in the GEO count table.")

    count_cols = sample_table["count_column"].tolist()
    missing = set(count_cols) - set(raw_table.columns)
    if missing:
        raise KeyError(f"GEO count table is missing sample columns: {sorted(missing)}")

    gene_by_sample = raw_table.set_index("Ensembl_ID")[count_cols].apply(
        pd.to_numeric, errors="raise"
    )
    gene_by_sample.columns = sample_table.index

    values = gene_by_sample.to_numpy()
    if (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError("DESeq2 requires non-negative integer raw counts.")

    gene_by_sample = gene_by_sample.astype("int64")
    if not gene_by_sample.index.is_unique:
        gene_by_sample = gene_by_sample.groupby(level=0, sort=False).sum()

    counts = gene_by_sample.T
    counts.index.name = "sample"
    metadata = sample_table.loc[counts.index, ["condition"]].copy()
    metadata["condition"] = pd.Categorical(
        metadata["condition"], categories=["Control", "Propionate"], ordered=True
    )
    return counts, metadata


def filter_counts(
    counts: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Save library-size diagnostics and filter out low-count genes."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    library_sizes = pd.DataFrame(
        {
            "condition": metadata.loc[counts.index, "condition"].astype(str),
            "raw_library_size": counts.sum(axis=1).astype("int64"),
        },
        index=counts.index,
    )
    library_sizes.to_csv(TABLES_DIR / "GSE224740_3day_library_sizes.csv", index_label="sample")
    print("Raw library sizes before filtering:\n", library_sizes.to_string())

    colors = ["#4C78A8" if c == "Control" else "#F58518" for c in library_sizes["condition"]]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(library_sizes.index, library_sizes["raw_library_size"], color=colors)
    ax.set(
        ylabel="Raw counts per sample",
        xlabel="GEO sample",
        title="GSE224740: A549 3-day library sizes",
    )
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "GSE224740_3day_library_sizes.png", dpi=300)
    plt.close(fig)

    keep_genes = (counts >= 10).sum(axis=0) >= 3
    filtered = counts.loc[:, keep_genes].copy()
    if filtered.empty:
        raise ValueError("The low-expression filter removed every gene.")
    print(f"Low-count filter retained {filtered.shape[1]:,} of {counts.shape[1]:,} genes.")
    return filtered, library_sizes


def run_deseq2(
    counts_filtered: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[DeseqDataSet, pd.DataFrame]:
    """Fit PyDESeq2, compute Propionate vs Control contrast, and shrink LFCs."""
    inference = DefaultInference(n_cpus=4)
    dds = DeseqDataSet(
        counts=counts_filtered,
        metadata=metadata,
        design="~condition",
        refit_cooks=False,
        inference=inference,
    )
    if "condition[T.Propionate]" not in dds.obsm["design_matrix"].columns:
        raise ValueError("Control was not set as the reference condition.")

    dds.deseq2()
    stats = DeseqStats(
        dds, contrast=("condition", "Propionate", "Control"), alpha=ALPHA, inference=inference
    )
    stats.summary()

    results = stats.results_df.copy()
    results.index.name = "Ensembl_ID"
    sig = results["padj"].notna() & (results["padj"] < ALPHA)
    results["direction"] = np.select(
        [sig & (results["log2FoldChange"] > 0), sig & (results["log2FoldChange"] < 0)],
        ["increased_in_propionate", "decreased_in_propionate"],
        default="not_significant",
    )

    stats.lfc_shrink(coeff="condition[T.Propionate]")
    results["log2FoldChange_shrunk"] = stats.results_df["log2FoldChange"]
    results.to_csv(TABLES_DIR / "GSE224740_propionate_vs_control_deseq2.csv")

    print(f"Genes tested: {int(results['pvalue'].notna().sum()):,}")
    print(f"Increased (padj < {ALPHA}): {int((results['direction'] == 'increased_in_propionate').sum()):,}")
    print(f"Decreased (padj < {ALPHA}): {int((results['direction'] == 'decreased_in_propionate').sum()):,}")
    return dds, results


def save_ma_plot(results: pd.DataFrame) -> None:
    """Save an MA plot using shrunken log2 fold changes."""
    base_mean = results["baseMean"].to_numpy(dtype=float)
    lfc = results["log2FoldChange_shrunk"].to_numpy(dtype=float)
    finite = np.isfinite(base_mean) & np.isfinite(lfc)
    sig = (results["padj"] < ALPHA).to_numpy() & finite
    x_ma = np.log10(np.maximum(base_mean, 1e-8))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(x_ma[finite & ~sig], lfc[finite & ~sig], s=8, alpha=0.45, color="#888888", label="Not significant")
    ax.scatter(x_ma[sig], lfc[sig], s=10, alpha=0.8, color="#D62728", label=f"BH-adjusted p < {ALPHA}")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set(
        xlabel="log10 mean normalized count",
        ylabel="Shrunken log2 fold change (Propionate / Control)",
        title="GSE224740: Propionate vs Control at 3 days",
    )
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "GSE224740_propionate_vs_control_MA.png", dpi=300)
    plt.close(fig)


def _extract_ensembl_ids(hit: dict) -> list[str]:
    ensembl = hit.get("ensembl", {})
    raw_ids = (
        ensembl.get("gene", [])
        if isinstance(ensembl, dict)
        else [e.get("gene") for e in ensembl if isinstance(e, dict) and e.get("gene")]
        if isinstance(ensembl, list)
        else []
    )
    return [raw_ids] if isinstance(raw_ids, str) else list(raw_ids)


def fetch_gene_annotations(gene_ids: list[str]) -> pd.DataFrame:
    """Fetch symbols and Entrez IDs for Ensembl IDs from MyGene.info."""
    annotations: dict[str, dict[str, str | int | None]] = {}
    batch_set = set(gene_ids)
    with requests.Session() as session:
        for start in range(0, len(gene_ids), ANNOTATION_BATCH_SIZE):
            batch = gene_ids[start : start + ANNOTATION_BATCH_SIZE]
            resp = session.get(
                ANNOTATION_URL,
                params={
                    "q": " OR ".join(f"ensembl.gene:{gid}" for gid in batch),
                    "fields": "ensembl.gene,symbol,entrezgene",
                    "species": "human",
                    "size": len(batch),
                },
                timeout=30,
            )
            resp.raise_for_status()
            for hit in resp.json().get("hits", []):
                for gid in _extract_ensembl_ids(hit):
                    if gid in batch_set:
                        annotations[gid] = {
                            "symbol": hit.get("symbol"),
                            "entrezgene": hit.get("entrezgene"),
                        }

    missing = batch_set - annotations.keys()
    if missing:
        print(f"MyGene.info returned no annotation for {len(missing):,} genes.")
    return pd.DataFrame.from_dict(annotations, orient="index").rename_axis("Ensembl_ID")


def save_heatmap_matrix(
    matrix: pd.DataFrame,
    *,
    x_labels: list[str],
    y_labels: list[str],
    xlabel: str,
    ylabel: str,
    title: str,
    cbar_label: str,
    out_path: Path,
    vmin: float,
    vmax: float,
    annotate_cells: bool = False,
    group_divider_x: float | None = None,
    footnote: str | None = None,
    figsize: tuple[float, float] = (7, 6),
) -> None:
    """Plot a labeled matrix heatmap and save to disk."""
    fig, ax = plt.subplots(figsize=figsize)
    img = ax.imshow(matrix.to_numpy(), aspect="auto", interpolation="nearest", cmap="coolwarm", vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(len(x_labels)), labels=x_labels, rotation=45 if annotate_cells else 0, ha="right" if annotate_cells else "center")
    ax.set_yticks(np.arange(len(y_labels)), labels=y_labels)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    if group_divider_x is not None:
        ax.axvline(group_divider_x, color="black", linewidth=1)
    if annotate_cells:
        for r in range(matrix.shape[0]):
            for c in range(matrix.shape[1]):
                ax.text(c, r, f"{matrix.iloc[r, c]:.2f}", ha="center", va="center", color="black")
    fig.colorbar(img, ax=ax, label=cbar_label)
    if footnote:
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        fig.text(0.5, 0.02, footnote, ha="center", fontsize=9)
    else:
        fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_vst_and_heatmaps(
    dds: DeseqDataSet,
    results: pd.DataFrame,
    annotations: pd.DataFrame,
    metadata: pd.DataFrame,
    library_sizes: pd.DataFrame,
) -> None:
    """Compute VST expression and save PCA, sample correlation, and top-15 gene heatmaps."""
    if min(dds.shape) < 2:
        raise ValueError("VST diagnostics require at least 2 samples and 2 genes.")

    dds.vst(use_design=False, fit_type="parametric")
    vst = pd.DataFrame(
        np.asarray(dds.layers["vst_counts"]), index=dds.obs_names, columns=dds.var_names
    )
    vst.to_csv(TABLES_DIR / "GSE224740_3day_vst_expression.csv", index_label="sample")

    top_var = vst.var(axis=0).nlargest(min(500, vst.shape[1])).index
    pca = PCA(n_components=2)
    coords = pca.fit_transform(vst.loc[:, top_var])
    pca_scores = pd.DataFrame(
        {
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
            "condition": metadata.loc[vst.index, "condition"].astype(str).to_numpy(),
        },
        index=vst.index,
    )
    pca_scores.to_csv(TABLES_DIR / "GSE224740_3day_vst_pca_scores.csv", index_label="sample")
    print("VST-PCA scores and raw library sizes:\n", pca_scores.join(library_sizes["raw_library_size"]).to_string())
    print("PCA variance explained (%):", (100 * pca.explained_variance_ratio_).round(2).tolist())

    fig, ax = plt.subplots(figsize=(7, 5.5))
    for cond, color in [("Control", "#4C78A8"), ("Propionate", "#F58518")]:
        sub = pca_scores[pca_scores["condition"] == cond]
        ax.scatter(sub["PC1"], sub["PC2"], s=90, color=color, label=cond)
        for sample, row in sub.iterrows():
            ax.annotate(sample, (row["PC1"], row["PC2"]), xytext=(5, 5), textcoords="offset points", fontsize=8)
    ax.set(
        xlabel=f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}%)",
        ylabel=f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}%)",
        title="PCA of VST expression: GSE224740, 3 days",
    )
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "GSE224740_3day_vst_PCA.png", dpi=300)
    plt.close(fig)

    corr = vst.T.corr(method="pearson")
    print("VST sample correlation matrix:\n", corr.round(3).to_string())
    save_heatmap_matrix(
        corr,
        x_labels=list(corr.columns),
        y_labels=list(corr.index),
        xlabel="",
        ylabel="",
        title="VST sample correlation: GSE224740, 3 days",
        cbar_label="Pearson correlation",
        out_path=FIGURES_DIR / "GSE224740_3day_vst_sample_correlations.png",
        vmin=0.8,
        vmax=1.0,
        annotate_cells=True,
    )

    top_genes = results.dropna(subset=["padj"]).sort_values("padj").index[:15].tolist()
    if len(top_genes) < 10:
        raise ValueError(f"Expected 10-15 testable genes for the heatmap, got {len(top_genes)}.")

    gene_mat = vst.loc[:, top_genes].T
    row_sd = gene_mat.std(axis=1, ddof=0).replace(0, np.nan)
    z_scores = gene_mat.sub(gene_mat.mean(axis=1), axis=0).div(row_sd, axis=0)
    z_scores.to_csv(TABLES_DIR / "GSE224740_top15_vst_gene_heatmap_zscores.csv", index_label="Ensembl_ID")

    symbols = annotations.reindex(top_genes)["symbol"]
    y_labels = symbols.where(symbols.notna(), top_genes).tolist()
    x_labels = [f"{s}\n{metadata.loc[s, 'condition']}" for s in z_scores.columns]
    save_heatmap_matrix(
        z_scores,
        x_labels=x_labels,
        y_labels=y_labels,
        xlabel="Samples (Control followed by Propionate)",
        ylabel="Gene symbol",
        title="Top 15 differentially expressed genes: GSE224740, 3 days",
        cbar_label="Gene-wise z-score of VST expression",
        out_path=FIGURES_DIR / "GSE224740_top15_vst_gene_heatmap.png",
        vmin=-2,
        vmax=2,
        group_divider_x=2.5,
        footnote="PyDESeq2 VST values; each gene (row) is z-scored across the six samples.",
        figsize=(9, 7),
    )


def run_enrichment(results: pd.DataFrame, annotations: pd.DataFrame) -> None:
    """Run directional GO, KEGG, and Hallmark over-representation tests."""
    gene_sets_plan = [
        ("enrichment", ["GO_Biological_Process_2023", "KEGG_2021_Human"]),
        ("hallmark_enrichment", "MSigDB_Hallmark_2020"),
    ]
    annotated = results.join(annotations, how="left")

    for direction in ("increased_in_propionate", "decreased_in_propionate"):
        symbols = (
            annotated.loc[annotated["direction"] == direction, "symbol"]
            .dropna()
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
        if not symbols:
            print(f"No annotated {direction} genes; enrichment was not run.")
            continue

        for suffix, gene_sets in gene_sets_plan:
            enr = gp.enrichr(
                gene_list=symbols,
                gene_sets=gene_sets,
                organism="human",
                outdir=None,
                cutoff=1.0,
                no_plot=True,
            )
            enr.results.to_csv(TABLES_DIR / f"GSE224740_{direction}_{suffix}.csv", index=False)
            print(
                f"{suffix.replace('_', ' ').title()} for {direction}: "
                f"{len(symbols):,} annotated genes; {len(enr.results):,} terms."
            )


def main() -> None:
    """Run the full GSE224740 A549 3-day RNA-seq pipeline."""
    counts, metadata = load_inputs()
    counts_filtered, library_sizes = filter_counts(counts, metadata)
    dds, results = run_deseq2(counts_filtered, metadata)
    save_ma_plot(results)

    sig_gene_ids = results.index[results["direction"] != "not_significant"].tolist()
    annotations = fetch_gene_annotations(sig_gene_ids)
    annotations.to_csv(
        TABLES_DIR / "GSE224740_significant_gene_annotations.csv", index_label="Ensembl_ID"
    )

    save_vst_and_heatmaps(dds, results, annotations, metadata, library_sizes)
    run_enrichment(results, annotations)

    print(f"DE results: {TABLES_DIR / 'GSE224740_propionate_vs_control_deseq2.csv'}")
    print(f"Figures: {FIGURES_DIR}")


if __name__ == "__main__":
    main()