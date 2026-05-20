"""Data processing functions for isoform entropy analysis.

This module handles:
- Entropy and correlation calculations
- Gene ranking
- Table data preparation
"""
import logging
import pandas as pd
import numpy as np

log = logging.getLogger(__name__)
from itertools import combinations
from scipy.stats import spearmanr


def calculate_entropy_and_correlation(df: pd.DataFrame, sample_cols, global_col: str):
    """Calculate entropy and Spearman correlation for each gene.
    
    Args:
        df: DataFrame with gene_id, transcript_id, and sample columns
        sample_cols: List of sample column names
        global_col: Name of global aggregation column
        
    Returns:
        List of dictionaries with results per gene
    """
    pairs = list(combinations(sample_cols, 2))
    results = []

    for gene, sub in df.groupby("gene_id", sort=False):
        row = {"gene_id": gene, "n_isoforms": len(sub)}

        top_idx = sub[global_col].idxmax()
        top_isoform = sub.loc[top_idx, sample_cols].astype(float)
        probs = top_isoform / top_isoform.sum()
        entropy = -np.sum(probs * np.log(probs + 1e-12))
        row["top_isoform_entropy"] = round(float(entropy), 4)

        sample_sums = sub[sample_cols].astype(float).sum(axis=0)
        sample_probs = sample_sums / sample_sums.sum()
        summed_entropy = -np.sum(sample_probs * np.log(sample_probs + 1e-12))
        row["summed_isoform_entropy"] = round(float(summed_entropy), 4)

        # Average expression: mean of all isoform values across all sample columns
        row["mean_expression"] = round(float(sub[sample_cols].astype(float).values.mean()), 4)

        for a, b in pairs:
            x = sub[a].tolist()
            y = sub[b].tolist()
            if len(set(x)) == 1 or len(set(y)) == 1:
                row[f"{a} vs {b}"] = None
            else:
                corr = spearmanr(x, y).correlation
                row[f"{a} vs {b}"] = round(float(corr), 4)
        results.append(row)
    return results


def compute_min_spearman_per_gene(results_df):
    """Compute minimum Spearman correlation across all sample pairs per gene.
    
    Args:
        results_df: DataFrame with correlation results
        
    Returns:
        List of minimum Spearman values per gene
    """
    corr_cols = [c for c in results_df.columns if " vs " in c]
    min_spearman = []
    for _, row in results_df.iterrows():
        spearman_vals = [row[c] for c in corr_cols if isinstance(row.get(c), (int, float))]
        min_spearman.append(float(min(spearman_vals)) if spearman_vals else None)
    return min_spearman


def compute_gene_ranking(results_df):
    """Rank genes based on min Spearman and entropy values.
    
    Ranking criteria:
    1. Negative min Spearman genes come before positive ones
    2. Within each group, rank by minimum of (top_isoform_entropy, summed_isoform_entropy)
       Lower min entropy = higher rank (lower rank number)
    
    Args:
        results_df: DataFrame with gene results
        
    Returns:
        List of rank numbers (1-indexed) corresponding to each row
    """
    df_work = results_df.copy()
    df_work["min_spearman"] = compute_min_spearman_per_gene(results_df)
    
    # Compute minimum of the two entropies
    df_work["min_entropy"] = df_work[["top_isoform_entropy", "summed_isoform_entropy"]].min(axis=1)
    
    # Separate into negative and positive Spearman groups
    df_work["is_negative_spearman"] = df_work["min_spearman"] < 0
    
    # Sort: negative Spearman first (True=1, False=0, so descending puts True first)
    # Then by min_entropy ascending (lower entropy = higher priority)
    df_work_sorted = df_work.sort_values(
        by=["is_negative_spearman", "min_entropy"],
        ascending=[False, True]
    )
    
    # Assign ranks (1-indexed)
    df_work_sorted["rank"] = range(1, len(df_work_sorted) + 1)
    
    # Map ranks back to original order
    rank_map = df_work_sorted.set_index(df_work_sorted.index)["rank"].to_dict()
    ranks = [rank_map[i] for i in df_work.index]
    
    return ranks


def compute_gene_ranking_by_expression(results_df):
    """Rank genes by mean expression (highest expression = rank 1).

    Genes with the highest average isoform expression across all samples
    appear at the top of the ranking.

    Args:
        results_df: DataFrame with gene results, must contain 'mean_expression'

    Returns:
        List of rank numbers (1-indexed) corresponding to each row
    """
    df_work = results_df[["mean_expression"]].copy()
    df_work_sorted = df_work.sort_values("mean_expression", ascending=False)
    df_work_sorted["rank_expr"] = range(1, len(df_work_sorted) + 1)
    rank_map = df_work_sorted["rank_expr"].to_dict()
    return [rank_map[i] for i in df_work.index]


def gene_has_cds(gene_id, isoforms_by_gene):
    """Check if a gene has any coding exons in any of its isoforms.
    
    Args:
        gene_id: Gene ID to check
        isoforms_by_gene: Dictionary mapping gene_id -> transcript_id -> list of exons
    
    Returns:
        True if any exon in any isoform has CDS, False otherwise
    """
    if gene_id not in isoforms_by_gene:
        return False
    
    isoforms = isoforms_by_gene[gene_id]
    for transcript_id, exons in isoforms.items():
        for exon in exons:
            if (exon.get("cds_start") is not None and 
                exon.get("cds_end") is not None and 
                exon["cds_end"] > exon["cds_start"]):
                return True
    
    return False


def prepare_table_data(results_df, isoforms_by_gene=None, gene_names=None,
                       af_geometry_mapping=None, default_ranking="spearman",
                       protein_sequences=None, domain_mapping=None):
    """Prepare data for the data table with rankings.

    Args:
        results_df: DataFrame with gene results
        isoforms_by_gene: Optional dictionary of exon structures
        gene_names: Optional dictionary mapping gene_id to gene_name
        af_geometry_mapping: Optional output of build_alphafold_geometry_mapping()
        default_ranking: Which ranking to sort the table by initially.
            'spearman' sorts by 'rank' (Spearman/entropy-based, default),
            'expression' sorts by 'rank_by_expression' (highest mean expression first).
        protein_sequences: Optional dict mapping transcript_id -> amino-acid sequence.
            When provided, a 'Max AA Length' column is added showing the longest
            protein sequence among all transcripts of each gene.
        domain_mapping: Optional dict mapping transcript_id -> list of InterPro domains.
            When provided, a 'Has Domains' column is added for genes that have
            at least one transcript with domain data.

    Returns:
        DataFrame formatted for display in data table
    """
    df_table = results_df[["gene_id", "n_isoforms", "top_isoform_entropy", "summed_isoform_entropy", "mean_expression"]].copy()
    df_table["min_spearman"] = compute_min_spearman_per_gene(results_df)
    df_table["rank"] = compute_gene_ranking(results_df)
    df_table["rank_by_expression"] = compute_gene_ranking_by_expression(results_df)

    # Add Gene Name column if gene_names is available
    if gene_names is not None:
        df_table["gene_name"] = df_table["gene_id"].apply(lambda g: gene_names.get(g, ""))
        df_table["gene_name_link"] = df_table.apply(
            lambda row: f'[{row["gene_name"]}](https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={row["gene_id"]})' if row["gene_name"] else "",
            axis=1
        )

    # Add Has CDS column if exon data is available
    if isoforms_by_gene is not None:
        df_table["has_cds"] = df_table["gene_id"].apply(lambda g: gene_has_cds(g, isoforms_by_gene))

    # Add Max AA Length column: longest protein sequence among all transcripts of the gene
    if isoforms_by_gene is not None and protein_sequences:
        def gene_max_aa_length(gene_id):
            if gene_id not in isoforms_by_gene:
                return None
            lengths = [
                len(protein_sequences[tid])
                for tid in isoforms_by_gene[gene_id].keys()
                if tid in protein_sequences
            ]
            return max(lengths) if lengths else None

        df_table["max_aa_length"] = df_table["gene_id"].apply(gene_max_aa_length)

    # Add Has 3D column: True if any transcript has AlphaFold geometry
    if isoforms_by_gene is not None and af_geometry_mapping:
        af_geometry_mapping = af_geometry_mapping or {}

        def gene_has_3d(gene_id):
            """Check if any transcript of this gene has AlphaFold geometry."""
            if gene_id not in isoforms_by_gene:
                return False
            for transcript_id in isoforms_by_gene[gene_id].keys():
                key = transcript_id.replace(".", "").lower()
                if key in af_geometry_mapping:
                    return True
            return False

        df_table["has_3d"] = df_table["gene_id"].apply(gene_has_3d)
        n_with_3d = df_table["has_3d"].sum()
        log.info("prepare_table_data: %d / %d genes have Has 3D = True (af_geometry_mapping has %d entries)",
                 n_with_3d, len(df_table), len(af_geometry_mapping))
    elif not af_geometry_mapping:
        log.warning("prepare_table_data: af_geometry_mapping is empty – Has 3D column will be absent")

    # Add Has Domains column: True if any transcript has InterPro domains
    if isoforms_by_gene is not None and domain_mapping:
        domain_mapping = domain_mapping or {}

        def gene_has_domains(gene_id):
            if gene_id not in isoforms_by_gene:
                return False
            for transcript_id in isoforms_by_gene[gene_id].keys():
                domains = domain_mapping.get(transcript_id)
                if domains:
                    return True
            return False

        df_table["has_domains"] = df_table["gene_id"].apply(gene_has_domains)
        n_with_domains = df_table["has_domains"].sum()
        log.info("prepare_table_data: %d / %d genes have Has Domains = True (domain_mapping has %d entries)",
                 n_with_domains, len(df_table), len(domain_mapping))
    elif domain_mapping is not None and not domain_mapping:
        log.warning("prepare_table_data: domain_mapping is empty – Has Domains column will be absent")
    
    sort_col = "rank_by_expression" if default_ranking == "expression" else "rank"
    df_table = df_table.sort_values(sort_col)
    
    df_table = df_table.rename(columns={
        "rank": "Rank",
        "gene_id": "Gene ID",
        "gene_name": "Gene Name",
        "gene_name_link": "Gene Name",
        "n_isoforms": "# Isoforms",
        "max_aa_length": "Max AA Length",
        "top_isoform_entropy": "Top Entropy",
        "summed_isoform_entropy": "Summed Entropy",
        "mean_expression": "Mean Expression",
        "min_spearman": "Min Spearman",
        "rank_by_expression": "Rank (Expr)",
        "has_cds": "Has CDS",
        "has_3d": "Has 3D",
        "has_domains": "Has Domains"
    })
    
    base_cols = ["Rank"]
    if "Gene Name" in df_table.columns:
        base_cols.append("Gene Name")
    base_cols.append("Gene ID")
    base_cols.append("# Isoforms")
    if "Max AA Length" in df_table.columns:
        base_cols.append("Max AA Length")
    if "Has CDS" in df_table.columns:
        base_cols.append("Has CDS")
    if "Has 3D" in df_table.columns:
        base_cols.append("Has 3D")
    if "Has Domains" in df_table.columns:
        base_cols.append("Has Domains")
    if "Rank (Expr)" in df_table.columns:
        base_cols.append("Rank (Expr)")
    base_cols.extend(["Mean Expression", "Min Spearman", "Top Entropy", "Summed Entropy"])
    
    df_table = df_table[base_cols]
    
    return df_table

