"""
Shared utility functions for isoform distribution analysis.
Used by both table generation and plotting scripts.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

def get_multi_isoform_genes(df_isoform_matrix):
    gene_isoform_counts = df_isoform_matrix.groupby('gene_id').size()
    return gene_isoform_counts[gene_isoform_counts > 1].index.tolist()

def parse_sample_name(sample_name, sep='_', region_index=0, condition_index=1):
    """
    Parse sample name to extract region and condition tokens.

    Args:
        sample_name: Sample name string (format: region_condition_...)
        sep: Token separator used in sample names
        region_index: Token index for the region value
        condition_index: Token index for the condition value

    Returns:
        Tuple of (region, condition) or (None, None) if parsing fails
    """
    parts = sample_name.split(sep)
    if len(parts) > max(region_index, condition_index):
        region = parts[region_index]
        condition = parts[condition_index]
        return region, condition
    return None, None


def _normalize_gene_samples(gene_transcripts, sample_cols):
    """
    Return a DataFrame where each sample column is normalized within the gene:
    for the isoforms of the gene, values per sample sum to 1.
    
    Args:
        gene_transcripts: DataFrame of transcript expressions for a single gene
        sample_cols: List of sample column names
        
    Returns:
        DataFrame with normalized values (sum to 1 per sample)
    """
    mat = gene_transcripts[sample_cols].copy()
    # Ensure numeric dtype to avoid object downcasting warnings
    mat = mat.apply(pd.to_numeric, errors='coerce')
    col_sums = mat.sum(axis=0)

    # Avoid division by zero: if a sample has total 0 for this gene, keep zeros
    safe_sums = col_sums.replace(0, np.nan)
    normalized = mat.divide(safe_sums, axis=1)

    # Fill NaNs with 0 while keeping float dtype, then infer any remaining objects
    normalized = normalized.fillna(0.0).infer_objects(copy=False)
    return normalized


def aggregate_samples_by_group(df_isoform_matrix, sample_cols, group_by='region', stat='sum',
                               sample_name_sep='_', region_index=0, condition_index=1):
    """
    Aggregate samples by brain region or condition.
    
    Args:
        df_isoform_matrix: DataFrame with isoform expression data
        sample_cols: List of sample column names
    group_by: 'region' or 'condition'
        stat: 'sum', 'mean', or 'normalized'
    sample_name_sep: Separator used in sample names
    region_index: Token index for the region value
    condition_index: Token index for the condition value
        
    Returns:
        DataFrame with aggregated data
    """
    aggregated_data = defaultdict(list)

    # Group samples by region/condition
    for sample_col in sample_cols:
        region, condition = parse_sample_name(
            sample_col,
            sep=sample_name_sep,
            region_index=region_index,
            condition_index=condition_index
        )
        if region and condition:
            group_key = region if group_by == 'region' else condition
            aggregated_data[group_key].append(sample_col)

    # Aggregate data
    result_data = {}
    for group_key, group_sample_cols in aggregated_data.items():
        if stat == 'mean':
            result_data[group_key] = df_isoform_matrix[group_sample_cols].mean(axis=1)
        elif stat == 'normalized':
            # Compute per-gene normalized values per sample, then aggregate by mean
            # Build a series aligned to df index
            vals = []
            for gene_id, gene_block in df_isoform_matrix.groupby('gene_id'):
                gb = gene_block.drop(columns=['gene_id'])
                norm = _normalize_gene_samples(gb, group_sample_cols).mean(axis=1)
                vals.append(norm)
            result_data[group_key] = pd.concat(vals).reindex(df_isoform_matrix.index)
        else:
            result_data[group_key] = df_isoform_matrix[group_sample_cols].sum(axis=1)

    result_df = pd.DataFrame(result_data)
    result_df['gene_id'] = df_isoform_matrix['gene_id']
    return result_df


def get_filtered_isoforms(gene_transcripts, sample_cols, cutoff_pct=2, stat='sum'):
    """
    Filter isoforms for a gene based on their contribution to global expression.
    Uses sum or mean depending on stat.
    
    Args:
        gene_transcripts: DataFrame of transcript expressions for a single gene
        sample_cols: List of sample column names
        cutoff_pct: Minimum percentage contribution to keep an isoform
        stat: 'sum', 'mean', or 'normalized'
        
    Returns:
        List of filtered isoform IDs, or None if no multi-isoform gene after filtering
    """
    if stat == 'mean':
        global_expr = gene_transcripts[sample_cols].mean(axis=1)
    elif stat == 'normalized':
        norm = _normalize_gene_samples(gene_transcripts, sample_cols)
        global_expr = norm.mean(axis=1)  # mean normalized contribution across samples
    else:
        global_expr = gene_transcripts[sample_cols].sum(axis=1)

    total_gene_expr = global_expr.sum()
    if total_gene_expr == 0:
        return None

    # Filter out isoforms contributing less than cutoff_pct to global expression
    global_expr_pct = (global_expr / total_gene_expr * 100)
    filtered_isoforms = global_expr_pct[global_expr_pct >= cutoff_pct].index.tolist()
    if len(filtered_isoforms) <= 1:
        return None
    return filtered_isoforms


def load_confidence_intervals(ci_path):
    """
    Read confidence_intervals.tsv and return dict: isoform -> {ci_lower, ci_upper}.
    
    Args:
        ci_path: Path to confidence intervals file, or None
        
    Returns:
        Dictionary mapping isoform IDs to confidence interval data
    """
    if not ci_path:
        return {}
    ci_path = Path(ci_path)
    if not ci_path.exists():
        return {}
    ci_df = pd.read_csv(ci_path, sep="\t")
    # Accept either index or column named 'isoform'
    if 'isoform' in ci_df.columns:
        ci_df = ci_df.set_index('isoform')
    ci_map = ci_df[['ci_lower', 'ci_upper']].to_dict(orient='index')
    return ci_map


def prepare_gene_list_and_paths(df_isoform_matrix, out_dir):
    """
    Helper function to prepare gene list and create output directory structure.
    
    Args:
        df_isoform_matrix: DataFrame with isoform expression data
        out_dir: Output directory path
        
    Returns:
        Tuple of (out_dir_path, gene_list)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    all_genes = get_multi_isoform_genes(df_isoform_matrix)
    return out_dir, all_genes
