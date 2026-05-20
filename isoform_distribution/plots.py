"""
Isoform distribution plotting script.
Generates PDF visualizations of isoform expression patterns.
Uses shared utilities from utils.py.
"""

import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from pathlib import Path
import utilities.gene_expression_analysis as gene_expression_analysis
from .utils import (
    aggregate_samples_by_group,
    get_filtered_isoforms,
    load_confidence_intervals,
    prepare_gene_list_and_paths,
    _normalize_gene_samples
)

# Define base directory relative to this script
BASE_DIR = Path(__file__).parent.parent  # source/


def plot_gene_isoforms(axes, gene_transcripts, gene_id, sample_cols, start_row, ncols,
                       cutoff_pct=2, log_scale=False, stat='sum', ci_map=None):
    """
    Plot isoform distributions for a single gene across multiple samples.
    
    Args:
        axes: Matplotlib axes array
        gene_transcripts: DataFrame of transcript expressions for the gene
        gene_id: Gene identifier
        sample_cols: List of sample column names
        start_row: Starting row index in the axes array
        ncols: Number of columns in the subplot grid
        cutoff_pct: Minimum percentage contribution to keep an isoform
        log_scale: Whether to use log scale for y-axis
        stat: 'sum', 'mean', or 'normalized'
        ci_map: Optional confidence interval mapping
        
    Returns:
        True if gene was plotted, False otherwise
    """
    filtered_isoforms = get_filtered_isoforms(gene_transcripts, sample_cols, cutoff_pct, stat)
    if filtered_isoforms is None:
        return False

    # Global: sum, mean, or mean of normalized
    if stat == 'mean':
        global_expr = gene_transcripts[sample_cols].mean(axis=1).loc[filtered_isoforms]
    elif stat == 'normalized':
        norm = _normalize_gene_samples(gene_transcripts, sample_cols)
        global_expr = norm.mean(axis=1).loc[filtered_isoforms]
    else:
        global_expr = gene_transcripts[sample_cols].sum(axis=1).loc[filtered_isoforms]

    # Sort by global
    sorted_isoforms = [isoform for _, isoform in sorted(zip(global_expr, filtered_isoforms), reverse=True)]
    global_expr_sorted = global_expr.loc[sorted_isoforms]

    for sample_idx, sample_col in enumerate(sample_cols):
        row_offset = sample_idx // ncols
        col_idx = sample_idx % ncols
        row_idx = start_row + row_offset
        ax = axes[row_idx, col_idx]

        if stat == 'normalized':
            sample_matrix = _normalize_gene_samples(gene_transcripts, [sample_col])
            sample_expr = sample_matrix[sample_col].loc[sorted_isoforms]
        else:
            sample_expr = gene_transcripts[sample_col].loc[sorted_isoforms]

        if col_idx == 0 and row_offset == 0:
            ax.set_ylabel(gene_id, fontsize=9, fontweight='bold')
        ax.set_title(sample_col, fontsize=7)
        x = range(1, len(sorted_isoforms) + 1)

        # Global bars
        ax.bar(x, global_expr_sorted.values, color='lightgrey', label='Global', zorder=0)
        # Sample bars
        ax.bar(x, sample_expr.values, color='#4C78A8', label='Sample', zorder=1)
        ax.set_xticks([])
        ax.tick_params(axis='y', labelsize=6)

        # CI indicators on grey Global bars only for mean mode
        if ci_map and stat == 'mean':
            for i, iso in enumerate(sorted_isoforms, start=1):
                ci = ci_map.get(iso)
                if ci:
                    ax.vlines(i, ci['ci_lower'], ci['ci_upper'], colors='black', linewidth=1.0, alpha=0.8, zorder=2)
                    ax.plot([i], [ci['mean']], marker='_', color='black', markersize=8, zorder=3)

        if log_scale and stat != 'normalized':
            # Normalized values are in [0,1]; log scale is not meaningful here
            ax.set_yscale('log')

    return True


def turn_off_unused_axes(axes, start_row, last_sample_idx, ncols):
    """
    Turn off unused axes in the subplot grid.
    
    Args:
        axes: Matplotlib axes array
        start_row: Starting row index
        last_sample_idx: Index of the last sample plotted
        ncols: Number of columns in the grid
    """
    last_row_offset = last_sample_idx // ncols
    last_col = last_sample_idx % ncols
    last_row = start_row + last_row_offset
    for col_idx in range(last_col + 1, ncols):
        axes[last_row, col_idx].axis('off')


def create_pdf_pages(df_isoform_matrix, sample_cols, all_genes, pdf_path, ncols,
                     cutoff_pct=2, log_scale=False, stat='sum', ci_map=None):
    """
    Create multi-page PDF with isoform distribution plots.
    
    Args:
        df_isoform_matrix: DataFrame with isoform expression data
        sample_cols: List of sample column names
        all_genes: List of gene IDs to plot
        pdf_path: Output PDF file path
        ncols: Number of columns in the subplot grid
        cutoff_pct: Minimum percentage contribution to keep an isoform
        log_scale: Whether to use log scale for y-axis
        stat: 'sum', 'mean', or 'normalized'
        ci_map: Optional confidence interval mapping
    """
    rows_per_gene = math.ceil(len(sample_cols) / ncols)
    genes_per_page = max(1, 10 // rows_per_gene)
    pages = math.ceil(len(all_genes) / genes_per_page)
    Path(pdf_path).parent.mkdir(parents=True, exist_ok=True)
    
    with PdfPages(pdf_path) as pdf:
        for p in range(pages):
            gene_chunk = all_genes[p*genes_per_page:(p+1)*genes_per_page]
            nrows = len(gene_chunk) * rows_per_gene
            fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*2.2, nrows*2.2))
            if nrows == 1:
                axes = axes.reshape(1, -1)
            
            for gene_idx, gene_id in enumerate(gene_chunk):
                gene_transcripts = df_isoform_matrix[df_isoform_matrix['gene_id'] == gene_id].drop(columns=['gene_id'])
                start_row = gene_idx * rows_per_gene
                was_plotted = plot_gene_isoforms(
                    axes, gene_transcripts, gene_id, sample_cols, start_row, ncols,
                    cutoff_pct=cutoff_pct, log_scale=log_scale, stat=stat, ci_map=ci_map
                )
                if was_plotted:
                    turn_off_unused_axes(axes, start_row, len(sample_cols) - 1, ncols)
                else:
                    # Turn off all axes for genes that weren't plotted
                    for row_offset in range(rows_per_gene):
                        for col_idx in range(ncols):
                            ax = axes[start_row + row_offset, col_idx]
                            ax.axis('off')
                    # Add gene ID as grayed out text
                    ax = axes[start_row, 0]
                    ax.text(-0.2, 0.5, gene_id,
                           transform=ax.transAxes,
                           fontsize=9, fontweight='bold',
                           verticalalignment='center',
                           horizontalalignment='right',
                           color='gray', style='italic',
                           rotation=90)
            
            fig.suptitle(f'Isoform distributions by gene (page {p+1}/{pages})', fontsize=12)
            fig.tight_layout(rect=[0, 0.03, 1, 0.97])
            pdf.savefig(fig)
            plt.close(fig)


def generate_individual_plots(df_isoform_matrix, sample_cols, out_dir, ncols=13, 
                              cutoff_pct=2, log_scale=False, 
                              stat='sum', ci_map=None):
    """
    Generate individual sample plots for each gene.
    
    Args:
        df_isoform_matrix: DataFrame with isoform expression data
        sample_cols: List of sample column names
        out_dir: Output directory path
        ncols: Number of columns in the subplot grid
        cutoff_pct: Minimum percentage contribution to keep an isoform
        log_scale: Whether to use log scale for y-axis
        stat: 'sum', 'mean', or 'normalized'
        ci_map: Optional confidence interval mapping
    """
    out_dir, all_genes = prepare_gene_list_and_paths(df_isoform_matrix, out_dir)
    
    genes_suffix = f"_{len(all_genes)}genes"
    scale_suffix = "_log" if log_scale else ""
    pdf_path = out_dir / f'isoform_distributions_by_gene_{stat}_cutoff{cutoff_pct}pct{genes_suffix}{scale_suffix}.pdf'
    
    create_pdf_pages(df_isoform_matrix, sample_cols, all_genes, pdf_path, ncols, 
                     cutoff_pct=cutoff_pct, log_scale=log_scale, stat=stat, ci_map=ci_map)
    print(f"Wrote {pdf_path}")


def generate_aggregated_plots(df_isoform_matrix, sample_cols, out_dir, cutoff_pct=2, 
                              log_scale=False, stat='sum', ci_map=None):
    """
    Create two PDFs: one aggregated by brain region, one by condition.
    
    Args:
        df_isoform_matrix: DataFrame with isoform expression data
        sample_cols: List of sample column names
        out_dir: Output directory path
        cutoff_pct: Minimum percentage contribution to keep an isoform
        log_scale: Whether to use log scale for y-axis
        stat: 'sum', 'mean', or 'normalized'
        ci_map: Optional confidence interval mapping
    """
    out_dir, all_genes = prepare_gene_list_and_paths(df_isoform_matrix, out_dir)
    
    genes_suffix = f"_{len(all_genes)}genes"
    scale_suffix = "_log" if log_scale else ""
    
    for group_by in ['region', 'condition']:
        print(f"Aggregating by {group_by} with stat={stat}...")
        df_aggregated = aggregate_samples_by_group(df_isoform_matrix, sample_cols, 
                                                   group_by=group_by, stat=stat)
        agg_cols = [col for col in df_aggregated.columns if col != 'gene_id']
        
        pdf_path = out_dir / f'isoform_distributions_by_{group_by}_{stat}_cutoff{cutoff_pct}pct{genes_suffix}{scale_suffix}.pdf'
        create_pdf_pages(df_aggregated, agg_cols, all_genes, pdf_path, ncols=len(agg_cols), 
                        cutoff_pct=cutoff_pct, log_scale=log_scale, stat=stat, ci_map=ci_map)
        print(f"Wrote {pdf_path}")


def main():
    """Main entry point for plotting script."""
    parser = argparse.ArgumentParser(description='Generate isoform distribution plots')
    parser.add_argument('--cutoff-pct', type=float, default=1.5, 
                       help='Percentage cutoff for filtering isoforms')
    parser.add_argument('--plot-type', type=str, choices=['individual', 'aggregated', 'both'], 
                       default='aggregated', help='Type of plots to generate')
    parser.add_argument('--log-scale', action='store_true', 
                       help='Use logarithmic scale for y-axis')
    parser.add_argument('--stat', choices=['sum', 'mean', 'normalized'], default='sum', 
                       help='Use sum, mean, or normalized for distributions')
    parser.add_argument('--ci-file', type=str, default=None, 
                       help='Path to confidence_intervals.tsv (optional)')
    args = parser.parse_args()

    distributions_out_dir = BASE_DIR / 'data/neuro_project/output/isoform_distributions'

    # Load isoform matrix
    file_path_isoform_matrix = BASE_DIR / 'data/neuro_project/expressed_isoforms_matrix.txt'
    df_isoform_matrix = pd.read_csv(file_path_isoform_matrix, delimiter='\t', index_col=0)

    # Map transcripts to genes
    file_path_isoform_gtf = BASE_DIR / 'data/neuro_project/expressed_isoforms.gtf'
    transcript_id_to_gene = gene_expression_analysis.load_gtf_mapping(file_path_isoform_gtf)
    df_isoform_matrix = gene_expression_analysis.map_transcripts_to_genes(df_isoform_matrix, transcript_id_to_gene)

    # Exclude fetal samples
    all_sample_cols = df_isoform_matrix.columns.difference(['gene_id']).tolist()
    sample_cols = [c for c in all_sample_cols if 'fetal' not in c.lower()]
    df_isoform_matrix = df_isoform_matrix[['gene_id'] + sample_cols]

    # Load CI (optional; used for markers)
    ci_map = load_confidence_intervals(args.ci_file)

    # Generate plots
    if args.plot_type in ['individual', 'both']:
        generate_individual_plots(df_isoform_matrix, sample_cols, distributions_out_dir,
                                 ncols=13, cutoff_pct=args.cutoff_pct,
                                 log_scale=args.log_scale, stat=args.stat, ci_map=ci_map)
    
    if args.plot_type in ['aggregated', 'both']:
        generate_aggregated_plots(df_isoform_matrix, sample_cols, distributions_out_dir,
                                 cutoff_pct=args.cutoff_pct,
                                 log_scale=args.log_scale, stat=args.stat, ci_map=ci_map)


if __name__ == "__main__":
    main()
