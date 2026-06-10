"""
Isoform distribution table generation script.
Generates TSV tables of isoform expression data.
For plotting, use plots.py instead.
"""

import pandas as pd
import utilities.gene_expression_analysis as gene_expression_analysis
import argparse
from pathlib import Path
from collections import defaultdict
from .utils import (
    get_filtered_isoforms,
    prepare_gene_list_and_paths,
    _normalize_gene_samples
)

# Define base directory relative to this script
BASE_DIR = Path(__file__).parent.parent  # source/


def write_isoform_table(df_isoform_matrix, sample_cols, all_genes, output_dir, suffix, 
                       cutoff_pct=2, stat='sum'):
    """
    Write a single TSV with one row per retained transcript across all genes.
    
    Args:
        df_isoform_matrix: DataFrame with isoform expression data
        sample_cols: List of sample column names
        all_genes: List of gene IDs to include
        output_dir: Output directory path
        suffix: Suffix for the output filename (e.g., 'individual', 'region', 'condition')
        cutoff_pct: Minimum percentage contribution to keep an isoform
        stat: 'sum', 'mean', or 'normalized'
        
    Output columns:
        gene_id, transcript_id, global_sum or global_mean, 
        <sample/group columns...>
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for gene_id in all_genes:
        gene_transcripts = df_isoform_matrix[df_isoform_matrix['gene_id'] == gene_id].drop(columns=['gene_id'])

        # Compute normalized matrix once per gene when needed
        norm = None
        if stat == 'normalized':
            norm = _normalize_gene_samples(gene_transcripts, sample_cols)

        filtered_isoforms = get_filtered_isoforms(gene_transcripts, sample_cols, cutoff_pct, stat)
        if filtered_isoforms is None:
            continue

        if stat == 'mean':
            global_expr = gene_transcripts[sample_cols].mean(axis=1).loc[filtered_isoforms]
        elif stat == 'normalized':
            global_expr = norm.mean(axis=1).loc[filtered_isoforms]
        else:
            global_expr = gene_transcripts[sample_cols].sum(axis=1).loc[filtered_isoforms]

        sorted_isoforms = [isoform for _, isoform in sorted(zip(global_expr, filtered_isoforms), reverse=True)]
        global_expr_sorted = global_expr.loc[sorted_isoforms]

        for isoform in sorted_isoforms:
            row = {
                'gene_id': gene_id,
                'transcript_id': isoform,
                ('global_mean' if stat in ['mean', 'normalized'] else 'global_sum'): float(global_expr_sorted.loc[isoform]),
            }
            if stat == 'normalized':
                for sc in sample_cols:
                    row[sc] = float(norm.at[isoform, sc])
            else:
                for sc in sample_cols:
                    row[sc] = float(gene_transcripts.at[isoform, sc])
            rows.append(row)

    output_df = pd.DataFrame(rows)
    out_path = output_dir / f"distributions_{suffix}_{stat}.tsv"
    output_df.to_csv(out_path, sep='\t', index=False)
    print(f"Wrote {out_path}")


def generate_individual_tables(df_isoform_matrix, sample_cols, out_dir, cutoff_pct=2, stat='sum'):
    """
    Generate tables for individual samples.
    
    Args:
        df_isoform_matrix: DataFrame with isoform expression data
        sample_cols: List of sample column names
        out_dir: Output directory path
        cutoff_pct: Minimum percentage contribution to keep an isoform
        stat: 'sum', 'mean', or 'normalized'
    """
    out_dir, all_genes = prepare_gene_list_and_paths(df_isoform_matrix, out_dir)
    write_isoform_table(df_isoform_matrix, sample_cols, all_genes, out_dir, "individual", 
                       cutoff_pct=cutoff_pct, stat=stat)


def _load_metadata(meta_path, sample_col, group_col):
    meta_df = pd.read_csv(meta_path, sep=r'\s+', quotechar='"', engine='python', dtype=str)
    meta_df.columns = [c.strip().strip('"') for c in meta_df.columns]
    normalized = {c.lower(): c for c in meta_df.columns}
    required = [sample_col.lower(), group_col.lower()]
    if not all(col in normalized for col in required):
        raise ValueError(
            f"Metadata file missing required columns: {required}. Found: {meta_df.columns.tolist()}"
        )
    meta_df = meta_df.rename(columns={
        normalized[sample_col.lower()]: 'sample_id',
        normalized[group_col.lower()]: 'group'
    })
    for col in ['sample_id', 'group']:
        meta_df[col] = meta_df[col].astype(str).str.strip().str.strip('"')
    return meta_df


def _build_groups_from_metadata(meta_df, sample_cols):
    cleaned = meta_df[['sample_id', 'group']].dropna().copy()

    sample_id_map = cleaned.set_index('sample_id')['group'].to_dict()

    grouped = defaultdict(list)
    for sample_col in sample_cols:
        group = sample_id_map.get(sample_col)

        if group:
            grouped[group].append(sample_col)
    return grouped


def _aggregate_samples_by_mapping(df_isoform_matrix, grouped_cols, stat='sum'):
    result_data = {}
    for group_key, group_sample_cols in grouped_cols.items():
        if stat == 'mean':
            result_data[group_key] = df_isoform_matrix[group_sample_cols].mean(axis=1)
        elif stat == 'normalized':
            vals = []
            for _, gene_block in df_isoform_matrix.groupby('gene_id'):
                gb = gene_block.drop(columns=['gene_id'])
                norm = _normalize_gene_samples(gb, group_sample_cols).mean(axis=1)
                vals.append(norm)
            result_data[group_key] = pd.concat(vals).reindex(df_isoform_matrix.index)
        else:
            result_data[group_key] = df_isoform_matrix[group_sample_cols].sum(axis=1)

    result_df = pd.DataFrame(result_data)
    result_df['gene_id'] = df_isoform_matrix['gene_id']
    return result_df


def generate_metadata_group_tables(df_isoform_matrix, sample_cols, out_dir, meta_path,
                                   meta_sample_col='sample_id', meta_group_col='cell_type',
                                   group_label='group',
                                   cutoff_pct=2, stat='sum'):
    """
    Generate aggregated tables using metadata-provided grouping.

    Args:
        df_isoform_matrix: DataFrame with isoform expression data
        sample_cols: List of sample column names
        out_dir: Output directory path
        meta_path: Path to metadata file
        meta_sample_col: Column in metadata that identifies samples
        meta_group_col: Column in metadata to group by
        group_label: Label used in output filename
        cutoff_pct: Minimum percentage contribution to keep an isoform
        stat: 'sum', 'mean', or 'normalized'
    """
    out_dir, all_genes = prepare_gene_list_and_paths(df_isoform_matrix, out_dir)
    meta_df = _load_metadata(meta_path, meta_sample_col, meta_group_col)
    grouped_cols = _build_groups_from_metadata(
        meta_df,
        sample_cols
    )
    if not grouped_cols:
        raise ValueError("No groups found from metadata; check sample IDs and column names.")

    print(f"Generating {group_label} table with stat={stat}...")
    df_aggregated = _aggregate_samples_by_mapping(df_isoform_matrix, grouped_cols, stat=stat)
    agg_cols = [col for col in df_aggregated.columns if col != 'gene_id']
    write_isoform_table(df_aggregated, agg_cols, all_genes, out_dir, group_label,
                        cutoff_pct=cutoff_pct, stat=stat)


def main():
    """Main entry point for table generation script."""
    parser = argparse.ArgumentParser(
        description='Generate isoform distribution tables (TSV format only). '
                   'For plotting, use plots.py instead.'
    )
    parser.add_argument('--cutoff-pct', type=float, default=1.5, help='Percentage cutoff for filtering isoforms')
    parser.add_argument('--table-type', type=str, choices=['individual', 'aggregated', 'both'], default='aggregated', help='Type of tables to generate')
    parser.add_argument('--stat', choices=['sum', 'mean', 'normalized'], default='mean', help='Use sum, mean, or normalized for distributions')
    parser.add_argument('--matrix', type=str, required=True, help='Path to isoform expression matrix')
    parser.add_argument('--gtf', type=str, required=True, help='Path to GTF file for transcript->gene mapping')
    parser.add_argument('--meta-file', type=str, required=True, help='Path to metadata file for grouping')
    parser.add_argument('--output-dir', type=str, required=True, help='Output directory for distributions')
    parser.add_argument('--meta-sample-col', type=str, default='sample_id',
                        help='Metadata column that identifies samples')
    parser.add_argument('--meta-group-col', type=str, default=None,
                        help='Metadata column to group by')
    parser.add_argument('--exclude-sample-substr', action='append', default=[],
                        help='Exclude samples containing this substring (can be repeated)')
    args = parser.parse_args()

    file_path_isoform_matrix = Path(args.matrix)
    file_path_isoform_gtf = Path(args.gtf)
    meta_path = Path(args.meta_file)
    distributions_out_dir = Path(args.output_dir)

    df_isoform_matrix = pd.read_csv(file_path_isoform_matrix, delimiter='\t', index_col=0)
    transcript_id_to_gene = gene_expression_analysis.load_gtf_mapping(file_path_isoform_gtf)
    df_isoform_matrix = gene_expression_analysis.map_transcripts_to_genes(df_isoform_matrix, transcript_id_to_gene)

    all_sample_cols = df_isoform_matrix.columns.difference(['gene_id']).tolist()
    df_isoform_matrix[all_sample_cols] = (
        df_isoform_matrix[all_sample_cols].apply(pd.to_numeric, errors='coerce').fillna(0.0)
    )

    sample_cols = all_sample_cols
    if args.exclude_sample_substr:
        exclude_tokens = [tok.lower() for tok in args.exclude_sample_substr]
        sample_cols = [
            c for c in all_sample_cols
            if not any(tok in c.lower() for tok in exclude_tokens)
        ]
        df_isoform_matrix = df_isoform_matrix[['gene_id'] + sample_cols]

    if args.table_type in ['individual', 'both']:
        generate_individual_tables(df_isoform_matrix, sample_cols, distributions_out_dir,
                                   cutoff_pct=args.cutoff_pct, stat=args.stat)

    if args.table_type in ['aggregated', 'both']:
        generate_metadata_group_tables(
            df_isoform_matrix,
            sample_cols,
            distributions_out_dir,
            meta_path=meta_path,
            meta_sample_col=args.meta_sample_col,
            meta_group_col=args.meta_group_col,
            group_label='group',
            cutoff_pct=args.cutoff_pct,
            stat=args.stat
        )


if __name__ == "__main__":
    main()

