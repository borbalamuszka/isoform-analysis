#!/usr/bin/env python3
"""Dash dashboard for isoform entropy exploration.

Features:
- Dataset toggle: Switch between mean and sum datasets in real-time
- Split view: Left side shows visualizations, right side shows gene data table
- Single interactive scatter plot: Summed vs Top Isoform Entropy colored by min Spearman
- Click any point to update the isoform distribution panel:
  * Bar chart showing isoform expression per sample with global mean/sum background
  * Toggle to show/hide confidence intervals (optional, requires --ci-file)
  * Confidence intervals automatically use condition/region grouping from data
  * Exon structure visualization with orange coding exons and blue non-coding exons
- Data table showing gene rankings, correlations, and entropies
  * Click Gene ID to highlight rows (click again to unhighlight)
  * Export highlighted genes to a timestamped text file

Run:
  python dashboard_app.py \
    --input-mean data/distributions_condition_mean.tsv \
    --input-sum data/distributions_condition_sum.tsv \
    --ci-file data/bootstrap/confidence_intervals.tsv \
    --exons data/expressed_isoforms.gtf \
    --geometry-dir data/neuro_project/output/alphafold_geometry

Requires: dash, plotly, pandas, scipy
"""
import argparse
import logging
import os
import sys
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s | %(name)s | %(message)s",
)

from .data_processing import calculate_entropy_and_correlation
from .gtf_parser import parse_isoform_file, parse_gene_names
from .interpro_parser import build_domain_mapping
from .app_layout import create_app


def parse_args(argv=None):
    """Parse command line arguments."""
    p = argparse.ArgumentParser(description="Run Dash dashboard for isoform entropy/correlation exploration")
    p.add_argument("--input-mean", required=False, default=None,
                   help="Input gene/isoform expression TSV (mean values, optional)")
    p.add_argument("--input-sum",
                   help="Input gene/isoform expression TSV (sum values, optional)")
    p.add_argument("--gene-coexpression-dir",
                   help="Path to directory containing precomputed sparse coexpression matrices (gene_coexpression.npz, etc.)",
                   default=None)
    p.add_argument("--ci-file", 
                default=None,
                   help="Bootstrap confidence intervals TSV (optional)")
    p.add_argument("--exons", 
                default=None,
                   help="Input GTF file with exon and CDS annotations (optional)")
    p.add_argument("--geometry-dir",
                   help="Path to alphafold_geometry output directory produced by "
                        "extract_3d_geometry.py.",
                default=None)
    p.add_argument("--proteins",
                   help="FASTA file containing protein sequences (optional).",
                default=None)
    p.add_argument("--interpro-dir",
                   help="Path to InterPro scan results directory (optional).",
                default=None)
    p.add_argument("--interpro-evalue-threshold", type=float, default=1e-5,
                   help="E-value threshold for InterPro domain significance (default: 1e-5)")
    p.add_argument("--port", type=int, default=8050, help="Port to serve the dashboard")
    p.add_argument("--host", default=None, help="Host interface (use 0.0.0.0 for LAN access)")
    p.add_argument("--default-ranking", choices=["spearman", "expression"], default="spearman",
                   help="Default sort order for the gene table: "
                        "'spearman' (Spearman/entropy-based, default) or "
                        "'expression' (highest mean expression first)")
    return p.parse_args(argv)


def load_protein_sequences(fasta_file: str) -> dict:
    """Load protein sequences from a FASTA file.

    Args:
        fasta_file: Path to FASTA file containing protein sequences.

    Returns:
        Dictionary mapping transcript_id -> amino-acid sequence string.
    """
    protein_sequences: dict = {}

    if not fasta_file or not os.path.exists(fasta_file):
        if fasta_file:
            print(f"Warning: Protein FASTA file not found: {fasta_file}", file=sys.stderr)
        return protein_sequences

    try:
        current_id: str | None = None
        current_seq: list[str] = []

        with open(fasta_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id and current_seq:
                        protein_sequences[current_id] = "".join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_id and current_seq:
            protein_sequences[current_id] = "".join(current_seq)

        if protein_sequences:
            print(f"Loaded {len(protein_sequences)} protein sequences from {fasta_file}")
        else:
            print(f"No protein sequences found in {fasta_file}", file=sys.stderr)
    except Exception as e:
        print(f"Error loading protein sequences: {e}", file=sys.stderr)

    return protein_sequences


def main():
    """Main entry point for the dashboard application."""
    args = parse_args()

    df_mean = None
    results_df_mean = pd.DataFrame()
    sample_cols = []
    global_col_mean = ""

    # Load mean dataset if provided
    if args.input_mean:
        if not os.path.exists(args.input_mean):
            print(f"Input file not found: {args.input_mean}", file=sys.stderr)
            sys.exit(1)

        try:
            df_mean = pd.read_csv(args.input_mean, sep="\t")
        except Exception as e:
            print(f"Failed to read input TSV (mean): {e}", file=sys.stderr)
            sys.exit(1)

        if "gene_id" not in df_mean.columns or "transcript_id" not in df_mean.columns:
            print("Input must contain gene_id and transcript_id columns", file=sys.stderr)
            sys.exit(1)

        if len(df_mean.columns) < 3:
            print("Input must have at least 3 columns (gene_id, transcript_id, and a global column)", file=sys.stderr)
            sys.exit(1)
        
        global_col_mean = df_mean.columns[2]
        print(f"Using '{global_col_mean}' as the global aggregation column for mean dataset")

        # Identify sample columns
        meta_cols_mean = {"gene_id", "transcript_id", global_col_mean}
        sample_cols = [c for c in df_mean.columns if c not in meta_cols_mean and pd.api.types.is_numeric_dtype(df_mean[c])]
        
        if len(sample_cols) < 2:
            print("Need at least two numeric sample columns to compute correlations", file=sys.stderr)
            sys.exit(1)

        # Calculate entropy and correlations for mean dataset
        results_mean = calculate_entropy_and_correlation(df_mean, sample_cols, global_col_mean)
        results_df_mean = pd.DataFrame(results_mean)

    # Load sum dataset if provided
    has_sum = False
    df_sum = None
    results_df_sum = pd.DataFrame()
    global_col_sum = ""

    if args.input_sum:
        if not os.path.exists(args.input_sum):
            print(f"Warning: Sum input file not found: {args.input_sum}", file=sys.stderr)
        else:
            try:
                df_sum = pd.read_csv(args.input_sum, sep="\t")
                if "gene_id" not in df_sum.columns or "transcript_id" not in df_sum.columns:
                    print("Sum input must contain gene_id and transcript_id columns", file=sys.stderr)
                else:
                    has_sum = True
                    print(f"Loaded sum dataset from: {args.input_sum}")
                    global_col_sum = df_sum.columns[2]
                    print(f"Using '{global_col_sum}' as the global aggregation column for sum dataset")
                    
                    if not sample_cols:
                        meta_cols_sum = {"gene_id", "transcript_id", global_col_sum}
                        sample_cols = [c for c in df_sum.columns if c not in meta_cols_sum and pd.api.types.is_numeric_dtype(df_sum[c])]
                        if len(sample_cols) < 2:
                            print("Need at least two numeric sample columns to compute correlations", file=sys.stderr)
                            has_sum = False
                            df_sum = None
                            
                    if has_sum:
                        results_sum = calculate_entropy_and_correlation(df_sum, sample_cols, global_col_sum)
                        results_df_sum = pd.DataFrame(results_sum)
            except Exception as e:
                print(f"Failed to read sum TSV: {e}", file=sys.stderr)

    if not has_sum:
        # Create empty placeholder
        global_col_sum = global_col_mean
        results_df_sum = pd.DataFrame()

    # Verify that at least one primary widget source is supplied
    if not (args.input_mean or args.input_sum or args.exons or args.gene_coexpression_dir):
        print("Warning: Insufficient data to render any widget. Please supply at least one primary source: --input-mean, --input-sum, --exons, or --gene-coexpression-dir.", file=sys.stderr)

    # Load exon data
    isoforms_by_gene = {}
    gene_names = {}
    if args.exons:
        if os.path.exists(args.exons):
            print(f"Loading exon structures from GTF file: {args.exons}")
            try:
                isoforms_by_gene = parse_isoform_file(args.exons)
                print(f"Loaded exon data for {len(isoforms_by_gene)} genes")
                
                # Also load gene names
                print(f"Loading gene names from GTF file: {args.exons}")
                gene_names = parse_gene_names(args.exons)
                print(f"Loaded gene names for {len(gene_names)} genes")
            except Exception as e:
                print(f"Failed to parse GTF: {e}", file=sys.stderr)
        else:
            print(f"Warning: GTF file not found: {args.exons}", file=sys.stderr)

    # Load confidence intervals from bootstrap output
    ci_df = None
    ci_columns = []
    if args.ci_file:
        if os.path.exists(args.ci_file):
            print(f"Loading confidence intervals from: {args.ci_file}")
            try:
                ci_df = pd.read_csv(args.ci_file, sep="\t")
                if 'isoform' in ci_df.columns:
                    ci_df = ci_df.set_index('isoform')
                elif ci_df.index.name != 'isoform':
                    ci_df.index.name = 'isoform'
                ci_columns = [col for col in ci_df.columns if col.startswith('ci_')]
                print(f"Loaded CI data with {len(ci_columns)} CI columns for {len(ci_df)} isoforms")
            except Exception as e:
                print(f"Failed to read CI file: {e}", file=sys.stderr)
        else:
            print(f"Warning: CI file not found: {args.ci_file}", file=sys.stderr)

    # Load precomputed coexpression matrices
    gene_coexpression = None
    gene_coexpression_idx = None
    isoform_coexpression = None
    isoform_coexpression_idx = None
    gene_iso_mapping = {}
    
    if args.gene_coexpression_dir:
        import pickle
        from scipy.sparse import load_npz
        print(f"Loading coexpression matrices from: {args.gene_coexpression_dir}")
        try:
            gene_mat_path = os.path.join(args.gene_coexpression_dir, "gene_coexpression.npz")
            gene_idx_path = os.path.join(args.gene_coexpression_dir, "gene_index.pkl")
            if os.path.exists(gene_mat_path) and os.path.exists(gene_idx_path):
                gene_coexpression = load_npz(gene_mat_path)
                with open(gene_idx_path, 'rb') as f:
                    gene_coexpression_idx = pickle.load(f)
                    
            iso_mat_path = os.path.join(args.gene_coexpression_dir, "isoform_coexpression.npz")
            iso_idx_path = os.path.join(args.gene_coexpression_dir, "isoform_index.pkl")
            if os.path.exists(iso_mat_path) and os.path.exists(iso_idx_path):
                isoform_coexpression = load_npz(iso_mat_path)
                with open(iso_idx_path, 'rb') as f:
                    isoform_coexpression_idx = pickle.load(f)
                    
            # Load precomputed gene-isoform mapping table
            mapping_path = os.path.join(args.gene_coexpression_dir, "gene_iso_mapping.pkl")
            if os.path.exists(mapping_path):
                print(f"Loading gene-isoform mapping table from: {mapping_path}")
                with open(mapping_path, 'rb') as f:
                    gene_iso_mapping = pickle.load(f)
                    
        except Exception as e:
            print(f"Failed to read coexpression matrices: {e}", file=sys.stderr)

    # Create the Dash app
    protein_sequences = load_protein_sequences(args.proteins)
    
    # Load InterPro domain data
    domain_mapping = {}
    if args.interpro_dir:
        if os.path.isdir(args.interpro_dir):
            print(f"Loading InterPro domain data from: {args.interpro_dir}")
            try:
                domain_mapping = build_domain_mapping(
                    args.interpro_dir,
                    isoforms_by_gene,
                    min_evalue=args.interpro_evalue_threshold
                )
                print(f"Loaded InterPro data for {len(domain_mapping)} transcripts")
            except Exception as e:
                print(f"Warning: Failed to load InterPro data: {e}", file=sys.stderr)
        else:
            print(f"Warning: InterPro directory not found: {args.interpro_dir}", file=sys.stderr)

    app = create_app(df_mean, df_sum if has_sum else df_mean,
                    results_df_mean, results_df_sum,
                    sample_cols, ci_df, ci_columns,
                    global_col_mean, global_col_sum,
                    isoforms_by_gene, gene_names, has_sum,
                    geometry_dir=args.geometry_dir,
                    protein_sequences=protein_sequences,
                    domain_mapping=domain_mapping,
                    default_ranking=args.default_ranking,
                    gene_coexpression=gene_coexpression,
                    gene_coexpression_idx=gene_coexpression_idx,
                    isoform_coexpression=isoform_coexpression,
                    isoform_coexpression_idx=isoform_coexpression_idx,
                    gene_iso_mapping=gene_iso_mapping,
                    path_mean=args.input_mean,
                    path_sum=args.input_sum,
                    path_gtf=args.exons,
                    path_proteins=args.proteins,
                    path_interpro_dir=args.interpro_dir,
                    path_coexpression_dir=args.gene_coexpression_dir,
                    path_ci=args.ci_file)
    
    # Get port from environment variable (for Render) or args or default
    port = args.port if args.port is not None else int(os.environ.get('PORT', 8050))
    
    # Get host from args or environment (use 0.0.0.0 for cloud deployment, 127.0.0.1 for local)
    host = args.host if args.host is not None else os.environ.get('HOST', '127.0.0.1')
    
    print(f"* Serving dashboard on http://{host}:{port}")
    app.run(debug=False, host=host, port=port)


if __name__ == "__main__":
    main()
