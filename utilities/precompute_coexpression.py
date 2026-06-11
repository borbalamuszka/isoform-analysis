import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.sparse import csr_matrix, save_npz
import pickle

# Ensure we can import from utilities
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utilities.gtf_utils import parse_gtf_attributes

def parse_gtf_mapping_local(gtf_path):
    """Parse GTF to build transcript_id -> gene_id mapping."""
    mapping = {}
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            transcript_id, gene_id = parse_gtf_attributes(parts[8])
            if transcript_id and gene_id:
                mapping[transcript_id] = gene_id
    return mapping

def compute_correlation_chunked(matrix, threshold=1e-3, chunk_size=250):
    """
    Computes the pairwise correlation matrix in chunks to display a progress bar.
    Takes a 2D numpy array where rows are variables (genes/isoforms) and cols are observations (samples).
    Returns a scipy sparse CSR matrix.
    """
    n_vars = matrix.shape[0]
    
    # Pre-rank the data for Spearman correlation
    print("Ranking data for Spearman correlation...")
    ranked = rankdata(matrix, axis=1)
    
    # Normalize ranked data for dot product correlation
    ranked = ranked - ranked.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(ranked, axis=1, keepdims=True)
    norms[norms == 0] = 1e-10  # Prevent division by zero
    ranked = (ranked / norms).astype(np.float32)
    
    # We will build lists of rows, cols, and data for COO matrix
    data = []
    rows = []
    cols = []
    
    total_chunks = (n_vars + chunk_size - 1) // chunk_size
    
    print("Computing correlation matrix in chunks...")
    for i in range(total_chunks):
        start_i = i * chunk_size
        end_i = min((i + 1) * chunk_size, n_vars)
        
        # Display simple progress bar
        progress = (i / total_chunks) * 100
        sys.stdout.write(f"\rProgress: [{int(progress):3d}%] {'=' * int(progress // 2)}{' ' * (50 - int(progress // 2))}")
        sys.stdout.flush()
        
        chunk = ranked[start_i:end_i]
        
        # Compute correlation of this chunk with the entire matrix
        corr_chunk = np.dot(chunk, ranked.T)
        
        # Ensure self-correlation is 0
        for j in range(end_i - start_i):
            corr_chunk[j, start_i + j] = 0.0
            
        # Apply threshold
        valid_mask = np.abs(corr_chunk) >= threshold
        
        # Get coordinates of valid elements
        r, c = np.where(valid_mask)
        
        # Downcast indices and values to save memory
        r = r.astype(np.int32)
        c = c.astype(np.int32)
        vals = corr_chunk[r, c].astype(np.float32)
        
        # Add to COO lists
        rows.append(r + start_i)
        cols.append(c)
        data.append(vals)
        
    sys.stdout.write(f"\rProgress: [100%] {'=' * 50}\n")
    sys.stdout.flush()
    
    if len(data) == 0:
        return csr_matrix((n_vars, n_vars), dtype=np.float32)
        
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    data = np.concatenate(data)
    
    print(f"Creating sparse matrix with {len(data)} non-zero elements...")
    return csr_matrix((data, (rows, cols)), shape=(n_vars, n_vars), dtype=np.float32)

def main():
    parser = argparse.ArgumentParser(description="Precompute co-expression matrices")
    parser.add_argument("--expression", required=True, help="Input expression matrix TSV")
    parser.add_argument("--gtf", required=True, help="Input GTF file for transcript to gene mapping")
    parser.add_argument("--outdir", default="data/co-expression", help="Output directory for sparse matrices")
    parser.add_argument("--threshold", type=float, default=1e-3, help="Correlation magnitude threshold (default: 1e-3)")
    parser.add_argument("--chunk-size", type=int, default=250, help="Chunk size for correlation computation")
    
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    print(f"Loading expression matrix from {args.expression}...")
    df = pd.read_csv(args.expression, sep='\t')
    
    if 'transcript_id' not in df.columns:
        print("Error: Expression matrix must contain a 'transcript_id' column.")
        sys.exit(1)
        
    print(f"Loading GTF mapping from {args.gtf}...")
    gene_map = parse_gtf_mapping_local(args.gtf)
    
    # Map transcript to gene
    df['gene_id'] = df['transcript_id'].map(gene_map)
    
    unmapped = df['gene_id'].isna().sum()
    if unmapped > 0:
        print(f"Warning: {unmapped} transcripts could not be mapped to a gene_id. They will be ignored.")
        df = df.dropna(subset=['gene_id'])
        
    sample_cols = [c for c in df.columns if c not in ('transcript_id', 'gene_id')]
    
    print(f"Found {len(sample_cols)} sample columns and {len(df)} mapped transcripts.")
    
    # 1. Isoform level
    print("\n--- Isoform-level Co-expression ---")
    isoform_ids = df['transcript_id'].values
    isoform_matrix = df[sample_cols].values
    
    iso_sparse = compute_correlation_chunked(isoform_matrix, threshold=args.threshold, chunk_size=args.chunk_size)
    
    iso_out_mat = os.path.join(args.outdir, "isoform_coexpression.npz")
    iso_out_idx = os.path.join(args.outdir, "isoform_index.pkl")
    save_npz(iso_out_mat, iso_sparse)
    with open(iso_out_idx, 'wb') as f:
        pickle.dump(isoform_ids, f)
    print(f"Saved isoform sparse matrix to {iso_out_mat}")
    
    # 2. Gene level
    print("\n--- Gene-level Co-expression ---")
    print("Aggregating isoform expression to gene level (sum)...")
    gene_df = df.groupby('gene_id')[sample_cols].sum()
    gene_ids = gene_df.index.values
    gene_matrix = gene_df.values
    
    gene_sparse = compute_correlation_chunked(gene_matrix, threshold=args.threshold, chunk_size=args.chunk_size)
    
    gene_out_mat = os.path.join(args.outdir, "gene_coexpression.npz")
    gene_out_idx = os.path.join(args.outdir, "gene_index.pkl")
    save_npz(gene_out_mat, gene_sparse)
    with open(gene_out_idx, 'wb') as f:
        pickle.dump(gene_ids, f)
    print(f"Saved gene sparse matrix to {gene_out_mat}")
    
    print("\nPrecomputation complete.")

if __name__ == "__main__":
    main()