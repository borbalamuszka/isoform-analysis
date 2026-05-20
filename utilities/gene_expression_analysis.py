import pandas as pd
from .gtf_utils import parse_gtf_attributes, UNKNOWN_GENE


def load_and_prepare_data(file_path):
    df = pd.read_csv(file_path, delimiter='\t')
    df['gene_id'] = df['gene_id'].fillna(UNKNOWN_GENE).astype(str)
    df['transcript_id'] = df['transcript_id'].astype(str)
    return df


def load_gtf_mapping(gtf_path):
    """Load GTF file and extract transcript_id to gene_id mapping."""
    mapping = {}
    
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            # Only process transcript entries to get the mapping
            if feature_type == 'transcript':
                attributes = fields[8]
                transcript_id, gene_id = parse_gtf_attributes(attributes)
                
                if transcript_id:
                    mapping[transcript_id] = gene_id
    
    return pd.Series(mapping)

def build_mapping(df):
    dup_df = df[['gene_id', 'transcript_id']].drop_duplicates()
    transcript_id_to_gene = dup_df.set_index('transcript_id')['gene_id']
    return transcript_id_to_gene

def map_transcripts_to_genes(df_isoform_matrix, transcript_id_to_gene):
    """
    Add a 'gene_id' column to the isoform matrix by mapping transcript IDs to gene IDs.
    Returns modified DataFrame.
    """
    df_isoform_matrix['gene_id'] = df_isoform_matrix.index.to_series().map(transcript_id_to_gene)
    return df_isoform_matrix

def compute_gene_expression_counts(df_isoform_matrix, sample_cols):
    """
    Count transcripts per gene per sample where expression > 0.
    Returns DataFrame indexed by gene_id with sample columns.
    """
    bool_df = df_isoform_matrix[sample_cols] > 0
    bool_df['gene_id'] = df_isoform_matrix['gene_id']
    df_gene_counts = bool_df.groupby('gene_id').sum()
    return df_gene_counts

def count_transcripts_for_90_percent(group, sample_col):
    """
    For a group of transcripts (same gene), count how many transcripts
    (sorted by expression descending) are needed to reach 90% of total expression.
    """
    values = group[sample_col]
    sorted_values = values.sort_values(ascending=False)
    total = sorted_values.sum()
    if total == 0:
        return 0
    cumsum = sorted_values.cumsum()
    threshold = 0.9 * total
    count = (cumsum < threshold).sum() + 1
    return min(count, len(sorted_values))

def compute_gene_90percent_counts(df_isoform_matrix, df_gene_counts, sample_cols):
    """
    For each gene and sample, calculate how many transcripts contribute to 90%
    of total expression. Returns DataFrame with same structure as df_gene_counts.
    """
    df_gene_90percent = pd.DataFrame(0, index=df_gene_counts.index, columns=sample_cols)
    
    for gene_id in df_gene_90percent.index:
        gene_transcripts = df_isoform_matrix[df_isoform_matrix['gene_id'] == gene_id]
        if len(gene_transcripts) == 0:
            continue
        for sample_col in sample_cols:
            count = count_transcripts_for_90_percent(gene_transcripts, sample_col)
            df_gene_90percent.loc[gene_id, sample_col] = count
    
    return df_gene_90percent

def main():
    # File paths
    file_path_isoform_gtf = './data/neuro_project/expressed_isoforms.gtf'
    file_path_isoform_matrix = './data/neuro_project/expressed_isoforms_matrix.txt'
    output_path_counts = './data/neuro_project/output/result_tables/gene_isoform_counts.txt'
    output_path_combined = './data/neuro_project/output/result_tables/gene_isoform_combined.txt'
    
    # Load and prepare data
    transcript_id_to_gene = load_gtf_mapping(file_path_isoform_gtf)
    
    df_isoform_matrix = pd.read_csv(file_path_isoform_matrix, delimiter='\t', index_col=0)
    df_isoform_matrix = map_transcripts_to_genes(df_isoform_matrix, transcript_id_to_gene)
    
    # Get sample columns
    sample_cols = df_isoform_matrix.columns.difference(['gene_id']).tolist()
    
    # Compute gene-level expression counts
    df_gene_counts = compute_gene_expression_counts(df_isoform_matrix, sample_cols)
    
    # Create output directory if it doesn't exist
    from pathlib import Path
    Path(output_path_counts).parent.mkdir(parents=True, exist_ok=True)
    
    df_gene_counts.to_csv(output_path_counts, sep='\t')
    print(f"Saved gene-level isoform counts to {output_path_counts}")
    
    # Compute 90% transcript counts (for combined table)
    df_gene_90percent = compute_gene_90percent_counts(df_isoform_matrix, df_gene_counts, sample_cols)
    
    # Create combined table with two rows per gene (count, then 90%)
    df_gene_counts_labeled = df_gene_counts.copy()
    df_gene_counts_labeled['metric'] = 'count'
    
    df_gene_90percent_labeled = df_gene_90percent.copy()
    df_gene_90percent_labeled['metric'] = '90pct'
    
    # Concatenate vertically and set multi-index (gene_id, metric)
    df_combined = pd.concat([df_gene_counts_labeled, df_gene_90percent_labeled])
    df_combined = df_combined.set_index('metric', append=True)
    
    # Sort so count comes before 90pct for each gene
    df_combined = df_combined.sort_index(level=[0, 1])
    
    df_combined.to_csv(output_path_combined, sep='\t')
    print(f"Saved combined table to {output_path_combined}")
    
    return df_gene_counts, df_combined

if __name__ == "__main__":
    main()
