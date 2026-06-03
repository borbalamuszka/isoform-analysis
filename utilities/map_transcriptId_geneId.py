"""Transcript/gene mapping and expression stats utility.

This script can:
- Build transcript↔gene mappings from a GTF file and report ambiguous mappings.
- Compute expression and CDS statistics from an expression matrix + GTF.

Usage examples:
    Mapping report (default):
        python -m utilities.map_transcriptId_geneId --gtf /path/to/expressed_isoforms.gtf

    Expression stats:
        python -m utilities.map_transcriptId_geneId --mode stats \
            --gtf /path/to/expressed_isoforms.gtf \
            --matrix /path/to/expressed_isoforms_matrix.txt

Notes:
    Explicit --gtf and --matrix paths are recommended unless you intend to use
    the default neuro_project files.
"""

# %%

from pathlib import Path
import argparse

import pandas as pd

from isoform_distribution.utils import get_filtered_isoforms
from isoform_dashboard.gtf_parser import parse_isoform_file
from .gene_expression_analysis import load_gtf_mapping, map_transcripts_to_genes
from .gtf_utils import parse_gtf_attributes, UNKNOWN_GENE

def load_and_prepare_data(file_path):
    """Load GTF file and prepare data with transcript and gene IDs."""
    data = []
    
    with open(file_path, 'r') as f:
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
                    data.append({
                        'transcript_id': transcript_id,
                        'gene_id': gene_id
                    })
    
    df = pd.DataFrame(data)
    df['gene_id'] = df['gene_id'].fillna(UNKNOWN_GENE).astype(str)
    df['transcript_id'] = df['transcript_id'].astype(str)
    df['transcript_prefix'] = df['transcript_id'].str.split('.', n=1).str[0]
    return df

def print_data_statistics(df):
    print(f"Total rows in data: {len(df)}")
    print(f"Number of unique transcript IDs (e.g. G249.3.nnc): {df['transcript_id'].nunique()}")
    print(f"Number of unique gene IDs (e.g. ENSG00000162426.16): {df['gene_id'].nunique()}")
    print(f"Number of unique transcript prefixes (e.g. G249): {df['transcript_prefix'].nunique()}")
    print(f"Number of unique transcript IDs associated with {UNKNOWN_GENE}: {df.loc[df['gene_id'] == UNKNOWN_GENE, 'transcript_id'].nunique()}")
    print(f"Number of unique transcript prefixes associated with {UNKNOWN_GENE}: {df.loc[df['gene_id'] == UNKNOWN_GENE, 'transcript_prefix'].nunique()}\n")

def build_mappings(df):
    dup_df = df[['gene_id', 'transcript_prefix']].drop_duplicates()
    gene_to_prefixes = dup_df.groupby('gene_id')['transcript_prefix'].agg(lambda s: sorted(s))
    prefix_to_genes = dup_df.groupby('transcript_prefix')['gene_id'].agg(lambda s: sorted(s))
    return gene_to_prefixes, prefix_to_genes

def report_mappings(gene_to_prefixes, prefix_to_genes):
    print("# Genes linked to multiple transcript prefixes")
    for gene in sorted(gene_to_prefixes.index):
        prefixes = gene_to_prefixes.loc[gene]
        if len(prefixes) > 1:
            print(f"{gene}: {', '.join(prefixes)}")

    print("\n# Transcript prefixes associated with multiple genes")
    lengths = prefix_to_genes.str.len()
    multi_gene_prefixes = prefix_to_genes.index[lengths > 1]
    for prefix in multi_gene_prefixes:
        genes = prefix_to_genes.loc[prefix]
        print(f"{prefix}: {', '.join(genes)}")
    multi_gene_count = len(multi_gene_prefixes)

    # vectorized detection of prefixes that include UNKNOWN_GENE
    prefixes_with_unknown_gene = prefix_to_genes.index[
        prefix_to_genes.apply(lambda gs: UNKNOWN_GENE in gs)
    ].tolist()

    print(f"Number of transcript prefixes associated with multiple genes: {multi_gene_count}")

    print("\nTranscript prefixes associated with UNKNOWN_GENE:")
    for prefix in prefixes_with_unknown_gene:
        genes = prefix_to_genes.loc[prefix]
        print(f"{prefix}: {', '.join(genes)}")


def transcript_has_cds(exons):
    """Return True if any exon has a valid CDS interval."""
    return any(
        exon.get("cds_start") is not None
        and exon.get("cds_end") is not None
        and exon["cds_end"] > exon["cds_start"]
        for exon in exons
    )


def _load_expression_matrix(matrix_path, gtf_path, exclude_fetal=True):
    df_raw = pd.read_csv(matrix_path, sep="\t", index_col=0)
    if "gene_id" not in df_raw.columns:
        transcript_to_gene = load_gtf_mapping(gtf_path)
        df_raw = map_transcripts_to_genes(df_raw, transcript_to_gene)

    sample_cols = [
        c for c in df_raw.columns
        if c != "gene_id" and (not exclude_fetal or "fetal" not in c.lower())
    ]
    return df_raw, sample_cols


def _collect_cutoff_transcripts(df_raw, sample_cols, cutoff_pct):
    allowed_transcripts = set()
    for gene_id, gene_block in df_raw.groupby("gene_id"):
        gene_t = gene_block.drop(columns=["gene_id"])
        kept = get_filtered_isoforms(gene_t, sample_cols, cutoff_pct, stat="sum")
        if kept:
            allowed_transcripts.update(kept)
    return allowed_transcripts


def report_expression_stats(gtf_path, matrix_path, cutoff_pct=1.5):
    df_raw_all, _ = _load_expression_matrix(matrix_path, gtf_path, exclude_fetal=False)
    total_genes = df_raw_all["gene_id"].dropna().nunique()
    total_transcripts = df_raw_all.index.nunique()

    expr_sums_all = df_raw_all.drop(columns=["gene_id"]).sum(axis=1)
    expressed_transcripts_all = expr_sums_all[expr_sums_all > 0].index
    expressed_genes_all = df_raw_all.loc[expressed_transcripts_all, "gene_id"].dropna().unique()

    df_raw, sample_cols = _load_expression_matrix(matrix_path, gtf_path)

    expr_sums = df_raw[sample_cols].sum(axis=1)
    expressed_transcripts = expr_sums[expr_sums > 0].index
    expressed_genes = df_raw.loc[expressed_transcripts, "gene_id"].dropna().unique()

    print("# Expression matrix stats")
    print(f"Unique genes (including fetal): {total_genes}")
    print(f"Unique transcripts (including fetal): {total_transcripts}")
    print(f"Unique genes with expression (including fetal): {len(expressed_genes_all)}")
    print(f"Unique transcripts with expression (including fetal): {len(expressed_transcripts_all)}")
    print(f"Unique genes with expression (excluding fetal): {len(expressed_genes)}")
    print(f"Unique transcripts with expression (excluding fetal): {len(expressed_transcripts)}")

    allowed_transcripts = _collect_cutoff_transcripts(df_raw, sample_cols, cutoff_pct)
    print(f"Transcripts contributing >= {cutoff_pct}%: {len(allowed_transcripts)}")

    isoforms_by_gene = parse_isoform_file(gtf_path)
    cds_transcripts = set()
    for gene_id, isoforms in isoforms_by_gene.items():
        for transcript_id, exons in isoforms.items():
            if transcript_id in allowed_transcripts and transcript_has_cds(exons):
                cds_transcripts.add(transcript_id)

    genes_with_cds = set()
    for gene_id, isoforms in isoforms_by_gene.items():
        if any(
            transcript_id in allowed_transcripts and transcript_has_cds(exons)
            for transcript_id, exons in isoforms.items()
        ):
            genes_with_cds.add(gene_id)

    print(f"\n# CDS stats (cutoff-filtered, excluding fetal)")
    print(f"Genes with >=1 CDS transcript: {len(genes_with_cds)}")
    print(f"Transcripts with CDS in those genes: {len(cds_transcripts)}")

    genes_with_2_cds = set()
    cds_transcripts_in_2 = set()
    for gene_id, isoforms in isoforms_by_gene.items():
        cds_ids = [
            transcript_id
            for transcript_id, exons in isoforms.items()
            if transcript_id in allowed_transcripts and transcript_has_cds(exons)
        ]
        if len(cds_ids) >= 2:
            genes_with_2_cds.add(gene_id)
            cds_transcripts_in_2.update(cds_ids)

    print("\n# CDS stats (>=2 CDS isoforms)")
    print(f"Genes with >=2 CDS transcripts: {len(genes_with_2_cds)}")
    print(f"Transcripts with CDS in those genes: {len(cds_transcripts_in_2)}")


def parse_args():
    base_dir = Path(__file__).resolve().parents[1]
    default_gtf = base_dir / "data/neuro_project/expressed_isoforms.gtf"
    default_matrix = base_dir / "data/neuro_project/expressed_isoforms_matrix.txt"
    parser = argparse.ArgumentParser(description="Transcript/gene mapping and stats")
    parser.add_argument(
        "--gtf",
        default=default_gtf,
        type=Path,
        help="GTF file path",
    )
    parser.add_argument(
        "--matrix",
        default=default_matrix,
        type=Path,
        help="Expression matrix TSV path",
    )
    parser.add_argument(
        "--mode",
        choices=["mapping", "stats"],
        default="mapping",
        help="Run mapping output or expression/CDS stats",
    )
    parser.add_argument(
        "--cutoff-pct",
        type=float,
        default=1.5,
        help="Transcript contribution cutoff (percent) for stats mode",
    )
    args = parser.parse_args()
    args._default_gtf = default_gtf
    args._default_matrix = default_matrix
    return args


def main():
    args = parse_args()

    if args.mode == "stats":
        report_expression_stats(
            args.gtf,
            args.matrix,
            cutoff_pct=args.cutoff_pct
        )
        return

    df = load_and_prepare_data(args.gtf)
    print_data_statistics(df)
    gene_to_prefix_series, prefix_to_genes = build_mappings(df)
    report_mappings(gene_to_prefix_series, prefix_to_genes)

if __name__ == "__main__":
    main()

# %%
