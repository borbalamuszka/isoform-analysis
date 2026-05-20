#!/usr/bin/env python3
"""Generate AlphaFold Server job JSON files from protein sequences.

For each gene (ordered by ranking from data_processing.py), this script:
- Identifies transcripts that have at least one CDS exon
- Deduplicates transcripts with identical protein sequences
- Creates job entries with name = transcriptid_geneid (dots removed)
- Splits jobs into batches of 30 per JSON file

Usage:
    python create_alphafold_jobs.py \
        --input data/neuro_project/output/isoform_distributions/tables/distributions_condition_sum.tsv \
        --gtf data/neuro_project/expressed_isoforms.gtf \
        --fasta data/neuro_project/proteins.fasta \
        --output data/jobs

Outputs:
    data/jobs/alphafold_jobs_1.json
    data/jobs/alphafold_jobs_2.json
    ...
"""
import argparse
import json
import os
import sys

import pandas as pd

from isoform_dashboard.data_processing import calculate_entropy_and_correlation, compute_gene_ranking
from isoform_dashboard.gtf_parser import parse_isoform_file
from isoform_dashboard.dashboard_app import load_protein_sequences
from .gene_expression_analysis import load_gtf_mapping, map_transcripts_to_genes

# ---------------------------------------------------------------------------
# Gene ordering
# ---------------------------------------------------------------------------

def ranked_gene_order(df):
    """Return gene_ids sorted by the ranking logic in data_processing.py.

    Args:
        df: DataFrame with gene_id, transcript_id, and sample columns
            (same format as the distributions TSV files)

    Returns:
        List of gene_id strings in rank order (rank 1 first)
    """
    global_col = df.columns[2]
    meta = {"gene_id", "transcript_id", global_col}
    sample_cols = [
        c for c in df.columns
        if c not in meta and pd.api.types.is_numeric_dtype(df[c])
    ]
    if len(sample_cols) < 2:
        # Not enough samples for correlations — fall back to input order
        return list(df["gene_id"].unique())

    results = calculate_entropy_and_correlation(df, sample_cols, global_col)
    results_df = pd.DataFrame(results)

    ranks = compute_gene_ranking(results_df)
    results_df["rank"] = ranks
    results_df_sorted = results_df.sort_values("rank")
    return results_df_sorted["gene_id"].tolist()


# ---------------------------------------------------------------------------
# CDS helpers  (mirrors data_processing.gene_has_cds at transcript level)
# ---------------------------------------------------------------------------

def transcript_has_cds(exons):
    """Return True if any exon in the list contains a CDS."""
    return any(
        e.get("cds_start") is not None
        and e.get("cds_end") is not None
        and e["cds_end"] > e["cds_start"]
        for e in exons
    )


# ---------------------------------------------------------------------------
# Job building
# ---------------------------------------------------------------------------

def make_job_name(transcript_id, gene_id):
    """Return job name: transcriptid_geneid with all dots removed."""
    return f"{transcript_id}_{gene_id}".replace(".", "")


def collect_submitted_names(jobs_dir: str, max_index: int) -> set:
    """Return the set of job names already submitted (files 1..max_index).

    Handles both plain ``alphafold_jobs_N.json`` and date-suffixed
    ``alphafold_jobs_N-doneMMDD.json`` variants.

    Args:
        jobs_dir: Directory containing the job JSON files.
        max_index: Highest batch index that was submitted (inclusive).

    Returns:
        Set of job name strings (e.g. ``'G18635.44nnc_ENSG00000102781'``).
    """
    submitted: set = set()
    if not os.path.isdir(jobs_dir):
        print(f"Warning: jobs directory not found: {jobs_dir}", file=sys.stderr)
        return submitted
    for fname in os.listdir(jobs_dir):
        if not fname.endswith(".json") or not fname.startswith("alphafold_jobs_"):
            continue
        stem = fname[len("alphafold_jobs_"):-len(".json")]
        num_str = stem.split("-")[0]
        try:
            num = int(num_str)
        except ValueError:
            continue
        if num > max_index:
            continue
        path = os.path.join(jobs_dir, fname)
        try:
            with open(path) as fh:
                for job in json.load(fh):
                    submitted.add(job["name"])
        except Exception as e:
            print(f"Warning: could not read {fname}: {e}", file=sys.stderr)
    return submitted


def build_jobs(gene_order, isoforms_by_gene, protein_sequences,
               excluded_names: set = None,
               allowed_transcripts: set = None):
    """Build a flat list of AlphaFold job dicts ordered by gene ranking.

    For each gene:
      - Keep only transcripts that have at least one CDS exon
      - Optionally restrict to ``allowed_transcripts`` (cutoff-filtered set)
      - Deduplicate by protein sequence (keep first occurrence per gene)
      - Skip jobs whose name is in ``excluded_names`` (already submitted)

    Args:
        gene_order: Gene IDs in priority order.
        isoforms_by_gene: GTF exon structures.
        protein_sequences: Dict transcript_id -> amino-acid sequence.
        excluded_names: Set of job names to skip (already submitted).
        allowed_transcripts: If provided, only these transcript IDs are
            considered (use to enforce a cutoff-pct filter).
    """
    excluded_names = excluded_names or set()
    jobs = []
    for gene_id in gene_order:
        if gene_id not in isoforms_by_gene:
            continue

        seen_sequences = set()
        for transcript_id, exons in isoforms_by_gene[gene_id].items():
            if allowed_transcripts is not None and transcript_id not in allowed_transcripts:
                continue
            if not transcript_has_cds(exons):
                continue

            seq = protein_sequences.get(transcript_id)
            if seq is None:
                continue

            if seq in seen_sequences:
                continue
            seen_sequences.add(seq)

            job_name = make_job_name(transcript_id, gene_id)
            if job_name in excluded_names:
                continue

            jobs.append({
                "name": job_name,
                "modelSeeds": [],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": seq,
                            "count": 1,
                            "useStructureTemplate": True
                        }
                    }
                ],
                "dialect": "alphafoldserver",
                "version": 1
            })

    return jobs


# ---------------------------------------------------------------------------
# File writing
# ---------------------------------------------------------------------------

def write_job_files(jobs, output_dir, batch_size=30, start_index=1):
    os.makedirs(output_dir, exist_ok=True)
    num_files = (len(jobs) + batch_size - 1) // batch_size
    written = []
    for i in range(num_files):
        batch = jobs[i * batch_size: (i + 1) * batch_size]
        filename = os.path.join(output_dir, f"alphafold_jobs_{start_index + i}.json")
        with open(filename, "w") as f:
            json.dump(batch, f, indent=1)
        written.append(filename)
        print(f"  Wrote {filename}  ({len(batch)} jobs)")
    return written


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    # Default paths are relative to the workspace root (two levels up)
    root = os.path.join(os.path.dirname(__file__), "../..")
    p = argparse.ArgumentParser(description="Create AlphaFold Server job JSON files")
    p.add_argument(
        "--input",
        default=os.path.join(root, "data/neuro_project/output/isoform_distributions/tables/distributions_condition_sum.tsv"),
        help="Isoform distributions TSV (gene_id, transcript_id, global_col, sample cols…)"
    )
    p.add_argument(
        "--gtf",
        default=os.path.join(root, "data/neuro_project/expressed_isoforms.gtf"),
        help="GTF file with exon/CDS annotations"
    )
    p.add_argument(
        "--fasta",
        default=os.path.join(root, "data/neuro_project/proteins.fasta"),
        help="Protein FASTA file"
    )
    p.add_argument(
        "--output",
        default=os.path.join(os.path.dirname(__file__), "../data/jobs"),
        help="Output directory for job JSON files"
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=30,
        help="Number of jobs per JSON file (default: 30)"
    )
    p.add_argument(
        "--cutoff-pct",
        type=float,
        default=0.0,
        help="Only include transcripts contributing >= this %% of their gene's total "
             "expression (default: 0 = no filter). Use 1.5 to match the dashboard default."
    )
    p.add_argument(
        "--exclude-submitted",
        type=int,
        default=None,
        metavar="MAX_INDEX",
        help="Skip any job whose name already appears in alphafold_jobs_1.json … "
             "alphafold_jobs_MAX_INDEX.json (including -doneMMDD variants). "
             "New files are numbered starting at MAX_INDEX + 1."
    )
    return p.parse_args()


def main():
    args = parse_args()

    # --- Validate inputs ---
    for path, label in [(args.input, "input TSV"), (args.gtf, "GTF"), (args.fasta, "FASTA")]:
        if not os.path.exists(path):
            print(f"Error: {label} file not found: {path}", file=sys.stderr)
            sys.exit(1)

    # --- Load distributions table ---
    print(f"Loading input TSV: {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    if "gene_id" not in df.columns or "transcript_id" not in df.columns:
        print("Error: input must have 'gene_id' and 'transcript_id' columns", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(df)} rows, {df['gene_id'].nunique()} genes")

    # --- Compute gene ranking ---
    print("Computing gene ranking...")
    gene_order = ranked_gene_order(df)
    print(f"  Ranked {len(gene_order)} genes")

    # --- Parse GTF ---
    print(f"Parsing GTF: {args.gtf}")
    isoforms_by_gene = parse_isoform_file(args.gtf)
    print(f"  Found {len(isoforms_by_gene)} genes with exon data")

    # --- Load protein sequences ---
    print(f"Loading protein sequences: {args.fasta}")
    protein_sequences = load_protein_sequences(args.fasta)
    print(f"  Loaded {len(protein_sequences)} sequences")

    # --- Apply cutoff filter ---
    allowed_transcripts = None
    if args.cutoff_pct > 0:
        from isoform_distribution.utils import get_filtered_isoforms as _get_filtered_isoforms
        print(f"Applying cutoff: keeping transcripts with >= {args.cutoff_pct}% of gene expression...")
        transcript_to_gene = load_gtf_mapping(args.gtf)
        df_raw = pd.read_csv(
            os.path.join(os.path.dirname(__file__), "../../",
                         "data/neuro_project/expressed_isoforms_matrix.txt"),
            sep="\t", index_col=0
        )
        df_raw = map_transcripts_to_genes(df_raw, transcript_to_gene)
        sample_cols = [c for c in df_raw.columns if c != 'gene_id' and 'fetal' not in c.lower()]
        allowed_transcripts = set()
        for gene_id, gene_block in df_raw.groupby('gene_id'):
            gene_t = gene_block.drop(columns=['gene_id'])
            kept = _get_filtered_isoforms(gene_t, sample_cols, args.cutoff_pct, stat='sum')
            if kept:
                allowed_transcripts.update(kept)
        print(f"  {len(allowed_transcripts)} transcripts pass the cutoff")

    # --- Collect already-submitted job names ---
    excluded_names: set = set()
    start_index = 1
    if args.exclude_submitted is not None:
        print(f"Collecting submitted job names from files 1–{args.exclude_submitted}...")
        excluded_names = collect_submitted_names(args.output, args.exclude_submitted)
        start_index = args.exclude_submitted + 1
        print(f"  {len(excluded_names)} job names will be skipped")
        print(f"  New files will be numbered starting at {start_index}")

    # --- Build jobs ---
    print("Building job entries...")
    jobs = build_jobs(gene_order, isoforms_by_gene, protein_sequences,
                      excluded_names=excluded_names,
                      allowed_transcripts=allowed_transcripts)
    print(f"  Total new jobs: {len(jobs)}")

    if not jobs:
        print("No jobs to write. All qualifying transcripts may already have been submitted.")
        sys.exit(0)

    # --- Write output ---
    print(f"Writing to: {args.output}  (batch size = {args.batch_size}, starting at index {start_index})")
    written = write_job_files(jobs, args.output, batch_size=args.batch_size, start_index=start_index)
    print(f"Done. {len(written)} file(s), {len(jobs)} job(s) total.")


if __name__ == "__main__":
    main()
