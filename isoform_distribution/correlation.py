#!/usr/bin/env python3
"""
Isoform entropy, TVD, and between-sample correlation per gene.

Summary:
- For each gene (across its isoforms), computes:
  - top_isoform_entropy: Shannon entropy across samples of the most abundant isoform (by global_sum).
  - summed_isoform_entropy: Shannon entropy across samples of the gene-level sum across isoforms.
  - Pairwise Spearman correlations over isoform-level abundances (all sample pairs).
  - Pairwise Total Variation Distance (TVD) over per-gene, per-sample normalized isoform distributions (normalization done internally).
  - Aggregates per gene: min_spearman (minimum across pairs) and max_tvd (maximum across pairs).
- Writes a TSV containing per-gene metrics (pairwise + min/max).
- Writes interactive HTML scatter plots (axes use min/max where applicable):
  - interactive_entropy_vs_min_spearman.html            (Top isoform entropy vs Min Spearman)
  - interactive_summed_entropy_vs_min_spearman.html     (Summed isoform entropy vs Min Spearman)
  - interactive_entropy_vs_max_tvd.html                 (Top isoform entropy vs Max TVD)
  - interactive_summed_entropy_vs_max_tvd.html          (Summed isoform entropy vs Max TVD)
  - interactive_min_spearman_vs_max_tvd.html            (Min Spearman vs Max TVD)
  - interactive_summed_vs_top_entropy_min_spearman_colored.html  (Summed vs Top entropy; color=Min Spearman, red→grey→green)
  - interactive_summed_vs_top_entropy_max_tvd_colored.html       (Summed vs Top entropy; color=Max TVD, green→grey→red; high TVD shown in red)

Input:
- A TSV with columns:
  - gene_id, transcript_id, global_sum
  - >= 2 numeric sample columns (e.g., S1, S2, S3, ...)
- Input sample columns need not be normalized; TVD normalizes per gene per sample internally.

Notes:
- Requires: pandas, numpy, scipy, plotly.
- Needs at least two numeric sample columns.
- Correlations on constant arrays are marked "ConstantInput" and excluded from interactive plots.
"""

import sys
import os
import pandas as pd
import numpy as np
from itertools import combinations
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import plotly.express as px
from matplotlib.backends.backend_pdf import PdfPages

def calculate_entropy_and_correlation(df, sample_cols):
    pairs = list(combinations(sample_cols, 2))
    results = []

    for gene, sub in df.groupby("gene_id", sort=False):
        row = {"gene_id": gene, "n_isoforms": len(sub)}

        # Top isoform entropy (across samples)
        top_idx = sub["global_sum"].idxmax()
        top_isoform = sub.loc[top_idx, sample_cols].astype(float)
        # Ensure normalized across samples for entropy of the top isoform row
        top_sum = top_isoform.sum()
        if top_sum <= 0:
            probs = np.zeros_like(top_isoform, dtype=float)
        else:
            probs = top_isoform / top_sum
        entropy = -np.sum(probs * np.log(probs + 1e-12))  # avoid log(0)
        row["top_isoform_entropy"] = round(entropy, 2)

        # Summed isoform entropy (across samples, summed over isoforms)
        sample_sums = sub[sample_cols].astype(float).sum(axis=0)
        total_sum = sample_sums.sum()
        if total_sum <= 0:
            sample_probs = np.zeros_like(sample_sums, dtype=float)
        else:
            sample_probs = sample_sums / total_sum
        summed_entropy = -np.sum(sample_probs * np.log(sample_probs + 1e-12))
        row["summed_isoform_entropy"] = round(summed_entropy, 2)

        # Per-sample isoform distributions within this gene for TVD/correlation
        # Each column becomes a probability vector over isoforms of the gene.
        sub_probs = sub[sample_cols].astype(float).copy()
        col_sums = sub_probs.sum(axis=0)
        # Normalize each sample column within the gene; if zero, keep zeros
        safe_sums = col_sums.replace(0, np.nan)
        sub_probs = sub_probs.divide(safe_sums, axis=1).fillna(0.0)

        for a, b in pairs:
            p = sub_probs[a].to_numpy()
            q = sub_probs[b].to_numpy()

            # TVD: 0.5 * sum |p_i - q_i|
            tvd = 0.5 * float(np.sum(np.abs(p - q)))
            row[f"TVD {a} vs {b}"] = round(tvd, 4)

            # Spearman correlation on raw values (not probabilities) across isoforms
            x = sub[a].tolist()
            y = sub[b].tolist()
            if len(set(x)) == 1 or len(set(y)) == 1:
                row[f"{a} vs {b}"] = "ConstantInput"
            else:
                corr = spearmanr(x, y).correlation
                row[f"{a} vs {b}"] = round(corr, 2)
        results.append(row)
    return results

def _write_scatter(df_plot, x, y, title, out_dir, file_name, color=None, hover_data=None, color_continuous_scale=None, range_color=None):
    if df_plot.empty:
        return
    if hover_data is None:
        hover_data = ["Number of Isoforms"]

    # Ensure Gene ID is present
    if "Gene ID" not in df_plot.columns:
        if "gene_id" in df_plot.columns:
            df_plot["Gene ID"] = df_plot["gene_id"]
        else:
            df_plot["Gene ID"] = ""

    fig = px.scatter(
        df_plot,
        x=x,
        y=y,
        color=color,
        hover_name="Gene ID",
        size="Number of Isoforms",
        hover_data=hover_data,
        title=title,
        color_continuous_scale=color_continuous_scale,
        range_color=range_color,
    )
    
    # Add Gene ID as a text attribute that JavaScript can access
    # Use mode='markers' to hide text from display but keep it accessible
    fig.update_traces(text=df_plot["Gene ID"].values, mode='markers')
    
    # Make hover readable
    fig.update_traces(hoverlabel=dict(bgcolor="#ADD8E6", font_color="black"))

    # Read external JS post_script from assets folder
    script_dir = os.path.dirname(os.path.abspath(__file__))
    js_path = os.path.join(script_dir, "assets", "plotly_click_copy_gene_id.js")
    try:
        with open(js_path, "r", encoding="utf-8") as f:
            post_script = f.read()
    except Exception:
        post_script = None

    os.makedirs(out_dir, exist_ok=True)
    fig.write_html(
        os.path.join(out_dir, file_name),
        include_plotlyjs="cdn",
        full_html=True,
        post_script=post_script,
    )

def _summed_vs_top_entropy_colored(results, out_dir, color_field, title, outfile):
    if not results:
        return

    rows = _collect_min_spearman_max_tvd(results)
    plot_rows = []
    for r in rows:
        plot_rows.append({
            "Top Isoform Entropy": r.get("top_isoform_entropy"),
            "Summed Isoform Entropy": r.get("summed_isoform_entropy"),
            "Min Spearman": r.get("min_spearman"),
            "Max TVD": r.get("max_tvd"),
            "Gene ID": r["gene_id"],
            "Number of Isoforms": r["n_isoforms"],
        })

    df_plot = pd.DataFrame(plot_rows)
    if df_plot.empty:
        return

    # Colorscale
    if color_field == "Min Spearman":
        # Spearman: red (low) -> grey (mid) -> green (high)
        colorscale = [[0.0, "red"], [0.5, "grey"], [1.0, "green"]]
        range_color = (-1.0, 1.0)
    else:
        # TVD reversed: green (low) -> grey (mid) -> red (high)
        colorscale = [[0.0, "green"], [0.5, "grey"], [1.0, "red"]]
        range_color = (0.0, 1.0)

    _write_scatter(
        df_plot,
        x="Summed Isoform Entropy",
        y="Top Isoform Entropy",
        title=title,
        out_dir=out_dir,
        file_name=outfile,
        color=color_field,
        hover_data=["Number of Isoforms", "Min Spearman", "Max TVD"],
        color_continuous_scale=colorscale,
        range_color=range_color,
    )

def interactive_summed_vs_top_entropy_colored_by_min_spearman(results, out_dir):
    _summed_vs_top_entropy_colored(
        results,
        out_dir,
        color_field="Min Spearman",
        title="Summed vs Top Isoform Entropy (colored by min Spearman per gene)",
        outfile="interactive_summed_vs_top_entropy_min_spearman_colored.html",
    )

def interactive_summed_vs_top_entropy_colored_by_max_tvd(results, out_dir):
    _summed_vs_top_entropy_colored(
        results,
        out_dir,
        color_field="Max TVD",
        title="Summed vs Top Isoform Entropy (colored by max TVD per gene)",
        outfile="interactive_summed_vs_top_entropy_max_tvd_colored.html",
    )

def _collect_min_spearman_max_tvd(results):
    """
    For each gene result dict, compute and return:
    - min_spearman: minimum Spearman across all pairs (float or None)
    - max_tvd: maximum TVD across all pairs (float or None)
    Returns a list of dicts with keys:
      gene_id, n_isoforms, top_isoform_entropy, summed_isoform_entropy,
      min_spearman, max_tvd
    """
    if not results:
        return []

    # Identify column sets once
    spearman_cols = [col for col in results[0] if " vs " in col and not col.startswith("TVD ")]
    tvd_cols = [col for col in results[0] if col.startswith("TVD ") and " vs " in col]

    out = []
    for row in results:
        s_vals = [row.get(c) for c in spearman_cols if isinstance(row.get(c), (int, float))]
        t_vals = [row.get(c) for c in tvd_cols if isinstance(row.get(c), (int, float))]
        out.append({
            "gene_id": row["gene_id"],
            "n_isoforms": row["n_isoforms"],
            "top_isoform_entropy": row.get("top_isoform_entropy"),
            "summed_isoform_entropy": row.get("summed_isoform_entropy"),
            "min_spearman": float(min(s_vals)) if s_vals else None,
            "max_tvd": float(max(t_vals)) if t_vals else None,
        })
    return out

def _entropy_vs_metric(results, out_dir, entropy_col, metric_col, x_label, title, filename):
    """
    Generic entropy vs metric scatter using min Spearman or max TVD.
    metric_col must be one of: 'min_spearman', 'max_tvd'
    """
    rows = _collect_min_spearman_max_tvd(results)
    plot_rows = []
    for r in rows:
        entropy_val = r.get(entropy_col)
        metric_val = r.get(metric_col)
        if isinstance(entropy_val, (int, float)) and isinstance(metric_val, (int, float)):
            plot_rows.append({
                x_label: metric_val,
                "Entropy": entropy_val,
                "Gene ID": r["gene_id"],
                "Number of Isoforms": r["n_isoforms"],
            })
    df_plot = pd.DataFrame(plot_rows)
    _write_scatter(
        df_plot,
        x=x_label,
        y="Entropy",
        title=title,
        out_dir=out_dir,
        file_name=filename,
        color=None,
    )

def interactive_entropy_vs_correlation(results, out_dir):
    _entropy_vs_metric(
        results,
        out_dir,
        entropy_col="top_isoform_entropy",
        metric_col="min_spearman",
        x_label="Min Spearman",
        title="Top Isoform Entropy vs. Min Spearman",
        filename="interactive_entropy_vs_min_spearman.html",
    )

def interactive_summed_entropy_vs_correlation(results, out_dir):
    _entropy_vs_metric(
        results,
        out_dir,
        entropy_col="summed_isoform_entropy",
        metric_col="min_spearman",
        x_label="Min Spearman",
        title="Summed Isoform Entropy vs. Min Spearman",
        filename="interactive_summed_entropy_vs_min_spearman.html",
    )

def interactive_entropy_vs_tvd(results, out_dir):
    _entropy_vs_metric(
        results,
        out_dir,
        entropy_col="top_isoform_entropy",
        metric_col="max_tvd",
        x_label="Max TVD",
        title="Top Isoform Entropy vs. Max TVD",
        filename="interactive_entropy_vs_max_tvd.html",
    )

def interactive_summed_entropy_vs_tvd(results, out_dir):
    _entropy_vs_metric(
        results,
        out_dir,
        entropy_col="summed_isoform_entropy",
        metric_col="max_tvd",
        x_label="Max TVD",
        title="Summed Isoform Entropy vs. Max TVD",
        filename="interactive_summed_entropy_vs_max_tvd.html",
    )

def interactive_min_spearman_vs_max_tvd(results, out_dir):
    if not results:
        return

    # Identify column sets
    spearman_cols = [col for col in results[0] if " vs " in col and not col.startswith("TVD ")]
    tvd_cols = [col for col in results[0] if col.startswith("TVD ") and " vs " in col]
    if not spearman_cols or not tvd_cols:
        return

    rows = []
    for row in results:
        gene_id = row["gene_id"]
        n_isoforms = row["n_isoforms"]

        # Min Spearman across pairs (ignore non-floats)
        spearman_vals = [row[c] for c in spearman_cols if isinstance(row.get(c), (int, float))]
        # Max TVD across pairs (ignore non-floats)
        tvd_vals = [row[c] for c in tvd_cols if isinstance(row.get(c), (int, float))]

        if not spearman_vals or not tvd_vals:
            continue

        rows.append({
            "Min Spearman": float(min(spearman_vals)),
            "Max TVD": float(max(tvd_vals)),
            "Gene ID": gene_id,
            "Number of Isoforms": n_isoforms,
        })

    df_plot = pd.DataFrame(rows)
    _write_scatter(
        df_plot,
        x="Min Spearman",
        y="Max TVD",
        title="Minimum Spearman vs. Maximum TVD per Gene",
        out_dir=out_dir,
        file_name="interactive_min_spearman_vs_max_tvd.html",
        color=None,
    )

def compute_gene_ranking(results):
    """
    Rank genes based on:
    1. Negative min Spearman genes come before positive ones
    2. Within each group, rank by minimum of (top_isoform_entropy, summed_isoform_entropy)
       Lower min entropy = higher rank (lower rank number)
    
    Returns a dict mapping gene_id to rank number (1-indexed)
    """
    if not results:
        return {}
    
    # Identify Spearman correlation columns
    spearman_cols = [col for col in results[0] if " vs " in col and not col.startswith("TVD ")]
    
    # Build a list of dicts with necessary info
    gene_data = []
    for row in results:
        # Compute min Spearman
        s_vals = [row.get(c) for c in spearman_cols if isinstance(row.get(c), (int, float))]
        min_spearman = float(min(s_vals)) if s_vals else None
        
        # Compute min entropy
        top_entropy = row.get("top_isoform_entropy")
        summed_entropy = row.get("summed_isoform_entropy")
        min_entropy = None
        if isinstance(top_entropy, (int, float)) and isinstance(summed_entropy, (int, float)):
            min_entropy = min(top_entropy, summed_entropy)
        elif isinstance(top_entropy, (int, float)):
            min_entropy = top_entropy
        elif isinstance(summed_entropy, (int, float)):
            min_entropy = summed_entropy
        
        gene_data.append({
            "gene_id": row["gene_id"],
            "min_spearman": min_spearman,
            "min_entropy": min_entropy,
            "is_negative_spearman": min_spearman < 0 if min_spearman is not None else False
        })
    
    # Sort: negative Spearman first (True before False), then by min_entropy ascending
    # Handle None values in sorting
    def sort_key(item):
        neg_spearman = item["is_negative_spearman"]
        min_ent = item["min_entropy"]
        # Put items with None min_entropy at the end
        return (not neg_spearman, min_ent if min_ent is not None else float('inf'))
    
    gene_data_sorted = sorted(gene_data, key=sort_key)
    
    # Assign ranks (1-indexed)
    rank_map = {item["gene_id"]: rank + 1 for rank, item in enumerate(gene_data_sorted)}
    
    return rank_map


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {os.path.basename(sys.argv[0])} <input_table.tsv>", file=sys.stderr)
        sys.exit(1)
    in_path = sys.argv[1]
    try:
        df = pd.read_csv(in_path, sep="\t")
    except Exception as e:
        print(f"Failed to read input: {e}", file=sys.stderr)
        sys.exit(1)

    if "gene_id" not in df.columns:
        print("Expected column 'gene_id' not found.", file=sys.stderr)
        sys.exit(1)

    meta = {"gene_id", "transcript_id", "global_sum"}
    sample_cols = [c for c in df.columns
                   if c not in meta and pd.api.types.is_numeric_dtype(df[c])]
    if len(sample_cols) < 2:
        print("Need at least two sample columns for correlations.", file=sys.stderr)
        sys.exit(1)

    results = calculate_entropy_and_correlation(df, sample_cols)

    # Compute gene ranking
    rank_map = compute_gene_ranking(results)
    
    # Add rank to each result row
    for row in results:
        row["rank"] = rank_map.get(row["gene_id"], None)

    # Filter out rows where any correlation column is not a float
    corr_cols = [col for col in results[0] if " vs " in col]
    filtered_results = [
        row for row in results
        if all(isinstance(row[col], float) for col in corr_cols)
    ]

    out_path = os.path.splitext(in_path)[0] + "_gene_correlations.tsv"
    pd.DataFrame(results).to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {out_path}")

    # Use the same directory as the input table for all outputs
    input_dir = os.path.dirname(os.path.abspath(in_path))
    out_dir = input_dir

    # interactive_entropy_vs_correlation(filtered_results, out_dir)
    # interactive_summed_entropy_vs_correlation(filtered_results, out_dir)
    # interactive_entropy_vs_tvd(filtered_results, out_dir)
    # interactive_summed_entropy_vs_tvd(filtered_results, out_dir)
    # interactive_min_spearman_vs_max_tvd(filtered_results, out_dir)

    # Color-mapped plots
    interactive_summed_vs_top_entropy_colored_by_min_spearman(filtered_results, out_dir)
    interactive_summed_vs_top_entropy_colored_by_max_tvd(filtered_results, out_dir)

if __name__ == "__main__":
    main()