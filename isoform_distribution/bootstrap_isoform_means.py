#!/usr/bin/env python3
"""
Bootstrap mean expression per isoform and 95% confidence intervals.

Steps:
1) Load isoform expression matrix (TSV). Rows = isoforms/transcripts, columns = samples.
2) Keep only samples whose column name contains the provided include key (--include-key, default "adult", case-insensitive).
3) For global bootstrapping:
   - For each iteration (N): sample S columns with replacement, where S is the number of columns matching --include-key.
   - Compute the mean expression per isoform across the resampled columns.
   - After N iterations, sort the bootstrap means for each isoform, drop the lowest/highest
     tail_count values (default for N=1000: 25 at each end) to form a 95% CI.
4) For grouped bootstrapping:
   - Group samples by brain region (e.g., Caudate, DLPFC, Hippocampus)
   - Group samples by condition (e.g., Control, MDD, Bipolar, Schizo)
   - For each group, perform bootstrapping using only samples within that group
   - Resample with replacement using the same number of samples as in the original group
   - Compute 95% CIs for each group separately
5) Save results: global CIs plus region-specific and condition-specific CIs to a single TSV file.

Output columns:
  - isoform: Transcript/isoform ID
  - ci_lower, ci_upper: Global 95% CI bounds across all samples
  - ci_lower_{region}, ci_upper_{region}: 95% CI bounds for each brain region
  - ci_lower_{condition}, ci_upper_{condition}: 95% CI bounds for each condition

Usage:
  python bootstrap_isoform_means.py \
    --input data/neuro_project/expressed_isoforms_matrix.tsv \
    --output-dir data/neuro_project/output/bootstrap/ \
    --iterations 1000 \
    --seed 42 \
    --include-key adult

Notes:
- Input is expected to have isoform IDs in the first column or as the dataframe index.
- If the first column is unnamed or named like an ID column, it will be set as index.
- Sample names should follow the format: Region_Condition_Adult_ID (e.g., DLPFC_Control_Adult_R12345)
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Tuple, Dict
from collections import defaultdict
import numpy as np
import pandas as pd
from .utils import parse_sample_name

# Define base directory relative to this script
BASE_DIR = Path(__file__).parent.parent  # source/


def read_matrix(path: str) -> pd.DataFrame:
    """Read expression matrix as DataFrame with isoforms as index.
    Input is TSV-only. If the first column looks like IDs, sets it as index.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    # TSV-only
    df = pd.read_csv(path, sep="\t")  # do NOT force dtype=float here

    # If first column is non-numeric or named like an id, treat it as index
    first_col = df.columns[0]
    if (df[first_col].dtype == object) or (first_col.lower() in {"id", "isoform", "transcript"}):
        df = df.set_index(first_col)

    # Ensure remaining columns are numeric
    for c in df.columns:
        if not np.issubdtype(df[c].dtype, np.number):
            # Try convert
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def bootstrap_means(
    df: pd.DataFrame,
    iterations: int = 1000,
    sample_size: int = 45,
    random_state: int | None = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute bootstrap means per isoform.
    Returns
    - sorted_means: numpy array (n_isoforms, iterations) sorted along axis=1
    - overall_mean: numpy array (n_isoforms,) mean across iterations
    """
    rng = np.random.default_rng(random_state)

    n_samples = df.shape[1]
    if n_samples == 0:
        raise ValueError("No columns remain after filtering.")

    # Pre-allocate results: rows = isoforms, cols = iterations
    means = np.empty((df.shape[0], iterations), dtype=float)

    # Convert to numpy for speed
    X = df.to_numpy()  # shape (n_isoforms, n_samples)

    for b in range(iterations):
        # Sample column indices with replacement
        idx = rng.integers(0, n_samples, size=sample_size)
        # Take subset and compute mean along columns
        subset = X[:, idx]
        means[:, b] = subset.mean(axis=1)

    # Also produce sorted means for CI computation
    sorted_means = np.sort(means, axis=1)
    # Overall mean across iterations
    overall_mean = means.mean(axis=1)
    return sorted_means, overall_mean


def compute_ci(sorted_means: np.ndarray, iterations: int, drop_each_side: int) -> Tuple[np.ndarray, np.ndarray]:
    """Compute CI bounds by dropping `drop_each_side` from each end on sorted bootstrap means."""
    if drop_each_side * 2 >= iterations:
        raise ValueError("drop_each_side is too large relative to iterations")
    lower_idx = drop_each_side
    upper_idx = iterations - drop_each_side - 1
    lower = sorted_means[:, lower_idx]
    upper = sorted_means[:, upper_idx]
    return lower, upper


def group_samples_by_attribute(df: pd.DataFrame, attribute: str = 'region') -> Dict[str, list]:
    """
    Group sample columns by brain region or condition.
    
    Args:
        df: DataFrame with sample columns
        attribute: 'region' or 'condition'
        
    Returns:
        Dictionary mapping group name to list of column names
    """
    groups = defaultdict(list)
    for col in df.columns:
        brain_region, condition = parse_sample_name(col)
        if brain_region and condition:
            group_key = brain_region if attribute == 'region' else condition
            groups[group_key].append(col)
    return dict(groups)


def bootstrap_grouped(
    df: pd.DataFrame,
    groups: Dict[str, list],
    iterations: int = 1000,
    random_state: int | None = 42,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Perform bootstrapping separately for each group.
    
    Args:
        df: DataFrame with isoform expression data
        groups: Dictionary mapping group name to list of column names
        iterations: Number of bootstrap iterations
        random_state: Random seed
        
    Returns:
        Dictionary mapping group name to (ci_lower, ci_upper) arrays
    """
    rng = np.random.default_rng(random_state)
    results = {}
    
    # Tail drop count for 95% CI
    tail_drop = int(np.ceil(iterations * 0.025))
    
    for group_name, group_cols in groups.items():
        if len(group_cols) < 1:
            continue
            
        # Get data for this group
        group_df = df[group_cols]
        n_samples = len(group_cols)
        
        # Pre-allocate results: rows = isoforms, cols = iterations
        means = np.empty((group_df.shape[0], iterations), dtype=float)
        
        # Convert to numpy for speed
        X = group_df.to_numpy()  # shape (n_isoforms, n_samples)
        
        for b in range(iterations):
            # Sample column indices with replacement
            idx = rng.integers(0, n_samples, size=n_samples)
            # Take subset and compute mean along columns
            subset = X[:, idx]
            means[:, b] = subset.mean(axis=1)
        
        # Sort means for CI computation
        sorted_means = np.sort(means, axis=1)
        
        # Compute CI
        lower, upper = compute_ci(sorted_means, iterations=iterations, drop_each_side=tail_drop)
        results[group_name] = (lower, upper)
    
    return results


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Bootstrap isoform mean expression and confidence intervals")
    p.add_argument("--input", required=False, default=str(BASE_DIR / "data/neuro_project/expressed_isoforms_matrix.tsv"), help="Path to input expression matrix (TSV)")
    p.add_argument("--output-dir", required=False, default=str(BASE_DIR / "data/neuro_project/output/isoform_distributions"), help="Directory to write outputs")
    p.add_argument("--iterations", type=int, default=1000, help="Number of bootstrap iterations N")
    p.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    p.add_argument("--include-key", default="adult", help="Keyword to include columns (case-insensitive contains)")
    args = p.parse_args(argv)

    df = read_matrix(args.input)
    # Filter to only columns containing the include key
    include_key = args.include_key
    mask = df.columns.str.contains(include_key, case=False, regex=False)
    filtered_df = df.loc[:, mask]

    if filtered_df.shape[1] < 1:
        print(f"No columns available after filtering for '{include_key}'. Check input file.", file=sys.stderr)
        return 2

    # Adjust sample size to the number of matched columns
    sample_size = filtered_df.shape[1]

    # Bootstrap means
    sorted_means, overall_mean = bootstrap_means(filtered_df, iterations=args.iterations, sample_size=sample_size, random_state=args.seed)

    # Tail drop count for 95% CI: 2.5% each side => 1000*0.025 = 25
    tail_drop = int(np.ceil(args.iterations * 0.025))

    lower, upper = compute_ci(sorted_means, iterations=args.iterations, drop_each_side=tail_drop)

    # Grouped bootstrapping by region
    print("Computing bootstrapped CIs by brain region...")
    region_groups = group_samples_by_attribute(filtered_df, attribute='region')
    region_cis = bootstrap_grouped(filtered_df, region_groups, iterations=args.iterations, random_state=args.seed)
    
    # Grouped bootstrapping by condition
    print("Computing bootstrapped CIs by condition...")
    condition_groups = group_samples_by_attribute(filtered_df, attribute='condition')
    condition_cis = bootstrap_grouped(filtered_df, condition_groups, iterations=args.iterations, random_state=args.seed)

    os.makedirs(args.output_dir, exist_ok=True)

    # Write confidence intervals with grouped CIs
    ci_out = os.path.join(args.output_dir, "confidence_intervals.tsv")

    ci_data = {
        "isoform": filtered_df.index,
        "ci_lower": lower,
        "ci_upper": upper,
    }
    
    # Add region-specific CI columns
    for region_name in sorted(region_cis.keys()):
        ci_lower, ci_upper = region_cis[region_name]
        ci_data[f"ci_lower_{region_name}"] = ci_lower
        ci_data[f"ci_upper_{region_name}"] = ci_upper
    
    # Add condition-specific CI columns
    for condition_name in sorted(condition_cis.keys()):
        ci_lower, ci_upper = condition_cis[condition_name]
        ci_data[f"ci_lower_{condition_name}"] = ci_lower
        ci_data[f"ci_upper_{condition_name}"] = ci_upper
    
    ci_df = pd.DataFrame(ci_data).set_index("isoform")
    ci_df.to_csv(ci_out, sep="\t")

    print(f"Wrote confidence intervals to: {ci_out}")
    print(f"Included CIs for {len(region_cis)} regions: {sorted(region_cis.keys())}")
    print(f"Included CIs for {len(condition_cis)} conditions: {sorted(condition_cis.keys())}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
