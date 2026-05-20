#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd

def read_ci_table(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    # Try with first column as index (expected "isoform")
    df = pd.read_csv(path, sep="\t", index_col=0)
    # Require columns
    for col in ("ci_lower", "ci_upper"):
        if col not in df.columns:
            raise ValueError(f"Missing required column '{col}' in {path}")
    # Compute width
    df = df.assign(width=df["ci_upper"] - df["ci_lower"])
    return df[["width"]]

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Compare CI widths between two result tables")
    p.add_argument("--all-ci", required=True, help="Path to CI table computed using all samples (n=45)")
    p.add_argument("--control-ci", required=True, help="Path to CI table computed using control-only samples (n=11)")
    p.add_argument("--out", default="ci_width_diff.tsv", help="Output TSV with width differences")
    args = p.parse_args(argv)

    all_df = read_ci_table(args.all_ci).rename(columns={"width": "width_all"})
    ctrl_df = read_ci_table(args.control_ci).rename(columns={"width": "width_control"})

    # Inner join on isoform index
    merged = all_df.join(ctrl_df, how="inner")
    if merged.empty:
        print("No overlapping isoforms between the two tables.", file=sys.stderr)
        return 2

    merged["diff"] = merged["width_control"] - merged["width_all"]          # positive => control CI wider
    merged["abs_diff"] = merged["diff"].abs()
    merged["pct_change"] = merged["diff"] / merged["width_all"]             # relative to all-sample width

    # Sort by largest widening under control-only
    merged = merged.sort_values("diff", ascending=False)

    merged.to_csv(args.out, sep="\t", index=True, header=True)
    print(f"Wrote width comparison for {len(merged)} isoforms to: {args.out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())