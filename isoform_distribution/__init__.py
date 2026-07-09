"""Isoform distribution analysis package.

Modules:
- distributions: Main script for generating distribution tables
- plots: Visualization of distributions
- correlation: Correlation and entropy analysis
- utils: Common utility functions
- bootstrap_isoform_means: Bootstrap confidence intervals
"""

from .utils import (
    parse_sample_name,
    aggregate_samples_by_group,
    get_filtered_isoforms,
    prepare_gene_list_and_paths,
)

__all__ = [
    "parse_sample_name",
    "aggregate_samples_by_group",
    "get_filtered_isoforms",
    "prepare_gene_list_and_paths",
]
