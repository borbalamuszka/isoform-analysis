"""Isoform dashboard package.

A Dash web application for exploring isoform entropy, correlations, and gene structure.

Modules:
- dashboard_app: Main entry point
- app_layout: Dashboard layout and callbacks
- data_processing: Data analysis functions
- gtf_parser: GTF file parsing
- visualizations: Plotting functions
- alphafold_geometry: structure visualization
- config: Configuration and styling
"""

from .data_processing import (
    calculate_entropy_and_correlation,
    compute_gene_ranking,
    compute_min_spearman_per_gene,
    prepare_table_data,
)
from .gtf_parser import (
    parse_isoform_file,
    parse_gene_names,
    parse_parentless_transcripts,
)
from .app_layout import create_app
from .config import Colors, Dimensions, Styles

__all__ = [
    "calculate_entropy_and_correlation",
    "compute_gene_ranking",
    "compute_min_spearman_per_gene",
    "prepare_table_data",
    "parse_isoform_file",
    "parse_gene_names",
    "parse_parentless_transcripts",
    "create_app",
    "Colors",
    "Dimensions",
    "Styles",
]
