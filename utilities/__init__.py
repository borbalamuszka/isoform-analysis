"""Utility modules for gene expression analysis.

Modules:
- gtf_utils: GTF file parsing utilities
- compare_fasta_files: FASTA file comparison
- compare_confidence_tables: Confidence interval comparison
- convert_to_tsv: Format conversion utilities
- exon_visualisation: Exon visualization functions
- map_transcriptId_geneId: Transcript/gene ID mapping
"""

from .gtf_utils import (
    parse_gtf_attributes,
    UNKNOWN_GENE,
)

__all__ = [
    "parse_gtf_attributes",
    "UNKNOWN_GENE",
]
