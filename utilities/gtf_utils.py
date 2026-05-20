"""
Utility functions for parsing and processing GTF (Gene Transfer Format) files.
"""

import re

UNKNOWN_GENE = 'UNKNOWN_GENE'


def parse_gtf_attributes(attributes):
    """
    Parse GTF attributes field to extract transcript_id and gene_id.
    
    Args:
        attributes: GTF attributes string (9th field in GTF format)
        
    Returns:
        Tuple of (transcript_id, gene_id)
    """
    transcript_match = re.search(r'transcript_id "([^"]+)"', attributes)
    gene_match = re.search(r'gene_id "([^"]+)"', attributes)
    
    transcript_id = transcript_match.group(1) if transcript_match else None
    gene_id = gene_match.group(1) if gene_match else UNKNOWN_GENE
    
    return transcript_id, gene_id
