"""GTF file parsing and exon structure handling.

This module handles:
- GTF attribute parsing
- Exon and CDS extraction
- Compressed coordinate mapping for visualization
"""
import re
import numpy as np
from collections import defaultdict


def parse_gtf_attributes(attributes):
    """Parse GTF attributes field to extract transcript_id, gene_id, and gene_name.
    
    Args:
        attributes: GTF attributes string
        
    Returns:
        Tuple of (transcript_id, gene_id, gene_name)
    """
    transcript_match = re.search(r'transcript_id "([^"]+)"', attributes)
    gene_match = re.search(r'gene_id "([^"]+)"', attributes)
    gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
    
    transcript_id = transcript_match.group(1) if transcript_match else None
    gene_id = gene_match.group(1) if gene_match else None
    gene_name = gene_name_match.group(1) if gene_name_match else None
    
    return transcript_id, gene_id, gene_name


def parse_gene_names(filename):
    """Parse GTF file to extract gene names.
    
    Args:
        filename: Path to GTF file
        
    Returns:
        Dictionary mapping gene_id -> gene_name
    """
    gene_names = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            
            attributes = parts[8]
            _, gene_id, gene_name = parse_gtf_attributes(attributes)
            
            if gene_id and gene_name:
                gene_names[gene_id] = gene_name
    
    return gene_names


def parse_isoform_file(filename):
    """Parse GTF file to extract exon structures with CDS information.
    
    GTF files have separate rows for 'exon' and 'CDS' features. We need to:
    1. Collect all exon entries for each transcript
    2. Collect all CDS entries for each transcript
    3. Merge CDS info into corresponding exons
    
    Args:
        filename: Path to GTF file
        
    Returns:
        Dictionary mapping gene_id -> transcript_id -> list of exon dicts
    """
    # First pass: collect exons and CDS separately
    exons_by_transcript = defaultdict(list)
    cds_by_transcript = defaultdict(list)
    gene_for_transcript = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            
            chrom = parts[0]
            feature_type = parts[2]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                raise ValueError(
                    f"Selected exons file does not appear to be a valid GTF format. "
                    f"Invalid integer coordinate values ('{parts[3]}', '{parts[4]}') at line: '{line[:100]}...'"
                )
            strand = parts[6] if len(parts) > 6 else '+'
            attributes = parts[8]
            
            transcript_id, gene_id, _ = parse_gtf_attributes(attributes)
            
            if not transcript_id or not gene_id:
                continue
            
            gene_for_transcript[transcript_id] = gene_id
            
            if feature_type == 'exon':
                exons_by_transcript[transcript_id].append({
                    'start': start,
                    'end': end,
                    'strand': strand,
                })
            elif feature_type == 'CDS':
                cds_by_transcript[transcript_id].append({
                    'start': start,
                    'end': end
                })
    
    # Second pass: merge CDS into exons
    isoforms_by_gene = defaultdict(lambda: defaultdict(list))
    
    for transcript_id, exons in exons_by_transcript.items():
        gene_id = gene_for_transcript.get(transcript_id)
        if not gene_id:
            continue
        
        cds_regions = cds_by_transcript.get(transcript_id, [])
        
        for exon in exons:
            exon_start = exon['start']
            exon_end = exon['end']
            
            # Find overlapping CDS regions for this exon
            cds_start = None
            cds_end = None
            
            for cds in cds_regions:
                # Check if CDS overlaps with this exon
                if cds['start'] <= exon_end and cds['end'] >= exon_start:
                    # Find the intersection
                    overlap_start = max(cds['start'], exon_start)
                    overlap_end = min(cds['end'], exon_end)
                    
                    if cds_start is None:
                        cds_start = overlap_start
                        cds_end = overlap_end
                    else:
                        # Extend CDS range if multiple CDS regions overlap
                        cds_start = min(cds_start, overlap_start)
                        cds_end = max(cds_end, overlap_end)
            
            isoforms_by_gene[gene_id][transcript_id].append({
                "exon_start": exon_start,
                "exon_end": exon_end,
                "cds_start": cds_start,
                "cds_end": cds_end,
                "strand": exon.get("strand", "+"),
            })
    
    return isoforms_by_gene


def build_compressed_mapping(isoforms_dict, exon_scale_bp_per_unit=50.0, intron_width_units=5.0):
    """Build a non-linear mapping from genomic to compressed coordinates for a gene's isoforms.
    
    This creates a compressed visualization where exons are drawn to scale and introns
    are compressed to a fixed width.
    
    Args:
        isoforms_dict: Dict mapping transcript_id -> list of exon dicts
        exon_scale_bp_per_unit: How many bp per visual unit for exons
        intron_width_units: Fixed width for intron segments
    
    Returns:
        Dict with mapping function and metadata
    """
    # 1) Gather all exon boundaries
    all_boundaries = set()
    for transcript_exons in isoforms_dict.values():
        for exon in transcript_exons:
            all_boundaries.add(exon["exon_start"])
            all_boundaries.add(exon["exon_end"])
    
    breakpoints = np.array(sorted(all_boundaries))
    
    # 2) Define segments between breakpoints
    seg_starts = breakpoints[:-1]
    seg_ends = breakpoints[1:]
    seg_len_bp = seg_ends - seg_starts
    n_seg = len(seg_starts)
    
    # 3) Determine which segments are exonic (covered by at least one isoform)
    is_exon_seg = np.zeros(n_seg, dtype=bool)
    for transcript_exons in isoforms_dict.values():
        for exon in transcript_exons:
            s, e = exon["exon_start"], exon["exon_end"]
            # Find segment indices overlapped by this exon
            i0 = np.searchsorted(breakpoints, s, side='left')
            i1 = np.searchsorted(breakpoints, e, side='right') - 1
            i0 = max(min(i0, n_seg - 1), 0)
            i1 = max(min(i1, n_seg), 0)
            is_exon_seg[i0:i1] = True
    
    # 4) Assign visual widths
    widths = np.empty(n_seg, dtype=float)
    widths[is_exon_seg] = seg_len_bp[is_exon_seg] / exon_scale_bp_per_unit
    widths[~is_exon_seg] = intron_width_units
    
    # 5) Cumulative compressed coordinates
    cum = np.concatenate([[0.0], np.cumsum(widths)])
    
    def map_coord(x):
        """Map genomic coordinate(s) to compressed coordinates."""
        x = np.asarray(x)
        scalar_input = x.ndim == 0
        x = np.atleast_1d(x)
        
        idx = np.searchsorted(breakpoints, x, side='right') - 1
        idx = np.clip(idx, 0, n_seg - 1)
        
        # Fraction within segment
        with np.errstate(divide='ignore', invalid='ignore'):
            frac = (x - seg_starts[idx]) / seg_len_bp[idx]
            frac = np.nan_to_num(frac, nan=0.0)
        
        result = cum[idx] + frac * widths[idx]
        return result[0] if scalar_input else result
    
    return {
        "breakpoints": breakpoints,
        "seg_starts": seg_starts,
        "seg_ends": seg_ends,
        "seg_widths": widths,
        "cum": cum,
        "is_exon_seg": is_exon_seg,
        "map_coord": map_coord,
        "total_width": cum[-1]
    }


def parse_parentless_transcripts(filename):
    """Parse GTF file to find transcript IDs that have no gene_id."""
    import re
    seen_transcripts = set()
    mapped_transcripts = set()
    with open(filename, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 9:
                    attributes = parts[8]
                    transcript_match = re.search(r'transcript_id "([^"]+)"', attributes)
                    gene_match = re.search(r'gene_id "([^"]+)"', attributes)
                    if transcript_match:
                        tx_id = transcript_match.group(1)
                        seen_transcripts.add(tx_id)
                        if gene_match and gene_match.group(1):
                            mapped_transcripts.add(tx_id)
    return seen_transcripts - mapped_transcripts
