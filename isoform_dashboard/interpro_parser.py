"""InterPro domain parsing and visualization utilities.

This module handles:
- Loading InterPro JSON results
- Converting protein coordinates to genomic coordinates
- Extracting domain information for specific transcripts
- Preparing domain data for visualization overlay on exons
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

log = logging.getLogger(__name__)


def _find_interpro_result(transcript_id: str, interpro_data: Dict[str, dict]) -> Optional[dict]:
    """Find best matching InterPro result for a transcript ID."""
    candidate_keys = [
        transcript_id,  # Exact match (with version)
        transcript_id.replace(".", ""),  # No dots
    ]

    for key in candidate_keys:
        result = interpro_data.get(key)
        if result:
            log.debug("Found InterPro data for %s using key %s", transcript_id, key)
            return result

    # Fuzzy fallback by transcript prefix (kept conservative to avoid false matches)
    base_key = transcript_id.split(".")[0]
    if len(base_key) >= 10:
        for available_key, result in interpro_data.items():
            suffix = available_key[len(base_key):]
            if available_key.startswith(base_key) and all(c in "0123456789." for c in suffix):
                log.debug("Found InterPro data for %s using fuzzy match %s", transcript_id, available_key)
                return result

    return None


def _get_transcript_cds_segments(exons: List[Dict]) -> Tuple[str, List[Tuple[int, int]], int, int]:
    """Extract strand, reading-order CDS segments, and global CDS bounds from transcript exons."""
    cds_segments_unsorted = [
        (exon["cds_start"], exon["cds_end"])
        for exon in exons
        if exon.get("cds_start") is not None
        and exon.get("cds_end") is not None
        and exon["cds_end"] > exon["cds_start"]
    ]

    strand = exons[0].get("strand", "+") if exons else "+"
    cds_segments = sorted(cds_segments_unsorted, reverse=(strand != "+"))
    cds_start = min(start for start, _ in cds_segments_unsorted)
    cds_end = max(end for _, end in cds_segments_unsorted)

    return strand, cds_segments, cds_start, cds_end


def _aa_ranges_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    """Return True when two inclusive AA ranges overlap."""
    return start_a <= end_b and end_a >= start_b


def load_interpro_results(interpro_dir: str) -> Dict[str, dict]:
    """Load all InterPro JSON files from a directory.
    
    Args:
        interpro_dir: Path to directory containing InterPro JSON files
                      (e.g., data/neuro_project/output/interpro_results/)
    
    Returns:
        Dictionary mapping transcript_id (without version) -> full JSON data
    """
    interpro_data = {}
    interpro_path = Path(interpro_dir)
    
    if not interpro_path.exists() or not interpro_path.is_dir():
        log.warning(f"InterPro directory not found: {interpro_dir}")
        return interpro_data
    
    for json_file in interpro_path.glob("*.json"):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            # Extract transcript ID from filename (e.g., "G65356.7_ENST00000304501.2.json")
            # We store by the base ID without version for flexibility
            transcript_id = json_file.stem  # Remove .json
            
            interpro_data[transcript_id] = data
            log.debug(f"Loaded InterPro results for {transcript_id}")
        except json.JSONDecodeError as e:
            log.error(f"Failed to parse JSON {json_file}: {e}")
        except Exception as e:
            log.error(f"Error loading {json_file}: {e}")
    
    return interpro_data


def convert_protein_to_genomic(aa_start: int, aa_end: int, 
                               cds_start: int, cds_end: int,
                               strand: str = "+",
                               cds_segments: List[Tuple[int, int]] = None) -> Tuple[int, int]:
    """Convert protein (amino acid) coordinates to genomic coordinates.
    
    For single continuous CDS: genomic_bp = cds_start + (aa_position - 1) * 3
    For multi-exon CDS: converts using cumulative CDS position across segments
    
    Args:
        aa_start: Start position in protein (1-based, inclusive)
        aa_end: End position in protein (1-based, inclusive)
        cds_start: Genomic CDS start coordinate (used as reference only if cds_segments not provided)
        cds_end: Genomic CDS end coordinate (used as reference only if cds_segments not provided)
        strand: '+' or '-' (for forward/reverse strand)
        cds_segments: List of (start, end) tuples for CDS segments in logical order
                     (left-to-right for +strand, right-to-left for -strand)
    
    Returns:
        Tuple of (genomic_start, genomic_end) both inclusive, 1-based
    """
    if cds_segments:
        # Multi-/single-exon CDS with explicit segments.
        # Map amino-acid codon nucleotide offsets across CDS segments in reading order.
        # This is strand-aware and correctly handles codons split across exon junctions.
        total_cds_nt = sum((seg_end - seg_start + 1) for seg_start, seg_end in cds_segments)
        if total_cds_nt <= 0:
            return (cds_start, cds_end)

        def _map_nt_offset_to_genomic(nt_offset: int) -> int:
            # Clamp to valid CDS offset range
            remaining = max(0, min(int(nt_offset), total_cds_nt - 1))

            for seg_start, seg_end in cds_segments:
                seg_len = seg_end - seg_start + 1
                if remaining < seg_len:
                    if strand == "+":
                        return seg_start + remaining
                    # On reverse strand, translation proceeds from seg_end towards seg_start
                    return seg_end - remaining
                remaining -= seg_len

            # Fallback (shouldn't normally happen due clamping)
            last_start, last_end = cds_segments[-1]
            return last_end if strand == "+" else last_start

        start_nt_offset = (aa_start - 1) * 3
        end_nt_offset = (aa_end - 1) * 3 + 2

        g_start_raw = _map_nt_offset_to_genomic(start_nt_offset)
        g_end_raw = _map_nt_offset_to_genomic(end_nt_offset)

        return (min(g_start_raw, g_end_raw), max(g_start_raw, g_end_raw))
    
    elif strand == "+":
        # Forward strand: AA position directly translates to base position
        genomic_start = cds_start + (aa_start - 1) * 3
        genomic_end = cds_start + (aa_end - 1) * 3 + 2  # +2 to include all 3 bases of last codon
        return (genomic_start, genomic_end)
    else:
        # Reverse strand: positions are reversed
        # The first AA in the protein corresponds to the last codon genomically
        genomic_end = cds_end - (aa_start - 1) * 3 - 2
        genomic_start = cds_end - (aa_end - 1) * 3
        return (genomic_start, genomic_end)


def extract_domains_from_result(interpro_result: dict, 
                                cds_start: int, cds_end: int,
                                strand: str = "+",
                                min_evalue: float = 1e-5,
                                domain_types: list = None,
                                cds_segments: List[Tuple[int, int]] = None) -> List[Dict]:
    """Extract domain locations from a single InterPro result.
    
    Converts protein coordinates from the InterPro JSON to genomic coordinates
    and filters by e-value significance and domain type.
    
    Args:
        interpro_result: Single entry from InterPro JSON results list
        cds_start: Genomic CDS start coordinate for this transcript
        cds_end: Genomic CDS end coordinate for this transcript
        strand: '+' or '-'
        min_evalue: Minimum e-value to include (lower = more stringent)
        domain_types: List of domain types to include (e.g., ['DOMAIN']).
                     If None, defaults to ['DOMAIN'] only.
        cds_segments: List of (start, end) tuples for CDS segments in reading order
    
    Returns:
        List of domain dictionaries with keys:
        - name: Domain name
        - accession: Domain accession ID
        - type: Domain type (DOMAIN, FAMILY, REGION, SUPERFAMILY, etc.)
        - genomic_start: Converted genomic coordinate
        - genomic_end: Converted genomic coordinate
        - aa_start: Original amino acid start
        - aa_end: Original amino acid end
        - evalue: E-value score
        - score: InterPro score
        - library: Source library (PFAM, SMART, CDD, etc.)
    """
    if domain_types is None:
        domain_types = ['DOMAIN']  # Only include DOMAIN type by default
    
    domains = []
    
    if "matches" not in interpro_result:
        return domains
    
    for match in interpro_result["matches"]:
        try:
            signature = match.get("signature") or {}
            entry = signature.get("entry") or {}
            evalue = match.get("evalue")
            score = match.get("score")
            domain_type = entry.get("type") or signature.get("type", "UNKNOWN")

            if domain_type not in domain_types:
                log.debug("Skipping domain type %s (not in %s)", domain_type, domain_types)
                continue

            if evalue is not None and evalue > min_evalue:
                log.debug("Skipping domain (e-value=%s > %s)", evalue, min_evalue)
                continue

            domain_name = entry.get("name") or signature.get("name", "Unknown")
            domain_description = entry.get("description") or signature.get("description") or ""
            accession = signature.get("accession", "")
            interpro_id = entry.get("accession") or ""
            if not interpro_id and isinstance(accession, str) and accession.startswith("IPR"):
                interpro_id = accession
            library = signature.get("signatureLibraryRelease", {}).get("library", "")

            for location in match.get("locations", []):
                aa_start = location.get("start")
                aa_end = location.get("end")
                if aa_start is None or aa_end is None:
                    continue

                try:
                    genomic_start, genomic_end = convert_protein_to_genomic(
                        aa_start,
                        aa_end,
                        cds_start,
                        cds_end,
                        strand,
                        cds_segments,
                    )
                except Exception as e:
                    log.warning("Error converting coordinates for %s: %s", domain_name, e)
                    continue

                domains.append(
                    {
                        "name": domain_name,
                        "accession": accession,
                        "interpro_id": interpro_id,
                        "description": domain_description,
                        "type": domain_type,
                        "library": library,
                        "genomic_start": min(genomic_start, genomic_end),
                        "genomic_end": max(genomic_start, genomic_end),
                        "aa_start": aa_start,
                        "aa_end": aa_end,
                        "evalue": evalue,
                        "score": score,
                        "is_representative": location.get("representative", False),
                    }
                )
        except Exception as e:
            log.warning("Error processing match: %s", e)
    
    return domains


def get_domains_for_transcript(transcript_id: str,
                               interpro_data: Dict[str, dict],
                               isoforms_by_gene: Dict,
                               gene_id: str,
                               min_evalue: float = 1e-5) -> List[Dict]:
    """Get all domains for a specific transcript.
    
    Args:
        transcript_id: Transcript ID (e.g., "ENST00000304501.2")
        interpro_data: Output from load_interpro_results()
        isoforms_by_gene: GTF exon data mapping gene_id -> transcript_id -> exons
        gene_id: Gene ID (for lookup in isoforms_by_gene)
        min_evalue: Minimum e-value threshold
    
    Returns:
        List of domain dictionaries with genomic coordinates
    """
    domains = []

    interpro_result = _find_interpro_result(transcript_id, interpro_data)
    if not interpro_result or "results" not in interpro_result:
        log.debug(f"No InterPro results found for {transcript_id}")
        return domains
    
    # Get exons for this transcript
    if gene_id not in isoforms_by_gene or transcript_id not in isoforms_by_gene[gene_id]:
        log.warning(f"No exon data found for {transcript_id} in gene {gene_id}")
        return domains
    
    exons = isoforms_by_gene[gene_id][transcript_id]
    
    try:
        strand, cds_segments, cds_start, cds_end = _get_transcript_cds_segments(exons)
    except ValueError:
        log.warning(f"No CDS exons found for {transcript_id}")
        return domains
    
    # Extract domains from all results (only DOMAIN type by default)
    for result in interpro_result.get("results", []):
        result_domains = extract_domains_from_result(
            result, cds_start, cds_end, strand, min_evalue, 
            domain_types=['DOMAIN'],
            cds_segments=cds_segments
        )
        domains.extend(result_domains)
    
    # Remove duplicate/redundant domains
    domains = deduplicate_domains(domains)
    
    return domains




def deduplicate_domains(domains: List[Dict]) -> List[Dict]:
    """Remove redundant/overlapping domains, keeping the most significant ones.
    
    Groups domains by overlapping AA positions, keeping only the most significant
    in each group. Prefers manually-curated domains (evalue=None) over computed ones.
    
    Args:
        domains: List of domain dictionaries
        
    Returns:
        Deduplicated list of domains
    """
    if len(domains) <= 1:
        return domains
    
    # Sort by priority: curated first (evalue=None), then by e-value
    sorted_domains = sorted(
        domains,
        key=lambda d: (
            0 if d.get('interpro_id') else 1,  # Prefer entries with InterPro IDs
            0 if d['evalue'] is None else 1,  # Curated first
            d['evalue'] if d['evalue'] is not None else 0,  # Lower e-value better
            -(d['aa_end'] - d['aa_start'])  # Longer is better
        )
    )
    
    # Greedily select non-overlapping domains
    result = []
    for domain in sorted_domains:
        d_start = domain['aa_start']
        d_end = domain['aa_end']
        
        overlaps = any(
            _aa_ranges_overlap(d_start, d_end, selected['aa_start'], selected['aa_end'])
            for selected in result
        )
        
        # If no overlap, add it
        if not overlaps:
            result.append(domain)
    
    return result


def build_domain_mapping(interpro_dir: str,
                         isoforms_by_gene: Dict,
                         transcript_to_gene: Dict = None,
                         min_evalue: float = 1e-5) -> Dict[str, List[Dict]]:
    """Build a complete mapping of transcript_id -> domains across all transcripts.
    
    Args:
        interpro_dir: Path to InterPro results directory
        isoforms_by_gene: GTF exon data
        transcript_to_gene: Optional mapping of transcript_id -> gene_id
                           If not provided, will infer from isoforms_by_gene
        min_evalue: E-value threshold
    
    Returns:
        Dictionary mapping transcript_id -> list of domain dicts
    """
    interpro_data = load_interpro_results(interpro_dir)
    
    if not interpro_data:
        log.warning("No InterPro data loaded")
        return {}
    
    # Build transcript -> gene mapping if not provided
    if transcript_to_gene is None:
        transcript_to_gene = {}
        for gene_id, transcripts in isoforms_by_gene.items():
            for transcript_id in transcripts.keys():
                transcript_to_gene[transcript_id] = gene_id
    
    domain_mapping = {}
    for transcript_id, gene_id in transcript_to_gene.items():
        try:
            domains = get_domains_for_transcript(
                transcript_id, interpro_data, isoforms_by_gene,
                gene_id, min_evalue
            )
            if domains:
                domain_mapping[transcript_id] = domains
                log.debug(f"Found {len(domains)} domains for {transcript_id}")
        except Exception as e:
            log.warning(f"Error processing domains for {transcript_id}: {e}")
    
    return domain_mapping


def filter_overlapping_domains(domains: List[Dict], 
                               exons: List[Dict]) -> List[Dict]:
    """Filter domains to only those that overlap with exons.
    
    Args:
        domains: List of domain dictionaries with genomic_start/genomic_end
        exons: List of exon dictionaries with exon_start/exon_end
    
    Returns:
        Filtered list of domains that overlap with at least one exon
    """
    if not exons:
        return domains
    
    # Collect all exon regions
    exon_regions = [(e["exon_start"], e["exon_end"]) for e in exons]
    
    filtered = []
    for domain in domains:
        d_start = domain["genomic_start"]
        d_end = domain["genomic_end"]
        
        # Check if domain overlaps with any exon
        for e_start, e_end in exon_regions:
            if d_start <= e_end and d_end >= e_start:
                filtered.append(domain)
                break
    
    return filtered
