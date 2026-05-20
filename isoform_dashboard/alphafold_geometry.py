"""AlphaFold 3D geometry handling.

Manages:
- Building transcript → geometry file mappings
- Looking up geometry by transcript ID
- Discovering pre-generated exon viewers
- Handling sequence-based geometry aliases
"""
import logging
import os
import glob
import re
from typing import Dict

log = logging.getLogger(__name__)


def build_alphafold_geometry_mapping(geometry_dir: str) -> dict:
    """Build a mapping from transcript_nodots -> geometry file paths.

    Output folders are named after the transcript ID only (dots removed,
    lowercased), e.g.:
        "G186352nnc"              -> "g186352nnc"
        "G1863546_ENST000002190223" -> "g1863546_enst000002190223"

    Args:
        geometry_dir: Path to the alphafold_geometry output directory produced
                      by extract_3d_geometry.py.

    Returns:
        Dict  transcript_nodots (lowercase)  ->
              {'folder': str, 'csv': str, 'html': str}
    """
    mapping = {}
    if not geometry_dir or not os.path.isdir(geometry_dir):
        log.warning("build_alphafold_geometry_mapping: geometry_dir %r is missing or not a directory", geometry_dir)
        return mapping

    all_entries = os.listdir(geometry_dir)
    subdirs = [e for e in all_entries if os.path.isdir(os.path.join(geometry_dir, e)) and not e.startswith('.')]
    log.info("build_alphafold_geometry_mapping: scanning %d subdirectories in %s", len(subdirs), geometry_dir)

    skipped_no_csv = []
    for folder_name in subdirs:
        folder_path = os.path.join(geometry_dir, folder_name)

        csv_files = sorted(
            f for f in os.listdir(folder_path)
            if f.startswith("geometry_") and f.endswith(".csv")
        )
        if not csv_files:
            skipped_no_csv.append(folder_name)
            continue

        csv_path  = os.path.join(folder_path, csv_files[0])
        html_name = csv_files[0].replace("geometry_", "viewer_").replace(".csv", ".html")
        html_path = os.path.join(folder_path, html_name)
        has_html  = os.path.isfile(html_path)

        # Key is the folder name exactly as produced by extract_3d_geometry.py
        # (transcript only, lowercase, dots removed).
        key = folder_name.lower()
        mapping[key] = {
            "folder": folder_path,
            "csv":    csv_path,
            "html":   html_path if has_html else None,
        }

    if skipped_no_csv:
        log.warning("build_alphafold_geometry_mapping: %d folder(s) skipped (no geometry CSV): %s",
                    len(skipped_no_csv), skipped_no_csv)

    log.info("build_alphafold_geometry_mapping: built %d entries (keys: %s%s)",
             len(mapping),
             ", ".join(list(mapping.keys())[:5]),
             " ..." if len(mapping) > 5 else "")
    return mapping


def discover_exon_viewers(geometry_mapping: Dict[str, Dict]) -> int:
    """Discover and load per-exon HTML viewers from geometry folders.
    
    Scans each geometry folder for files matching viewer_model{N}_exon{i}.html
    and builds sparse lists so exon_htmls[i] is the path for exon i.
    
    Args:
        geometry_mapping: Dict built by build_alphafold_geometry_mapping(), modified in-place
    
    Returns:
        Total number of exon HTML files found
    """
    total_found = 0
    
    for geo in geometry_mapping.values():
        folder = geo.get("folder", "")
        if not folder or not os.path.isdir(folder):
            geo["exon_htmls"] = []
            continue
        
        # Extract base filename from viewer_modelN.html
        base_html = geo.get("html", "")
        base_stem = os.path.splitext(os.path.basename(base_html))[0] if base_html else "viewer_model0"
        
        # Find all matching exon viewer files
        pattern = os.path.join(folder, f"{base_stem}_exon*.html")
        raw_found = glob.glob(pattern)
        
        if not raw_found:
            geo["exon_htmls"] = []
            continue
        
        # Map exon index → path
        exon_path_map = {}
        for p in raw_found:
            m = re.search(r'_exon(\d+)\.html$', p)
            if m:
                exon_path_map[int(m.group(1))] = p
        
        # Build sparse list (None for missing indices)
        if exon_path_map:
            max_idx = max(exon_path_map.keys())
            sparse = [exon_path_map.get(i) for i in range(max_idx + 1)]
            geo["exon_htmls"] = sparse
            total_found += len(exon_path_map)
        else:
            geo["exon_htmls"] = []
    
    if total_found:
        log.info("discover_exon_viewers: found %d pre-generated exon viewer HTML files", total_found)
    else:
        log.info("discover_exon_viewers: no exon viewer HTML files found "
                 "(run extract_3d_geometry.py --gtf ... to generate them)")
    
    return total_found


def resolve_alphafold_geometry(transcript_id: str,
                               geometry_mapping: dict) -> dict | None:
    """Look up the AlphaFold geometry entry for a transcript.

    Normalises the transcript ID the same way extract_3d_geometry.py does:
    strip dots, lowercase.

    Args:
        transcript_id: Transcript ID as it appears in the expression matrix,
                       e.g. "G18635.2.nnc", "G1863546_ENST00000219022.3"
        geometry_mapping: Output of build_alphafold_geometry_mapping(),
                          optionally extended by
                          extend_geometry_mapping_by_sequence().

    Returns:
        Dict with 'folder', 'csv', 'html' keys, or None if not found.
    """
    if not geometry_mapping:
        return None

    key = transcript_id.replace(".", "").lower()
    return geometry_mapping.get(key)


def extend_geometry_mapping_by_sequence(geometry_mapping: dict,
                                        protein_sequences: dict) -> int:
    """Add alias entries to geometry_mapping for transcripts that share an
    identical amino acid sequence with a transcript that already has geometry.

    For each pair of transcript IDs with the same protein sequence, if one has
    an entry in *geometry_mapping* but the other does not, a new entry is
    inserted pointing to the same geometry data.  The entry is marked with
    ``'is_sequence_alias': True`` and ``'alias_of': <original_key>`` so
    callers can distinguish real geometry from borrowed geometry.

    This function mutates *geometry_mapping* in-place and also returns the
    number of alias entries added.

    Args:
        geometry_mapping: Dict built by build_alphafold_geometry_mapping().
                          Modified in-place.
        protein_sequences: Dict mapping transcript_id -> amino-acid sequence,
                           as returned by dashboard_app.load_protein_sequences().

    Returns:
        Number of new alias entries added.
    """
    if not protein_sequences or not geometry_mapping:
        return 0

    # Normalise protein_sequences keys the same way as geometry_mapping keys
    # (dots stripped, lowercase) so the lookup is consistent.
    normalised: dict[str, list[str]] = {}   # norm_key -> list of original transcript IDs
    seq_of: dict[str, str] = {}             # norm_key -> amino-acid sequence

    for tid, seq in protein_sequences.items():
        norm_key = tid.replace(".", "").lower()
        normalised.setdefault(norm_key, []).append(tid)
        seq_of[norm_key] = seq

    # Group normalised keys by their amino acid sequence
    seq_to_keys: dict[str, list[str]] = {}
    for norm_key, seq in seq_of.items():
        seq_to_keys.setdefault(seq, []).append(norm_key)

    aliases_added = 0
    for seq, keys in seq_to_keys.items():
        if len(keys) < 2:
            continue  # Unique sequence — nothing to share

        # Find which of the group already have geometry
        with_geo = [k for k in keys if k in geometry_mapping]
        without_geo = [k for k in keys if k not in geometry_mapping]

        if not with_geo or not without_geo:
            continue  # Either all have geometry already, or none do

        # Use the first available geometry entry as the source
        source_key = with_geo[0]
        source_entry = geometry_mapping[source_key]

        for missing_key in without_geo:
            geometry_mapping[missing_key] = {
                **source_entry,
                "is_sequence_alias": True,
                "alias_of": source_key,
            }
            aliases_added += 1

    if aliases_added:
        log.info(
            "extend_geometry_mapping_by_sequence: added %d alias entries "
            "(%d transcripts now have geometry via sequence identity)",
            aliases_added,
            len(geometry_mapping),
        )
    return aliases_added
