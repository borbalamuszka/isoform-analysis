#!/usr/bin/env python3
"""Test script for InterPro domain parsing and visualization.

This script tests the domain coordinate conversion with the provided sample JSON.
"""

import json
import sys
from pathlib import Path

# Add repo root to path
ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from isoform_dashboard.interpro_parser import (
    extract_domains_from_result,
    convert_protein_to_genomic,
    load_interpro_results,
    get_domains_for_transcript,
)

def test_coordinate_conversion():
    """Test protein to genomic coordinate conversion."""
    print("\n" + "="*80)
    print("TEST 1: Protein to Genomic Coordinate Conversion")
    print("="*80)
    
    # Example: CDS from 1000 to 2000 bp, 367 amino acids
    cds_start = 1000
    cds_end = 2000
    
    # Test case 1: Full protein (AA 1-367)
    print("\nTest case 1: Full protein (AA 1-367)")
    print(f"  CDS span: {cds_start}-{cds_end} bp")
    g_start, g_end = convert_protein_to_genomic(1, 367, cds_start, cds_end, "+")
    print(f"  AA 1-367 → genomic {g_start}-{g_end}")
    assert g_start == cds_start, f"Expected {cds_start}, got {g_start}"
    expected_full_end = cds_start + (367 - 1) * 3 + 2
    assert g_end == expected_full_end, f"Expected {expected_full_end}, got {g_end}"
    print("  ✓ PASS")
    
    # Test case 2: Domain from AA 40-121 (from the JSON example)
    print("\nTest case 2: HMG box domain (AA 40-121)")
    g_start, g_end = convert_protein_to_genomic(40, 121, cds_start, cds_end, "+")
    print(f"  AA 40-121 → genomic {g_start}-{g_end}")
    # AA 40 starts at position (40-1)*3 = 117 bp from CDS start
    expected_start = cds_start + (40 - 1) * 3
    expected_end = cds_start + (121 - 1) * 3 + 2
    assert g_start == expected_start, f"Expected {expected_start}, got {g_start}"
    assert g_end == expected_end, f"Expected {expected_end}, got {g_end}"
    print(f"  Expected: {expected_start}-{expected_end}")
    print("  ✓ PASS")

    # Test case 3: Multi-exon conversion on reverse strand
    print("\nTest case 3: Multi-exon reverse-strand conversion")
    # CDS segments in reading order for '-' strand (higher genomic segment first)
    cds_segments_minus = [(300, 329), (200, 229)]  # 60 nt total = 20 AA

    g_start, g_end = convert_protein_to_genomic(
        1, 1,
        cds_start=200,
        cds_end=329,
        strand="-",
        cds_segments=cds_segments_minus,
    )
    print(f"  AA 1-1 (minus strand) → genomic {g_start}-{g_end}")
    assert (g_start, g_end) == (327, 329), f"Expected (327, 329), got ({g_start}, {g_end})"

    g_start, g_end = convert_protein_to_genomic(
        9, 12,
        cds_start=200,
        cds_end=329,
        strand="-",
        cds_segments=cds_segments_minus,
    )
    print(f"  AA 9-12 (crossing exon junction, minus strand) → genomic {g_start}-{g_end}")
    assert (g_start, g_end) == (224, 305), f"Expected (224, 305), got ({g_start}, {g_end})"
    print("  ✓ PASS")


def test_domain_extraction():
    """Test domain extraction from InterPro JSON."""
    print("\n" + "="*80)
    print("TEST 2: Domain Extraction from InterPro JSON")
    print("="*80)
    
    json_path = ROOT_DIR / "data/neuro_project/output/interpro_results/old/G65356.7_ENST00000304501.2.json"
    
    if not json_path.exists():
        print(f"\n⚠ WARNING: Test file not found at {json_path}")
        print("Skipping domain extraction test.")
        return
    
    print(f"\nLoading {json_path.name}...")
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Extract protein from JSON
    if "results" not in data or len(data["results"]) == 0:
        print("No results in JSON")
        return
    
    result = data["results"][0]
    print(f"✓ Loaded InterPro result")
    print(f"  Protein length: {len(result['sequence'])} AA")
    
    # Example CDS span: full protein is 388 AA
    # Total length = 388 * 3 = 1164 bp
    cds_start = 10000
    cds_end = cds_start + 388 * 3 - 1
    
    print(f"\nTest CDS span: {cds_start}-{cds_end} bp (for 388 AA protein)")
    
    # Extract domains
    domains = extract_domains_from_result(result, cds_start, cds_end, "+", min_evalue=1e-5)
    
    print(f"\n✓ Extracted {len(domains)} domains")
    
    # Print first few domains
    for i, domain in enumerate(domains[:5]):
        print(f"\n  Domain {i+1}:")
        print(f"    Name: {domain['name']}")
        print(f"    Type: {domain['type']}")
        print(f"    AA range: {domain['aa_start']}-{domain['aa_end']}")
        print(f"    Genomic range: {domain['genomic_start']}-{domain['genomic_end']}")
        print(f"    E-value: {domain['evalue']}")


def test_full_pipeline():
    """Test the full domain loading and filtering pipeline."""
    print("\n" + "="*80)
    print("TEST 3: Full Domain Loading Pipeline")
    print("="*80)
    
    # Create mock isoforms data
    isoforms_by_gene = {
        "G65356": {
            "G65356.7_ENST00000304501.2": [
                {
                    "exon_start": 10000,
                    "exon_end": 10500,
                    "cds_start": 10050,
                    "cds_end": 10450,
                    "strand": "+",
                },
                {
                    "exon_start": 11000,
                    "exon_end": 11500,
                    "cds_start": 11000,
                    "cds_end": 11164,  # partial
                    "strand": "+",
                }
            ]
        }
    }
    
    interpro_dir = ROOT_DIR / "data/neuro_project/output/interpro_results/old"
    
    if not interpro_dir.exists():
        print(f"\n⚠ WARNING: InterPro directory not found at {interpro_dir}")
        print("Skipping full pipeline test.")
        return
    
    print(f"Loading InterPro data from: {interpro_dir}")
    interpro_data = load_interpro_results(str(interpro_dir))
    print(f"✓ Loaded {len(interpro_data)} InterPro result files")
    
    # Try to get domains for our transcript
    transcript_id = "G65356.7_ENST00000304501.2"
    gene_id = "G65356"
    
    print(f"\nRequesting domains for {transcript_id}...")
    domains = get_domains_for_transcript(
        transcript_id, interpro_data, isoforms_by_gene, gene_id, min_evalue=1e-5
    )
    
    print(f"✓ Found {len(domains)} domains")
    
    if domains:
        print("\nFirst 3 domains:")
        for domain in domains[:3]:
            print(f"  - {domain['name']}: AA {domain['aa_start']}-{domain['aa_end']} → "
                  f"bp {domain['genomic_start']}-{domain['genomic_end']}")


if __name__ == "__main__":
    print("\n" + "="*80)
    print("InterPro Domain Parser Test Suite")
    print("="*80)
    
    try:
        test_coordinate_conversion()
        test_domain_extraction()
        test_full_pipeline()
        
        print("\n" + "="*80)
        print("✓ All tests completed successfully!")
        print("="*80 + "\n")
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
