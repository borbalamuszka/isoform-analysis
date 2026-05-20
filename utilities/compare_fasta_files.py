#!/usr/bin/env python3
"""
Compare FASTA files to identify matches and differences.

This script compares:
1. expressed_isoforms_ORF.fasta with cds_output.fa (nucleotide sequences)
2. expressed_isoforms_PEP.fasta with proteins.fa (protein sequences)

It reports:
- Number of items in each file
- Number of matching isoforms (by ID)
- Detailed sequence comparison for 20 matching isoforms
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple, Set
from collections import defaultdict


def parse_fasta(file_path: str) -> Dict[str, str]:
    """
    Parse a FASTA file and return a dictionary of ID -> sequence.
    
    Args:
        file_path: Path to the FASTA file
        
    Returns:
        Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    
                    # Start new sequence
                    current_id = line[1:]  # Remove '>' character
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Don't forget the last sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
                
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        sys.exit(1)
    
    return sequences


def analyze_sequence_similarity(seq1: str, seq2: str) -> dict:
    """
    Analyze the type of difference between two sequences.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        
    Returns:
        Dictionary with similarity analysis:
        - 'identical': True if sequences are exactly identical
        - 'case_only': True if sequences differ only in case
        - 'single_char': True if sequences differ by exactly one character
        - 'diff_count': Number of differences
    """
    result = {
        'identical': False,
        'case_only': False,
        'single_char': False,
        'diff_count': 0
    }
    
    # Check for exact identity
    if seq1 == seq2:
        result['identical'] = True
        return result
    
    # Check for case-only differences
    if seq1.upper() == seq2.upper():
        result['case_only'] = True
        # Count case differences
        result['diff_count'] = sum(1 for c1, c2 in zip(seq1, seq2) if c1 != c2)
        return result
    
    # Count all differences
    min_len = min(len(seq1), len(seq2))
    diff_count = sum(1 for i in range(min_len) if seq1[i] != seq2[i])
    
    # Add length difference to diff count
    diff_count += abs(len(seq1) - len(seq2))
    
    result['diff_count'] = diff_count
    result['single_char'] = (diff_count == 1)
    
    return result


def compare_sequences(seq1: str, seq2: str, max_display: int = 100) -> Tuple[bool, List[dict]]:
    """
    Compare two sequences and identify differences.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        max_display: Maximum number of differences to display
        
    Returns:
        Tuple of (sequences_identical, list_of_differences)
    """
    if seq1 == seq2:
        return True, []
    
    differences = []
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    
    # Compare position by position
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            differences.append({
                'position': i + 1,
                'file1': seq1[i],
                'file2': seq2[i],
                'type': 'substitution'
            })
            
            if len(differences) >= max_display:
                break
    
    # Check for length differences
    if len(seq1) != len(seq2):
        if len(seq1) > len(seq2):
            differences.append({
                'position': len(seq2) + 1,
                'file1': f'{len(seq1) - len(seq2)} extra characters',
                'file2': 'END',
                'type': 'length_diff'
            })
        else:
            differences.append({
                'position': len(seq1) + 1,
                'file1': 'END',
                'file2': f'{len(seq2) - len(seq1)} extra characters',
                'type': 'length_diff'
            })
    
    return False, differences


def highlight_differences(seq1: str, seq2: str, context: int = 20) -> List[str]:
    """
    Create highlighted regions showing differences between sequences.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        context: Number of characters to show around each difference
        
    Returns:
        List of formatted difference regions
    """
    if seq1 == seq2:
        return ["Sequences are identical"]
    
    regions = []
    min_len = min(len(seq1), len(seq2))
    
    # Find all difference positions
    diff_positions = []
    for i in range(min_len):
        if seq1[i] != seq2[i]:
            diff_positions.append(i)
    
    # Group nearby differences
    if diff_positions:
        grouped = []
        start = diff_positions[0]
        end = diff_positions[0]
        
        for pos in diff_positions[1:]:
            if pos - end <= context * 2:
                end = pos
            else:
                grouped.append((start, end))
                start = pos
                end = pos
        grouped.append((start, end))
        
        # Create highlighted regions
        for start, end in grouped[:10]:  # Limit to 10 regions
            region_start = max(0, start - context)
            region_end = min(min_len, end + context + 1)
            
            region1 = seq1[region_start:region_end]
            region2 = seq2[region_start:region_end]
            
            # Mark differences with brackets
            marked1 = []
            marked2 = []
            for i, (c1, c2) in enumerate(zip(region1, region2)):
                actual_pos = region_start + i
                if actual_pos >= start and actual_pos <= end and c1 != c2:
                    marked1.append(f'[{c1}]')
                    marked2.append(f'[{c2}]')
                else:
                    marked1.append(c1)
                    marked2.append(c2)
            
            region_text = f"  Position {start+1}-{end+1}:\n"
            region_text += f"    File1: {''.join(marked1)}\n"
            region_text += f"    File2: {''.join(marked2)}"
            regions.append(region_text)
    
    # Check length differences
    if len(seq1) != len(seq2):
        if len(seq1) > len(seq2):
            extra = seq1[len(seq2):]
            if len(extra) > 60:
                extra = extra[:60] + "..."
            regions.append(f"  File1 has {len(seq1) - len(seq2)} extra characters at end: {extra}")
        else:
            extra = seq2[len(seq1):]
            if len(extra) > 60:
                extra = extra[:60] + "..."
            regions.append(f"  File2 has {len(seq2) - len(seq1)} extra characters at end: {extra}")
    
    return regions


def compare_fasta_files(file1_path: str, file2_path: str, file_type: str) -> None:
    """
    Compare two FASTA files and print detailed comparison.
    
    Args:
        file1_path: Path to first FASTA file
        file2_path: Path to second FASTA file
        file_type: Type description (e.g., "ORF/CDS" or "Protein")
    """
    print(f"\n{'=' * 80}")
    print(f"Comparing {file_type} sequences")
    print(f"{'=' * 80}")
    print(f"File 1: {Path(file1_path).name}")
    print(f"File 2: {Path(file2_path).name}")
    print()
    
    # Parse both files
    print("Parsing FASTA files...")
    seqs1 = parse_fasta(file1_path)
    seqs2 = parse_fasta(file2_path)
    
    # Report counts
    print(f"\n{'-' * 80}")
    print(f"ITEM COUNTS:")
    print(f"  File 1 ({Path(file1_path).name}): {len(seqs1):,} sequences")
    print(f"  File 2 ({Path(file2_path).name}): {len(seqs2):,} sequences")
    print(f"{'-' * 80}")
    
    # Find matching IDs
    ids1 = set(seqs1.keys())
    ids2 = set(seqs2.keys())
    
    matching_ids = ids1 & ids2
    only_in_file1 = ids1 - ids2
    only_in_file2 = ids2 - ids1
    
    print(f"\nMATCHING ANALYSIS:")
    print(f"  Matching IDs (present in both files): {len(matching_ids):,}")
    print(f"  Only in File 1: {len(only_in_file1):,}")
    print(f"  Only in File 2: {len(only_in_file2):,}")
    print(f"  Match percentage: {len(matching_ids) / max(len(seqs1), len(seqs2)) * 100:.2f}%")
    
    # Show some examples of non-matching IDs
    if only_in_file1:
        print(f"\n  Examples only in File 1 (up to 10):")
        for id_name in sorted(only_in_file1)[:10]:
            print(f"    - {id_name}")
        if len(only_in_file1) > 10:
            print(f"    ... and {len(only_in_file1) - 10} more")
    
    if only_in_file2:
        print(f"\n  Examples only in File 2 (up to 10):")
        for id_name in sorted(only_in_file2)[:10]:
            print(f"    - {id_name}")
        if len(only_in_file2) > 10:
            print(f"    ... and {len(only_in_file2) - 10} more")
    
    # Compare sequences for matching IDs
    if matching_ids:
        print(f"\n{'-' * 80}")
        print(f"SEQUENCE COMPARISON (20 matching isoforms):")
        print(f"{'-' * 80}")
        
        # Track overall statistics for ALL matching sequences
        total_identical = 0
        total_case_only = 0
        total_single_char = 0
        total_other_diff = 0
        
        identical_count = 0
        different_count = 0
        
        # First, analyze ALL matching sequences for statistics
        print(f"\nAnalyzing all {len(matching_ids):,} matching sequences...")
        for seq_id in matching_ids:
            seq1 = seqs1[seq_id]
            seq2 = seqs2[seq_id]
            similarity = analyze_sequence_similarity(seq1, seq2)
            
            if similarity['identical']:
                total_identical += 1
            elif similarity['case_only']:
                total_case_only += 1
            elif similarity['single_char']:
                total_single_char += 1
            else:
                total_other_diff += 1
        
        # Print overall statistics
        print(f"\n{'-' * 80}")
        print(f"SEQUENCE SIMILARITY STATISTICS (All {len(matching_ids):,} matching transcripts):")
        print(f"{'-' * 80}")
        print(f"  Identical sequences:                    {total_identical:>8,} ({total_identical/len(matching_ids)*100:>6.2f}%)")
        print(f"  Differ only in case (upper/lower):      {total_case_only:>8,} ({total_case_only/len(matching_ids)*100:>6.2f}%)")
        print(f"  Differ by single character:             {total_single_char:>8,} ({total_single_char/len(matching_ids)*100:>6.2f}%)")
        print(f"  Other differences (multiple changes):   {total_other_diff:>8,} ({total_other_diff/len(matching_ids)*100:>6.2f}%)")
        print(f"{'-' * 80}")
        
        # Now show detailed comparison for first 20
        print(f"\nDETAILED COMPARISON (First 20 matching isoforms):")
        print(f"{'-' * 80}")
        
        # Sort matching IDs for consistent output
        sorted_matching = sorted(matching_ids)[:20]
        
        for idx, seq_id in enumerate(sorted_matching, 1):
            seq1 = seqs1[seq_id]
            seq2 = seqs2[seq_id]
            
            is_identical, differences = compare_sequences(seq1, seq2)
            similarity = analyze_sequence_similarity(seq1, seq2)
            
            print(f"\n{idx}. {seq_id}")
            print(f"   Length: File1={len(seq1):,} bp/aa, File2={len(seq2):,} bp/aa")
            
            if is_identical:
                identical_count += 1
                print(f"   Status: ✓ IDENTICAL")
            else:
                different_count += 1
                
                # Provide more specific status
                if similarity['case_only']:
                    print(f"   Status: ⚠ DIFFERS ONLY IN CASE ({similarity['diff_count']} positions)")
                elif similarity['single_char']:
                    print(f"   Status: ⚠ DIFFERS BY SINGLE CHARACTER")
                else:
                    print(f"   Status: ✗ DIFFERENT ({len([d for d in differences if d['type'] == 'substitution'])} substitutions)")
                
                # Show highlighted differences
                regions = highlight_differences(seq1, seq2)
                if regions[0] != "Sequences are identical":
                    print(f"   Differences:")
                    for region in regions:
                        print(region)
                
                # Show summary of differences
                if len(differences) > 0:
                    print(f"\n   Summary of first {min(10, len(differences))} differences:")
                    for diff in differences[:10]:
                        if diff['type'] == 'substitution':
                            print(f"     - Pos {diff['position']}: {diff['file1']} → {diff['file2']}")
                        else:
                            print(f"     - Pos {diff['position']}: Length difference - "
                                  f"{diff['file1']} vs {diff['file2']}")
                    
                    if len(differences) > 10:
                        print(f"     ... and {len(differences) - 10} more differences")
        
        # Summary statistics
        print(f"\n{'-' * 80}")
        print(f"SUMMARY OF 20 COMPARED SEQUENCES:")
        print(f"  Identical sequences: {identical_count}")
        print(f"  Different sequences: {different_count}")
        if identical_count + different_count > 0:
            print(f"  Identity rate: {identical_count / (identical_count + different_count) * 100:.1f}%")
        print(f"{'-' * 80}")


def main():
    """Main function to run the comparison."""
    
    # Define file paths
    base_dir = Path(__file__).parent.parent
    
    # ORF/CDS comparison
    orf_file = base_dir / "data/neuro_project/expressed_isoforms_ORF.fasta"
    cds_file = base_dir / "data/FASTA/cds_output.fa"
    
    # Protein comparison
    pep_file = base_dir / "data/neuro_project/expressed_isoforms_PEP.fasta"
    protein_file = base_dir / "data/FASTA/proteins.fa"
    
    # Output file
    output_dir = base_dir / "data/neuro_project/output"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "fasta_comparison_report.txt"
    
    # Redirect print to file only
    import sys
    from datetime import datetime
    
    # Show minimal console output
    print(f"Starting FASTA comparison analysis...")
    print(f"Output will be saved to: {output_file}")
    
    # Open output file and redirect stdout
    original_stdout = sys.stdout
    f = open(output_file, 'w')
    sys.stdout = f
    
    print("=" * 80)
    print("FASTA FILE COMPARISON ANALYSIS")
    print("=" * 80)
    print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Output file: {output_file}")
    print("\nThis script compares:")
    print("1. ORF sequences: expressed_isoforms_ORF.fasta vs cds_output.fa")
    print("2. Protein sequences: expressed_isoforms_PEP.fasta vs proteins.fa")
    print()
    
    # Check if files exist
    files_to_check = [
        (orf_file, "ORF file"),
        (cds_file, "CDS file"),
        (pep_file, "PEP file"),
        (protein_file, "Protein file")
    ]
    
    missing_files = []
    for file_path, description in files_to_check:
        if not file_path.exists():
            missing_files.append(f"{description}: {file_path}")
    
    if missing_files:
        sys.stdout = original_stdout
        print("ERROR: The following files are missing:")
        for missing in missing_files:
            print(f"  - {missing}")
        f.close()
        sys.exit(1)
    
    # Run comparisons
    try:
        compare_fasta_files(str(orf_file), str(cds_file), "ORF/CDS (Nucleotide)")
        compare_fasta_files(str(pep_file), str(protein_file), "Protein (Amino Acid)")
        
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE")
        print("=" * 80)
        print(f"\nResults written to: {output_file}")
    finally:
        # Restore stdout and close file
        sys.stdout = original_stdout
        f.close()
        print(f"Analysis complete! Results saved to: {output_file}")


if __name__ == "__main__":
    main()
