#!/usr/bin/env python3
"""
Script to submit amino acid sequences to InterPro Scan REST API
and retrieve protein domain data.

Usage:
    python interpro_scan.py <transcript_id>
    python interpro_scan.py --batch --fasta <fasta_file>
    python interpro_scan.py --gene-id <gene_id> --gtf <gtf_file>
    
Example:
    python interpro_scan.py G10110.21.nnc
    python interpro_scan.py --batch --fasta data/neuro_project/expressed_isoforms_PEP.fasta
    python interpro_scan.py --batch --fasta data/neuro_project/FASTA/proteins.fasta --output-dir results
    python interpro_scan.py --gene-id ENSG00000162426.16 --gtf data/neuro_project/expressed_isoforms.gtf
"""

import argparse
import requests
import time
import sys
import json
import re
import pandas as pd
from pathlib import Path
from collections import defaultdict

from utilities.gene_expression_analysis import load_gtf_mapping, map_transcripts_to_genes


class InterProScanAPI:
    """Client for InterPro Scan REST API"""
    
    BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
    
    def __init__(self, email="bm708@cam.ac.uk"):
        """
        Initialize the InterPro Scan API client.
        
        Args:
            email: Contact email (required by EBI)
        """
        self.email = email
        
    def submit_job(self, sequence, title="protein_sequence"):
        """
        Submit a protein sequence for InterPro Scan analysis.
        
        Args:
            sequence: Amino acid sequence (or multi-line sequences)
            title: Job title/identifier
            
        Returns:
            job_id: Job identifier for tracking
        """
        url = f"{self.BASE_URL}/run"
        
        params = {
            'email': self.email,
            'title': title,
            'goterms': 'false',  # Include GO terms
            'pathways': 'false',  # Include pathway annotations
            'sequence': sequence
        }
        
        print(f"Submitting sequence to InterPro Scan...")
        # Count sequences (if multiline format)
        seq_count = sequence.count('\n') if '\n' in sequence else 1
        total_length = len(sequence.replace('\n', '').replace('>', ''))
        print(f"Sequence length: {total_length} amino acids")
        
        try:
            response = requests.post(url, data=params)
            response.raise_for_status()
            job_id = response.text
            print(f"Job submitted successfully! Job ID: {job_id}")
            print(f"\nTrack your job in browser:")
            print(f"  Status: https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}")
            print(f"  Web Interface: https://www.ebi.ac.uk/interpro/result/InterProScan/{job_id}/")
            return job_id
        except requests.exceptions.RequestException as e:
            print(f"Error submitting job: {e}")
            sys.exit(1)
    
    def check_status(self, job_id):
        """
        Check the status of a submitted job.
        
        Args:
            job_id: Job identifier
            
        Returns:
            status: Job status (RUNNING, FINISHED, FAILED, etc.)
        """
        url = f"{self.BASE_URL}/status/{job_id}"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response.text
        except requests.exceptions.RequestException as e:
            print(f"Error checking job status: {e}")
            return "ERROR"
    
    def wait_for_completion(self, job_id, check_interval=10, max_wait=3600):
        """
        Wait for job to complete, checking periodically.
        
        Args:
            job_id: Job identifier
            check_interval: Seconds between status checks
            max_wait: Maximum seconds to wait
            
        Returns:
            success: True if job completed successfully
        """
        print(f"\nWaiting for job to complete...")
        elapsed = 0
        
        while elapsed < max_wait:
            status = self.check_status(job_id)
            print(f"Status: {status} (elapsed: {elapsed}s)")
            
            if status == "FINISHED":
                print("Job completed successfully!")
                return True
            elif status in ["FAILED", "ERROR", "NOT_FOUND"]:
                print(f"Job failed with status: {status}")
                return False
            elif status == "RUNNING":
                time.sleep(check_interval)
                elapsed += check_interval
            else:
                # Unknown status, wait and check again
                time.sleep(check_interval)
                elapsed += check_interval
        
        print(f"Job did not complete within {max_wait} seconds")
        return False
    
    def get_results(self, job_id, result_type='json'):
        """
        Retrieve results from a completed job.
        
        Args:
            job_id: Job identifier
            result_type: Format of results (json, tsv, xml, gff)
            
        Returns:
            results: Job results in requested format
        """
        url = f"{self.BASE_URL}/result/{job_id}/{result_type}"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            
            if result_type == 'json':
                return response.json()
            else:
                return response.text
        except requests.exceptions.RequestException as e:
            print(f"Error retrieving results: {e}")
            return None


def read_fasta_sequence(fasta_file, transcript_id):
    """
    Read a specific sequence from a FASTA file.
    
    Args:
        fasta_file: Path to FASTA file
        transcript_id: Transcript identifier to find
        
    Returns:
        sequence: Amino acid sequence (or None if not found)
    """
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('>'):
                # If we were building a sequence and found a match, return it
                if current_id == transcript_id and current_seq:
                    return ''.join(current_seq)
                
                # Start new sequence
                current_id = line[1:]  # Remove '>'
                current_seq = []
            else:
                # Add to current sequence
                if current_id == transcript_id:
                    current_seq.append(line)
    
    # Check if last sequence in file matches
    if current_id == transcript_id and current_seq:
        return ''.join(current_seq)
    
    return None


def read_all_fasta_sequences(fasta_file):
    """
    Read all sequences from a FASTA file.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        sequences: Dictionary mapping transcript IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id and current_seq:
                    sequences[current_id] = ''.join(current_seq)
                
                # Start new sequence
                current_id = line[1:]  # Remove '>'
                current_seq = []
            else:
                # Add to current sequence
                if current_id:
                    current_seq.append(line)
    
    # Save last sequence
    if current_id and current_seq:
        sequences[current_id] = ''.join(current_seq)
    
    return sequences


def parse_gtf_ids(attributes):
    """Extract transcript_id and gene_id from a GTF attributes field."""
    transcript_match = re.search(r'transcript_id\s+"([^"]+)"', attributes)
    gene_match = re.search(r'gene_id\s+"([^"]+)"', attributes)

    transcript_id = transcript_match.group(1) if transcript_match else None
    gene_id = gene_match.group(1) if gene_match else None
    return transcript_id, gene_id


def filter_transcripts_by_expression(expression_matrix_path, gtf_file, gene_id,
                                     cutoff_pct=1.5, exclude_fetal=True):
    """Return transcript IDs contributing at least cutoff_pct to gene expression."""
    if not expression_matrix_path or cutoff_pct <= 0:
        return None

    matrix_path = Path(expression_matrix_path)
    if not matrix_path.exists():
        print(f"Warning: expression matrix not found: {expression_matrix_path}")
        return None

    df_raw = pd.read_csv(matrix_path, sep="\t", index_col=0)
    if 'gene_id' not in df_raw.columns:
        transcript_to_gene = load_gtf_mapping(gtf_file)
        df_raw = map_transcripts_to_genes(df_raw, transcript_to_gene)

    sample_cols = [
        c for c in df_raw.columns
        if c != 'gene_id' and (not exclude_fetal or 'fetal' not in c.lower())
    ]
    if not sample_cols:
        print("Warning: no sample columns found for expression filtering")
        return None

    gene_block = df_raw[df_raw['gene_id'] == gene_id]
    if gene_block.empty:
        return []

    gene_transcripts = gene_block.drop(columns=['gene_id'])
    global_expr = gene_transcripts[sample_cols].sum(axis=1)
    total_gene_expr = global_expr.sum()
    if total_gene_expr == 0:
        return []

    global_expr_pct = (global_expr / total_gene_expr * 100)
    return global_expr_pct[global_expr_pct >= cutoff_pct].index.tolist()


def get_transcripts_for_gene(gtf_file, gene_id, expression_matrix_path=None,
                             cutoff_pct=1.5, exclude_fetal=True):
    """Return unique transcript IDs associated with a given gene ID from a GTF file.

    If expression_matrix_path is provided, only transcripts contributing at least
    cutoff_pct of the gene's total expression are returned.
    """
    transcripts = []
    seen = set()

    with open(gtf_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            attributes = fields[8]
            transcript_id, row_gene_id = parse_gtf_ids(attributes)

            if row_gene_id != gene_id or not transcript_id:
                continue

            # Prefer transcript features to avoid duplicates from exon/CDS rows.
            if feature_type == 'transcript' and transcript_id not in seen:
                seen.add(transcript_id)
                transcripts.append(transcript_id)

    if expression_matrix_path:
        filtered = filter_transcripts_by_expression(
            expression_matrix_path=expression_matrix_path,
            gtf_file=gtf_file,
            gene_id=gene_id,
            cutoff_pct=cutoff_pct,
            exclude_fetal=exclude_fetal,
        )
        if filtered is not None:
            filtered_set = set(filtered)
            transcripts = [t for t in transcripts if t in filtered_set]

    return transcripts


def format_results(results):
    """
    Format InterPro Scan results for display.
    
    Args:
        results: JSON results from InterPro Scan
    """
    if not results or 'results' not in results:
        print("No results found")
        return
    
    print("\n" + "="*80)
    print("INTERPRO SCAN RESULTS")
    print("="*80)
    
    for result in results['results']:
        if 'matches' not in result:
            continue
            
        print(f"\nSequence: {result.get('xref', [{}])[0].get('id', 'Unknown')}")
        print(f"Length: {result.get('seqlen', 'Unknown')} aa")
        
        matches = result['matches']
        print(f"\nFound {len(matches)} domain matches:")
        print("-"*80)
        
        for match in matches:
            signature = match.get('signature', {})
            print(f"\nDatabase: {signature.get('signatureLibraryRelease', {}).get('library', 'Unknown')}")
            print(f"Accession: {signature.get('accession', 'Unknown')}")
            print(f"Name: {signature.get('name', 'Unknown')}")
            print(f"Description: {signature.get('description', 'No description')}")
            
            if 'entry' in signature and signature['entry'] is not None:
                entry = signature['entry']
                print(f"InterPro Entry: {entry.get('accession', 'Unknown')}")
                print(f"InterPro Name: {entry.get('name', 'Unknown')}")
                print(f"Type: {entry.get('type', 'Unknown')}")
            
            # Print location information
            locations = match.get('locations', [])
            if locations:
                print(f"Locations:")
                for loc in locations:
                    start = loc.get('start', '?')
                    end = loc.get('end', '?')
                    score = loc.get('score', 'N/A')
                    print(f"  - Position: {start}-{end}, Score: {score}")


def save_results(results, output_file):
    """
    Save results to a JSON file.
    
    Args:
        results: JSON results from InterPro Scan
        output_file: Path to output file
    """
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_file}")


def process_batch_sequences(fasta_file, api, email='bm708@cam.ac.uk', 
                           output_dir='data/neuro_project/output/interpro_results',
                           check_interval=10, max_wait=3600, skip_existing=False):
    """
    Process all sequences from a FASTA file in batch.
    
    Args:
        fasta_file: Path to FASTA file
        api: InterProScanAPI instance
        email: Contact email for EBI services
        output_dir: Directory to save results
        check_interval: Seconds between status checks
        max_wait: Maximum seconds to wait per job
        skip_existing: Skip sequences that already have result files
    """
    # Read all sequences
    print(f"Reading sequences from: {fasta_file}")
    sequences = read_all_fasta_sequences(fasta_file)
    print(f"Found {len(sequences)} sequences\n")
    
    if not sequences:
        print("No sequences found in FASTA file")
        return
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Track jobs
    job_tracking = {}  # Maps job_id to transcript_id
    submitted_jobs = {}  # Maps transcript_id to job_id
    
    print("="*80)
    print("SUBMITTING JOBS")
    print("="*80)
    
    # Submit jobs for all sequences
    total_seqs = len(sequences)
    for idx, (transcript_id, sequence) in enumerate(sequences.items()):
        pct = int((idx / total_seqs) * 30)
        print(f"[PROGRESS_PERCENT] {pct}", flush=True)
        
        output_file = output_path / f"{transcript_id}.json"
        
        # Skip if already exists and skip_existing is True
        if skip_existing and output_file.exists():
            print(f"⊘ {transcript_id}: Already processed, skipping")
            continue
        
        # Clean sequence (remove stop codons)
        clean_seq = sequence.replace('*', '')
        
        if len(clean_seq) == 0:
            print(f"✗ {transcript_id}: Sequence is empty after cleaning")
            continue
        
        print(f"\n⊙ Submitting: {transcript_id}")
        print(f"  Length: {len(clean_seq)} aa")
        
        try:
            job_id = api.submit_job(clean_seq, title=transcript_id)
            job_tracking[job_id] = transcript_id
            submitted_jobs[transcript_id] = job_id
            print(f"  ✓ Job ID: {job_id}")
            
            # Small delay between submissions to avoid overwhelming the API
            time.sleep(2)
        except Exception as e:
            print(f"  ✗ Failed to submit: {e}")
            continue
            
    print(f"\n{'-'*80}")
    print(f"Total jobs submitted: {len(submitted_jobs)}")
    print(f"{'-'*80}\n")
    
    if not submitted_jobs:
        print("[PROGRESS_PERCENT] 100", flush=True)
        print("No jobs were submitted")
        return
        
    # Wait for all jobs to complete
    print("="*80)
    print("WAITING FOR RESULTS")
    print("="*80)
    
    completed_jobs = {}
    failed_jobs = {}
    
    total_jobs = len(submitted_jobs)
    for idx, (transcript_id, job_id) in enumerate(submitted_jobs.items()):
        pct = int(30 + (idx / total_jobs) * 50)
        print(f"[PROGRESS_PERCENT] {pct}", flush=True)
        
        print(f"\nWaiting for: {transcript_id} (Job ID: {job_id})")
        success = api.wait_for_completion(
            job_id,
            check_interval=check_interval,
            max_wait=max_wait
        )
        
        if success:
            completed_jobs[transcript_id] = job_id
        else:
            failed_jobs[transcript_id] = job_id
            
    # Retrieve and save results
    print(f"\n{'='*80}")
    print("RETRIEVING AND SAVING RESULTS")
    print(f"{'='*80}\n")
    
    successful_results = 0
    failed_results = 0
    
    total_completed = len(completed_jobs)
    for idx, (transcript_id, job_id) in enumerate(completed_jobs.items()):
        pct = int(80 + (idx / (total_completed or 1)) * 20)
        print(f"[PROGRESS_PERCENT] {pct}", flush=True)
        
        output_file = output_path / f"{transcript_id}.json"
        
        print(f"⊙ Retrieving results for: {transcript_id}")
        try:
            results = api.get_results(job_id, result_type='json')
            if results:
                save_results(results, str(output_file))
                successful_results += 1
                print(f"  ✓ Saved to: {output_file}")
            else:
                failed_results += 1
                print(f"  ✗ Failed to retrieve results")
        except Exception as e:
            failed_results += 1
            print(f"  ✗ Error: {e}")
            
    print("[PROGRESS_PERCENT] 100", flush=True)
    
    # Print summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"Total sequences: {len(sequences)}")
    print(f"Jobs submitted: {len(submitted_jobs)}")
    print(f"Jobs completed: {len(completed_jobs)}")
    print(f"Jobs failed: {len(failed_jobs)}")
    print(f"Results saved: {successful_results}")
    print(f"Results failed: {failed_results}")
    print(f"\nResults directory: {output_path.absolute()}")
    
    if failed_jobs:
        print(f"\nFailed jobs (may retry manually):")
        for transcript_id, job_id in failed_jobs.items():
            print(f"  - {transcript_id}: {job_id}")


def process_gene_sequences(gene_id, gtf_file, fasta_file, api,
                          output_dir='data/neuro_project/output/interpro_results',
                          check_interval=10, max_wait=3600, skip_existing=False,
                          expression_matrix_path='data/neuro_project/expressed_isoforms_matrix.txt',
                          cutoff_pct=1.5, exclude_fetal=True):
    """Process all transcripts for a given gene_id and save InterPro JSON results."""
    gtf_path = Path(gtf_file)
    if not gtf_path.exists():
        print(f"Error: GTF file not found: {gtf_file}")
        return

    print(f"Finding transcripts for gene: {gene_id}")
    print(f"From GTF: {gtf_file}")
    transcript_ids = get_transcripts_for_gene(
        gtf_file,
        gene_id,
        expression_matrix_path=expression_matrix_path,
        cutoff_pct=cutoff_pct,
        exclude_fetal=exclude_fetal,
    )

    if not transcript_ids:
        print(f"No transcripts found for gene_id '{gene_id}'")
        return

    print(f"Found {len(transcript_ids)} transcript(s) in GTF")
    print(f"Reading sequences from: {fasta_file}")
    all_sequences = read_all_fasta_sequences(fasta_file)

    if not all_sequences:
        print("No sequences found in FASTA file")
        return

    selected_sequences = {}
    missing_in_fasta = []

    for transcript_id in transcript_ids:
        sequence = all_sequences.get(transcript_id)
        if sequence is None:
            missing_in_fasta.append(transcript_id)
            continue
        selected_sequences[transcript_id] = sequence

    print(f"Transcripts with FASTA sequences: {len(selected_sequences)}")
    if missing_in_fasta:
        print(f"Transcripts missing in FASTA: {len(missing_in_fasta)}")

    if not selected_sequences:
        print("No transcript sequences were found in FASTA for the requested gene")
        return

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print("="*80)
    print("SUBMITTING GENE-LEVEL JOBS")
    print("="*80)

    submitted_jobs = {}
    for transcript_id, sequence in selected_sequences.items():
        output_file = output_path / f"{transcript_id}.json"

        if skip_existing and output_file.exists():
            print(f"⊘ {transcript_id}: Existing JSON found, skipping")
            continue

        clean_seq = sequence.replace('*', '')
        if len(clean_seq) == 0:
            print(f"✗ {transcript_id}: Sequence is empty after cleaning")
            continue

        print(f"\n⊙ Submitting: {transcript_id}")
        print(f"  Length: {len(clean_seq)} aa")

        try:
            job_id = api.submit_job(clean_seq, title=transcript_id)
            submitted_jobs[transcript_id] = job_id
            print(f"  ✓ Job ID: {job_id}")
            time.sleep(2)
        except Exception as e:
            print(f"  ✗ Failed to submit: {e}")

    if not submitted_jobs:
        print("No jobs were submitted")
        return

    print(f"\n{'='*80}")
    print("WAITING FOR GENE-LEVEL RESULTS")
    print(f"{'='*80}")

    completed_jobs = {}
    failed_jobs = {}

    for transcript_id, job_id in submitted_jobs.items():
        print(f"\nWaiting for: {transcript_id} (Job ID: {job_id})")
        success = api.wait_for_completion(
            job_id,
            check_interval=check_interval,
            max_wait=max_wait,
        )

        if success:
            completed_jobs[transcript_id] = job_id
        else:
            failed_jobs[transcript_id] = job_id

    print(f"\n{'='*80}")
    print("RETRIEVING AND SAVING GENE-LEVEL RESULTS")
    print(f"{'='*80}\n")

    successful_results = 0
    failed_results = 0

    for transcript_id, job_id in completed_jobs.items():
        output_file = output_path / f"{transcript_id}.json"
        print(f"⊙ Retrieving results for: {transcript_id}")
        try:
            results = api.get_results(job_id, result_type='json')
            if results:
                save_results(results, str(output_file))
                successful_results += 1
                print(f"  ✓ Saved to: {output_file}")
            else:
                failed_results += 1
                print(f"  ✗ Failed to retrieve results")
        except Exception as e:
            failed_results += 1
            print(f"  ✗ Error: {e}")

    print(f"\n{'='*80}")
    print("GENE-LEVEL SUMMARY")
    print(f"{'='*80}")
    print(f"Gene ID: {gene_id}")
    print(f"Transcripts in GTF: {len(transcript_ids)}")
    print(f"Transcripts in FASTA: {len(selected_sequences)}")
    print(f"Jobs submitted: {len(submitted_jobs)}")
    print(f"Jobs completed: {len(completed_jobs)}")
    print(f"Jobs failed: {len(failed_jobs)}")
    print(f"Results saved: {successful_results}")
    print(f"Results failed: {failed_results}")
    print(f"\nResults directory: {output_path.absolute()}")

    if missing_in_fasta:
        print(f"\nTranscripts missing from FASTA:")
        for transcript_id in missing_in_fasta:
            print(f"  - {transcript_id}")

    if failed_jobs:
        print(f"\nFailed jobs:")
        for transcript_id, job_id in failed_jobs.items():
            print(f"  - {transcript_id}: {job_id}")



def main():
    parser = argparse.ArgumentParser(
        description='Submit protein sequences to InterPro Scan REST API',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single sequence
  python interpro_scan.py G10110.21.nnc
  python interpro_scan.py G10118.1_ENST00000368648.8 --email your@email.com
  
  # Process all sequences in a FASTA file (batch mode)
  python interpro_scan.py --batch --fasta data/neuro_project/FASTA/proteins.fasta
  python interpro_scan.py --batch --fasta data/neuro_project/expressed_isoforms_PEP.fasta --output-dir custom_results
  python interpro_scan.py --batch --fasta data/neuro_project/FASTA/proteins.fasta --skip-existing

    # Process all transcripts for a given gene_id
    python interpro_scan.py --gene-id ENSG00000162426.16
    python interpro_scan.py --gene-id ENSG00000162426.16 --gtf data/neuro_project/expressed_isoforms.gtf --fasta data/neuro_project/FASTA/proteins.fasta
        """
    )
    
    parser.add_argument(
        'transcript_id',
        nargs='?',
        help='Transcript ID from the FASTA file (required for single-sequence mode)'
    )
    
    parser.add_argument(
        '--batch',
        action='store_true',
        help='Process all sequences in FASTA file (batch mode)'
    )

    parser.add_argument(
        '--gene-id',
        help='Gene ID to process all associated transcripts from GTF'
    )

    parser.add_argument(
        '--gtf',
        default='data/neuro_project/expressed_isoforms.gtf',
        help='Path to GTF file for gene_id mapping (default: data/neuro_project/expressed_isoforms.gtf)'
    )
    
    parser.add_argument(
        '--fasta',
        default='data/neuro_project/FASTA/proteins.fasta',
        help='Path to FASTA file (default: data/neuro_project/FASTA/proteins.fasta)'
    )

    parser.add_argument(
        '--expression-matrix',
        default='data/neuro_project/expressed_isoforms_matrix.txt',
        help='Isoform expression matrix for gene-level filtering '
             '(default: data/neuro_project/expressed_isoforms_matrix.txt)'
    )

    parser.add_argument(
        '--expression-cutoff',
        type=float,
        default=1.5,
        help='Minimum percent expression contribution per gene to keep a transcript '
             '(default: 1.5)'
    )
    
    parser.add_argument(
        '--email',
        default='bm708@cam.ac.uk',
        help='Contact email for EBI services (default: bm708@cam.ac.uk)'
    )
    
    parser.add_argument(
        '--output',
        help='Output file for results in single-sequence mode (JSON format)'
    )
    
    parser.add_argument(
        '--output-dir',
        default='data/neuro_project/output/interpro_results',
        help='Output directory for batch results (default: data/neuro_project/output/interpro_results)'
    )

    parser.add_argument(
        '--gene-output-dir',
        default='data/neuro_project/output/interpro_results',
        help='Output directory for gene_id mode JSON results (default: data/neuro_project/output/interpro_results)'
    )
    
    parser.add_argument(
        '--check-interval',
        type=int,
        default=10,
        help='Seconds between status checks (default: 10)'
    )
    
    parser.add_argument(
        '--max-wait',
        type=int,
        default=3600,
        help='Maximum seconds to wait for job completion (default: 3600)'
    )
    
    parser.add_argument(
        '--skip-existing',
        action='store_true',
        help='Skip sequences that already have result files (batch mode only)'
    )
    
    args = parser.parse_args()
    
    # Initialize API client
    api = InterProScanAPI(email=args.email)
    
    # BATCH MODE
    if args.batch:
        print(f"BATCH MODE: Processing all sequences in FASTA file")
        process_batch_sequences(
            fasta_file=args.fasta,
            api=api,
            email=args.email,
            output_dir=args.output_dir,
            check_interval=args.check_interval,
            max_wait=args.max_wait,
            skip_existing=args.skip_existing
        )

    # GENE MODE
    elif args.gene_id:
        print(f"GENE MODE: Processing transcripts for gene {args.gene_id}")
        process_gene_sequences(
            gene_id=args.gene_id,
            gtf_file=args.gtf,
            fasta_file=args.fasta,
            api=api,
            output_dir=args.gene_output_dir,
            check_interval=args.check_interval,
            max_wait=args.max_wait,
            skip_existing=args.skip_existing,
            expression_matrix_path=args.expression_matrix,
            cutoff_pct=args.expression_cutoff,
        )
    
    # SINGLE SEQUENCE MODE
    else:
        if not args.transcript_id:
            parser.print_help()
            print("\nError: transcript_id is required for single-sequence mode")
            print("Use --batch or --gene-id for alternative modes")
            sys.exit(1)
        
        # Read sequence from FASTA file
        print(f"Reading sequence for transcript: {args.transcript_id}")
        print(f"From file: {args.fasta}")
        
        sequence = read_fasta_sequence(args.fasta, args.transcript_id)
        
        if sequence is None:
            print(f"Error: Transcript ID '{args.transcript_id}' not found in {args.fasta}")
            sys.exit(1)
        
        # Remove any stop codons (*)
        sequence = sequence.replace('*', '')
        
        print(f"Found sequence with {len(sequence)} amino acids")
        
        # Submit job
        job_id = api.submit_job(sequence, title=args.transcript_id)
        
        # Wait for completion
        success = api.wait_for_completion(
            job_id, 
            check_interval=args.check_interval,
            max_wait=args.max_wait
        )
        
        if not success:
            print("Failed to complete InterPro Scan analysis")
            sys.exit(1)
        
        # Get results
        print("\nRetrieving results...")
        results = api.get_results(job_id, result_type='json')
        
        if results:
            # Display results
            format_results(results)
            
            # Save to file if requested
            if args.output:
                save_results(results, args.output)
            else:
                # Default output file in output directory
                output_dir = Path(args.output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)
                default_output = output_dir / f"{args.transcript_id}.json"
                save_results(results, str(default_output))
        else:
            print("Failed to retrieve results")
            sys.exit(1)


if __name__ == "__main__":
    main()
