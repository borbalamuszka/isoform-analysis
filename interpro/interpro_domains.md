# InterPro Domain Visualization

## Overview

This feature visualizes protein domains from InterPro scan results overlaid on exon structures. Domains appear as colored rectangles below exons, showing where functional domains align with genomic features.

## Quick Start

### Run Dashboard with Domains

```bash
python isoform_dashboard/dashboard_app.py \
    --input-mean "data/distributions_condition_mean.tsv" \
    --input-sum "data/distributions_condition_sum.tsv" \
    --exons "data/expressed_isoforms.gtf" \
    --interpro "data/interpro_results.json" \
    --port 8050
```

### Process Proteins with InterProScan

```bash
# Single sequence
python interpro_scan.py G249.3.nnc --fasta data/neuro_project/FASTA/proteins.fasta

# Batch processing (entire FASTA file)
python interpro_scan.py --batch --fasta data/neuro_project/FASTA/proteins.fasta
```

## Files

### Core Implementation
- **`isoform_dashboard/interpro_parser.py`** - Domain parsing, coordinate conversion, filtering
- **`isoform_dashboard/interpro_scan.py`** - InterProScan API interaction (single & batch)

### Modified Files
- **`isoform_dashboard/visualizations.py`** - Domain rendering below exons
- **`isoform_dashboard/dashboard_app.py`** - CLI arguments & domain loading
- **`isoform_dashboard/app_layout.py`** - Callback integration

## How It Works

### 1. Generate InterPro Results

Use the optimized `interpro_scan.py` script:

```bash
# For entire proteome (recommended for batch)
python interpro_scan.py --batch --fasta data/neuro_project/FASTA/proteins.fasta \
    --skip-existing  # Resume capability

# Results saved to: data/neuro_project/output/interpro_results/
```

The script provides:
- ✅ Automatic job tracking and error handling
- ✅ Progress reporting with summary statistics
- ✅ Resume from where you left off
- ✅ Backward compatible single-sequence mode

### 2. Coordinate Conversion

InterPro provides amino acid (AA) positions. The system converts to genomic coordinates:

```
AA position N → genomic bp = CDS_start + (N - 1) × 3
```

This uses GTF CDS boundaries to align domains with exons.

### 3. Visualization

Domains display as:
- **Position**: Below exon bars (smaller height)
- **Color**: By type
  - Green: DOMAIN
  - Orange: FAMILY
  - Blue: REGION
  - Purple: SUPERFAMILY
- **Hover Info**: Domain name, AA range, type, accession

### 4. Overlap Handling

Multiple overlapping domains are merged using sweep-line algorithm, reducing clutter while showing the union of annotated regions.

## Configuration

### Dashboard Parameters

```bash
--interpro FILE                    # Path to InterPro JSON
--interpro-evalue-threshold VALUE  # Default: 1e-5 (typical significance)
```

### InterProScan Parameters

```bash
--batch                    # Process entire FASTA file
--skip-existing           # Resume (skip already-processed)
--max-wait SECONDS        # Timeout per job (default: 3600)
--email USER@EXAMPLE.COM  # Required for EBI service
--output-dir DIRECTORY    # Custom results location
```

## Data Sources

### Input Requirements

1. **InterPro JSON** - From InterProScan output (standard format)
2. **GTF File** - With CDS annotations for coordinate conversion
3. **Protein FASTA** (optional) - For transcript ID matching

### Output

Results saved to: `data/neuro_project/output/interpro_results`

One JSON file per sequence:
- `G249.3.nnc.json`
- `G249.4.nnc.json`
- etc.

## Troubleshooting

### Domains Not Appearing

1. Check file path exists:
   ```bash
   ls -la path/to/interpro_results.json
   ```

2. Verify JSON format:
   ```bash
   python -m json.tool path/to/interpro_results.json | head -50
   ```

3. Check evalue threshold:
   ```bash
   python -c "from isoform_dashboard.interpro_parser import parse_interpro_json; \
   d = parse_interpro_json('path/to/file.json', evalue_threshold=None); \
   print(f'Found {sum(len(v) for v in d.values())} domains')"
   ```

4. Verify GTF has CDS info:
   ```bash
   grep CDS data/expressed_isoforms.gtf | head -5
   ```

### Job Failures

If some sequences fail during batch processing:

```bash
# Check failed jobs in summary report
tail -30 interpro_batch.log | grep "Failed jobs"

# Retry failed sequences
python interpro_scan.py --batch --fasta data/neuro_project/FASTA/proteins.fasta \
    --skip-existing
```

### Performance

- Parsing: One-time at startup
- Large files: Use higher `--interpro-evalue-threshold` to reduce domains
- InterPro service: Bottleneck is EBI server response time, not the script

## Workflow Example

```bash
# 1. Test single sequence
python interpro_scan.py G249.3.nnc --fasta data/neuro_project/FASTA/proteins.fasta

# 2. Verify output
cat data/neuro_project/output/interpro_results/G249.3.nnc.json | head -50

# 3. Process all sequences (takes hours, can run in background)
python interpro_scan.py --batch --fasta data/neuro_project/FASTA/proteins.fasta > batch.log 2>&1 &

# 4. Monitor progress
tail -f batch.log

# 5. Run dashboard with results
python isoform_dashboard/dashboard_app.py \
    --input-mean "..." \
    --exons "expressed_isoforms.gtf" \
    --interpro "data/neuro_project/output/interpro_results/G2479.1.nnc.json" \
    --port 8050
```

## Performance Notes

- **Memory**: ~50-100MB for loading full proteome
- **Rate Limiting**: 2-second delay between submissions (EBI-friendly)
- **Processing**: Hours for ~289k sequences (depends on EBI service load)
- **Bottleneck**: InterPro Scan service response time

## References

- **InterProScan**: https://www.ebi.ac.uk/interpro/search/sequence/
- **InterPro API**: https://www.ebi.ac.uk/interpro/api/
- **GTF Format**: https://www.ensembl.org/info/website/upload/gff.html
