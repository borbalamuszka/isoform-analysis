# Isoform Analysis Tools & Interactive Dashboard

- Analyze gene isoform expression and generate distribution tables/plots.
- Compute per‑gene metrics (entropy, Spearman) and bootstrap isoform means with CIs.
- Query protein domains via InterPro Scan.
- Create AlphaFold 3D geometry visualisations.
- Explore results in an interactive Dash dashboard with exon, domain, and protein structure visualization.
- Works with generic datasets; grouping and labels depend on the input data provided.

## Quick start: run the dashboard locally

1. **Clone and enter the repo**
    ```bash
    git clone https://github.com/borbalamuszka/isoform-analysis.git
    cd isoform-analysis
    ```

2. **Create and activate a virtual environment (recommended)**

  - macOS / Linux (bash or zsh):
    ```bash
    python3 -m venv isoform_dashboard_env
    source isoform_dashboard_env/bin/activate
    ```

  - Windows (PowerShell or cmd):
    ```powershell
    python -m venv isoform_dashboard_env
    isoform_dashboard_env\Scripts\activate
    ```

3. **Install dependencies** (from the repo root)
    ```bash
    pip install -r requirements.txt
    ```

4. **Run the dashboard** (from the repo root)

   ```bash
   python -m isoform_dashboard.dashboard_app \
     --input-mean   data/project/output/isoform_distributions/distributions_mean.tsv \
     --input-sum    data/project/output/isoform_distributions/distributions_sum.tsv \
     --ci-file      data/project/output/isoform_distributions/confidence_intervals.tsv \
     --exons        data/project/expressed_isoforms.gtf \
     --proteins     data/project/proteins.fasta \
     --geometry-dir data/project/output/alphafold_geometry \
     --interpro-dir data/project/output/interpro_results
   ```

   Windows PowerShell (line continuation with backticks):

   ```powershell
   python -m isoform_dashboard.dashboard_app `
     --input-mean   data/project/output/isoform_distributions/distributions_mean.tsv `
     --input-sum    data/project/output/isoform_distributions/distributions_sum.tsv `
     --ci-file      data/project/output/isoform_distributions/confidence_intervals.tsv `
     --exons        data/project/expressed_isoforms.gtf `
     --proteins     data/project/proteins.fasta `
     --geometry-dir data/project/output/alphafold_geometry `
     --interpro-dir data/project/output/interpro_results
   ```

   For Windows cmd.exe, replace the backticks with carets (`^`) or paste the command on one line.

   Then open `http://127.0.0.1:8050` in your browser.

**Troubleshooting:**

- If port 8050 is busy: add `--port 8051` to the command.
- If Python not found: ensure Python 3.8+ is installed and in PATH.
- To stop the dashboard: press Ctrl+C in the terminal.
- Windows note: run from the repo root (paths are relative), and prefer forward slashes in paths (Python accepts them on Windows). If you copy a multi-line command, use PowerShell backticks (`) or cmd carets (^) for line continuation, or paste it as a single line.

> **Alternatives to cloning the repo?**
> The simplest supported workflow is cloning this GitHub repo. You could instead download a ZIP from GitHub ("Code → Download ZIP") and unpack it, but you still need the same Python environment, data files, and commands as above.

<!-- ## Dashboard preview

![Isoform Entropy Dashboard](assets/isoform-entropy-dashboard.jpg) -->

## Input data

- **Expression data file**
  - E.g. `expressed_isoforms_matrix.tsv` / `.txt`
  - Rows: `transcript_id`.
  - Columns: sample names (e.g. `Sample_A`, `Group1_Sample3`).
  - Used by distribution scripts and dashboard.

- **Positions / GTF file**
  - E.g. `expressed_isoforms.gtf`
  - Genomic positions and metadata (chrom, exon_start/end, cds_start/end, gene_id, transcript_id, …).
  - Used by distribution scripts and dashboard.

- **Optional confidence intervals file**
  - E.g. `confidence_intervals.tsv`
  - Bootstrap mean and CI bounds per isoform, plus optional group‑specific CIs (based on your input data).
  - Used by the dashboard to overlay error bars.

## Main components

### Isoform distributions: `isoform_distribution/distributions.py`

Generates isoform distribution tables (one row per retained transcript).

- **Features**
  - Filters low‑contribution isoforms by percentage cutoff.
  - Aggregation statistics:
    - `sum`: summed expression across samples.
    - `mean`: mean expression across samples.
    - `normalized`: per‑gene per‑sample normalization (isoforms sum to 1); global value = mean normalized contribution.
  - Metadata grouping (default): use a metadata file to map samples to groups.
- **Inputs**
  - `--matrix`: isoform expression matrix (rows = transcript_id, columns = samples).
  - `--gtf`: GTF file for transcript → gene mapping.
  - `--output-dir`: output directory for distribution tables.
  - `--cutoff-pct`: minimum percent contribution to retain an isoform (default: 1.5).
  - `--stat`: aggregation mode (`sum`, `mean`, or `normalized`).
  - `--exclude-sample-substr`: exclude samples containing this substring.
  - `--meta-file`: metadata file mapping samples to groups (required for aggregated tables).
  - `--meta-sample-col`: metadata column with sample identifiers.
  - `--meta-group-col`: metadata column to group by.
  - `--sample-col-prefix`: prefix to strip from matrix sample column names.
  - `--sample-id-sep`: separator used to split metadata sample IDs.
  - `--normalize-groups`: optional group normalization.
- **Outputs**
  - TSV tables of retained isoforms per gene / aggregation.

**Example:**

```bash
python -m isoform_distribution.distributions \
  --matrix data/project/expressed_isoforms_matrix.tsv \
  --gtf data/project/expressed_isoforms.gtf \
  --meta-file data/project/meta_data.tsv \
  --output-dir data/project/output/isoform_distributions \
  --meta-sample-col sample_id \
  --meta-group-col cell_type  \
  --sample-col-prefix ENCFF \
  --sample-id-sep _ \
  --normalize-groups heart_brain \
  --stat sum \
  --cutoff-pct 1.5
```

### Dashboard: `isoform_dashboard/dashboard_app.py`
- **Features**
  - Sortable gene rankings table with statistics with filters.
  - Interactive scatter plot: summed vs top isoform entropy, colored by min Spearman.
  - Isoform distribution panel: per‑sample/group bar charts; optional bootstrap CIs (`--ci-file`); switch between **mean** and **sum** tables.
  - Exon structure: color‑coded exons (orange CDS, blue UTR) from GTF; clickable InterPro domains.
  - Support for AlphaFold 3D geometry visualization.
  - Exon-level highlighting in the 3D viewer and protein sequence display.
- **Inputs**
  - `--input-mean`: TSV with mean values.
  - `--input-sum`: TSV with sum values (optional).
  - `--ci-file`: bootstrap CI TSV (optional).
  - `--exons`: GTF file with exon/CDS annotations (optional).
  - `--proteins`: protein FASTA file (optional).
  - `--geometry-dir`: AlphaFold geometry output directory (optional).
  - `--interpro-dir`: InterPro results directory (optional).
- **Outputs**
  - Live app at `http://127.0.0.1:8050` (default).

**Example:**

```bash
python -m isoform_dashboard.dashboard_app \
  --input-mean data/project/output/isoform_distributions/distributions_mean.tsv \
  --input-sum  data/project/output/isoform_distributions/distributions_sum.tsv \
  --ci-file    data/project/output/isoform_distributions/confidence_intervals.tsv \
  --exons      data/project/expressed_isoforms.gtf \
  --proteins   data/project/proteins.fasta \
  --geometry-dir data/project/output/alphafold_geometry \
  --interpro-dir data/project/output/interpro_results
```

Then open `http://127.0.0.1:8050` in your browser.

### InterPro domains: `interpro/interpro_scan.py`

Submits protein sequences to the EBI InterPro Scan REST API.

- Reads peptide sequences from FASTA.
- Submits a transcript’s protein sequence, polls until done, saves JSON/TSV.

**Example:**

```bash
python3 -m interpro.interpro_scan \
  TRANSCRIPT_ID \
  --fasta data/project/expressed_isoforms_PEP.fasta \
  --output data/project/output/interpro_results/TRANSCRIPT_ID.json
```

### AlphaFold geometry conversion: `utilities/extract_3d_geometry.py`

Convert AlphaFold Server job outputs into per‑residue geometry summaries and
HTML viewers (including optional per‑exon highlight viewers).

- **Inputs**
  - `--alphafold-dir`: AlphaFold Server output folder containing job subfolders.
  - `--output-dir`: target folder for geometry outputs.
  - `--model-index`: which ranked model to use (default: 0).
  - `--all-models`: process model indices 0–4 for every job (optional).
  - `--gtf`: GTF file to generate per‑exon highlighted HTML viewers (optional).
- **Outputs**
  - `viewer_modelN.html` interactive HTML viewer per job.
  - `viewer_modelN_exonX.html` (when `--gtf` is provided).

**Example (batch mode, per‑exon viewers):**
```bash
python -m utilities.extract_3d_geometry \
  --alphafold-dir data/project/alphafold_jobs \
  --output-dir data/project/output/alphafold_geometry \
  --gtf data/project/expressed_isoforms.gtf
```

## Bootstrapping overview

**Script:** `isoform_distribution/bootstrap_isoform_means.py`

- Start from the isoform expression matrix.
- For each bootstrap iteration:
  - Resample columns (samples) with replacement.
  - Compute mean expression per isoform.
- For each isoform:
  - Sort its bootstrap means.
  - Drop the lowest and highest 2.5% (for N=1000, drop 25 values at each tail) → 95% CI.
- Generates global CIs plus optional group-specific CIs (if grouping metadata is provided).
- CIs are used as error bars on isoform distribution plots (mean dataset).

**Usage:**

```bash
python3 -m isoform_distribution.bootstrap_isoform_means \
  --input data/project/expressed_isoforms_matrix.tsv \
  --output-dir data/project/output/isoform_distributions \
  --iterations 1000 \
  --seed 42
```