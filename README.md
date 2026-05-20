# Quick start: run the dashboard locally

1. **Clone and enter the repo**
  ```bash
  git clone https://github.com/borbala19/isoform-analysis.git
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

   Minimal (use all built‑in defaults for inputs/paths):

   ```bash
   python -m isoform_dashboard.dashboard_app
   ```

   Explicit (override inputs or use custom files):

   ```bash
   python -m isoform_dashboard.dashboard_app \
     --input-mean   data/neuro_project/output/isoform_distributions/distributions_condition_mean.tsv \
     --input-sum    data/neuro_project/output/isoform_distributions/distributions_condition_sum.tsv \
     --ci-file      data/neuro_project/output/isoform_distributions/confidence_intervals.tsv \
     --exons        data/neuro_project/expressed_isoforms.gtf \
     --proteins     data/neuro_project/proteins.fasta \
     --geometry-dir data/neuro_project/output/alphafold_geometry
   ```

   Windows PowerShell (line continuation with backticks):

   ```powershell
   python -m isoform_dashboard.dashboard_app `
     --input-mean   data/neuro_project/output/isoform_distributions/distributions_condition_mean.tsv `
     --input-sum    data/neuro_project/output/isoform_distributions/distributions_condition_sum.tsv `
     --ci-file      data/neuro_project/output/isoform_distributions/confidence_intervals.tsv `
     --exons        data/neuro_project/expressed_isoforms.gtf `
     --proteins     data/neuro_project/proteins.fasta `
     --geometry-dir data/neuro_project/output/alphafold_geometry
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

# Gene Isoform Expression Analysis

Tools for analysing gene isoform expression across brain regions and conditions, generating isoform distribution tables/plots, computing per‑gene metrics (entropy, Spearman, TVD), bootstrapping isoform means with CIs, querying protein domains via InterPro Scan, and exploring everything in an interactive Dash dashboard.

## Main components

### Isoform distributions & correlations

**Distributions:** `isoform_distribution/distributions.py`

Generates per‑gene isoform distribution tables and multi‑page PDF plots.

- **Features**
  - Filters low‑contribution isoforms by percentage cutoff.
  - Aggregation modes:
    - `sum`: summed expression across samples.
    - `mean`: mean expression across samples.
    - `normalized`: per‑gene per‑sample normalization (isoforms sum to 1); global value = mean normalized contribution.
  - Individual‑sample and aggregated views (by brain region or condition).
  - Optional CI markers on global bars when a CI file is provided.
- **Inputs**
  - `expressed_isoforms_matrix.txt/tsv`: isoform expression matrix (rows = transcript_id, columns = samples).
  - Positions / GTF file for transcript → gene mapping.
  - Optional `confidence_intervals.tsv`.
- **Outputs**
  - TSV tables of retained isoforms per gene / aggregation.
  - PDF plots of isoform distributions.

**Correlations:** `isoform_distribution/correlation.py`

Computes per‑gene metrics (entropy, Spearman, TVD) and interactive scatter plots.

- **Metrics**
  - Top and summed isoform entropy.
  - Pairwise Spearman correlations.
  - Pairwise TVD on per‑gene, per‑sample normalized isoform profiles.
- **Outputs**
  - `<input>_gene_correlations.tsv` with metrics.
  - HTML scatter plots (summed vs top entropy, min Spearman, max TVD, etc.).

### InterPro domains: `interpro/interpro_scan.py`

Submits protein sequences to the EBI InterPro Scan REST API.

- Reads peptide sequences from FASTA.
- Submits a transcript’s protein sequence, polls until done, saves JSON/TSV.
- Default output under `data/neuro_project/output/interpro_results/`.

### Dash dashboard: `isoform_dashboard/`

Interactive web app for exploring isoform expression, entropy, correlations, and exon structure.

- **Entry point:** `isoform_dashboard/dashboard_app.py`
- **Features**
  - Dataset toggle: switch between **mean** and **sum** tables.
  - Interactive scatter plot: summed vs top isoform entropy, colored by min Spearman.
  - Isoform distribution panel: per‑sample/group bar charts; optional bootstrap CIs (`--ci-file`).
  - Exon structure: color‑coded exons (orange CDS, blue UTR) from GTF.
  - Sortable gene rankings table with statistics and export functionality.
  - Support for AlphaFold geometry visualization.
- **Inputs**
  - `--input-mean`: TSV with mean values.
  - `--input-sum`: TSV with sum values (optional).
  - `--ci-file`: bootstrap CI TSV (optional).
  - `--exons`: GTF file with exon/CDS annotations (optional).
- **Outputs**
  - Live app at `http://127.0.0.1:8050` (default).
  - `highlighted_genes_YYYYMMDD_HHMMSS.txt` on export.


## Input data

- **`expressed_isoforms_matrix.tsv` / `.txt`**
  - Rows: `transcript_id`.
  - Columns: sample names (e.g. `Region_Condition_Age_Subject`).
  - Used by distribution scripts, correlation, and dashboard.
  - Need not be pre‑normalized; TVD is computed on normalized profiles internally.

- **Positions / GTF file**
  - Genomic positions and metadata (chrom, exon_start/end, cds_start/end, gene_id, transcript_id, …).
  - Used to map transcripts to genes and to build exon diagrams.

- **Optional `confidence_intervals.tsv`**
  - Bootstrap mean and CI bounds per isoform, plus region/condition‑specific CIs.
  - Used by `distributions.py` / `plots.py` and the dashboard to overlay error bars.

## Typical usage

### Isoform distributions

```bash
python3 -m isoform_distribution.distributions \
  --table-type aggregated \
  --stat sum \
  --cutoff-pct 1.5
```

### InterPro Scan

```bash
python3 -m interpro.interpro_scan \
  G10110.21.nnc \
  --fasta data/neuro_project/expressed_isoforms_PEP.fasta \
  --output data/neuro_project/output/interpro_results/G10110.21.nnc.json
```

Windows PowerShell:

```powershell
python -m interpro.interpro_scan `
  G10110.21.nnc `
  --fasta data/neuro_project/expressed_isoforms_PEP.fasta `
  --output data/neuro_project/output/interpro_results/G10110.21.nnc.json
```

For Windows cmd.exe, replace the backticks with carets (`^`) or paste the command on one line.

### Dashboard

Minimal (use all built‑in defaults for inputs/paths):

```bash
python -m isoform_dashboard.dashboard_app
```

Explicit (override inputs or use custom files):

```bash
python -m isoform_dashboard.dashboard_app \
  --input-mean data/neuro_project/output/isoform_distributions/distributions_condition_mean.tsv \
  --input-sum  data/neuro_project/output/isoform_distributions/distributions_condition_sum.tsv \
  --ci-file    data/neuro_project/output/isoform_distributions/confidence_intervals.tsv \
  --exons      data/neuro_project/expressed_isoforms.gtf
```

Windows PowerShell:

```powershell
python -m isoform_dashboard.dashboard_app `
  --input-mean data/neuro_project/output/isoform_distributions/distributions_condition_mean.tsv `
  --input-sum  data/neuro_project/output/isoform_distributions/distributions_condition_sum.tsv `
  --ci-file    data/neuro_project/output/isoform_distributions/confidence_intervals.tsv `
  --exons      data/neuro_project/expressed_isoforms.gtf
```

For Windows cmd.exe, replace the backticks with carets (`^`) or paste the command on one line.

Then open `http://127.0.0.1:8050` in your browser.

## Bootstrapping overview

**Script:** `isoform_distribution/bootstrap_isoform_means.py`

- Start from the isoform expression matrix (adult samples only; fetal excluded via `--include-key adult`).
- For each bootstrap iteration:
  - Resample columns (samples) with replacement.
  - Compute mean expression per isoform.
- For each isoform:
  - Sort its bootstrap means.
  - Drop the lowest and highest 2.5% (for N=1000, drop 25 values at each tail) → 95% CI.
- Generates global CIs plus region-specific and condition-specific CIs.
- CIs are used as error bars on isoform distribution plots (mean dataset).

**Usage:**

```bash
python3 -m isoform_distribution.bootstrap_isoform_means \
  --input data/neuro_project/expressed_isoforms_matrix.tsv \
  --output-dir data/neuro_project/output/isoform_distributions \
  --iterations 1000 \
  --seed 42 \
  --include-key adult
```

Windows PowerShell:

```powershell
python -m isoform_distribution.bootstrap_isoform_means `
  --input data/neuro_project/expressed_isoforms_matrix.tsv `
  --output-dir data/neuro_project/output/isoform_distributions `
  --iterations 1000 `
  --seed 42 `
  --include-key adult
```

For Windows cmd.exe, replace the backticks with carets (`^`) or paste the command on one line.