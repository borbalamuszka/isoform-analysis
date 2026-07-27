# Isoform Analysis Tools & Interactive Dashboard

- Analyze gene isoform expression and generate distribution tables/plots.
- Compute per‑gene metrics (entropy, Spearman) and bootstrap isoform means with CIs.
- Query protein domains via EBI InterPro Scan REST API.
- Create 3D structure geometry visualisations from AlphaFold models.
- Explore results in an interactive Dash dashboard with exon, domain, and 3D protein structure visualization.
- Works with generic biological datasets; grouping and labels depend on the metadata provided (e.g., `condition`, `cell_type`, `tissue`).

---

## Data Flow & Architecture

The dashboard features a **modular, multi-stream architecture**. Input streams can be provided independently; missing optional inputs degrade gracefully while preserving core functionality.

```
  ┌─────────────────────────────────────────────────────────┐
  │ 1. EXPRESSION DATA STREAM (Primary Data)                │
  │    • Pre-computed TSVs (distributions_mean.tsv / sum.tsv) │ ──┐
  │                    -- OR --                             │   │
  │    • Raw Expression Matrix + Sample Metadata            │   │
  │      └──► Processed via Precomputation Suite (UI / CLI)  │   │
  └─────────────────────────────────────────────────────────┘   │
                                                                │
  ┌─────────────────────────────────────────────────────────┐   │
  │ 2. EXON STRUCTURE & GTF STREAM (Key Primary Input)      │   │
  │    • Exons GTF File (expressed_isoforms.gtf)            │   ├──► ┌──────────────────────────────────┐
  │      └──► Enables 2D Exon Structure Viewer, CDS/UTR   │   │    │                                  │
  │           Coordinates & Transcript-to-Gene Mappings    │   │    │  Interactive Dash Web Dashboard  │
  └─────────────────────────────────────────────────────────┘   │    │           (Port 8050)            │
                                                                │    │                                  │
  ┌─────────────────────────────────────────────────────────┐   │    │ • Graceful Degradation:          │
  │ 3. CO-EXPRESSION NETWORK STREAM (Optional Network Data) │   │    │   Launches with GTF annotations, │
  │    • Precomputed Sparse Matrices (coexpression/*.npz)    │───┼───►│   expression datasets, or        │
  └─────────────────────────────────────────────────────────┘   │    │   3D models independently!       │
                                                                │    └──────────────────────────────────┘
  ┌─────────────────────────────────────────────────────────┐   │
  │ 4. PROTEIN & 3D ASSETS STREAM (Optional Enhancements)   │   │
  │    • Protein FASTA (amino acid sequence viewer)         │   │
  │    • AlphaFold 3D Geometry Dir (3D molecular viewer)   │───┘
  │    • InterPro Scan Results Dir (functional domains)     │
  └─────────────────────────────────────────────────────────┘
```

### Input Optionality & Graceful Degradation
- **Exon Structure GTF (Key Primary Input)**: Loading an `expressed_isoforms.gtf` file enables 2D exon structures, CDS/UTR color-coding, transcript-to-gene resolution, and structural searches even if expression matrix datasets are not yet loaded.
- **Pre-computed vs. Raw Expression**: Supply pre-calculated `distributions_mean.tsv` / `distributions_sum.tsv` tables **OR** process raw count matrices using the built-in web precomputation tool.
- **Co-expression Networks**: Optional correlation matrices (`coexpression/*.npz` and `*.pkl`) can be supplied or calculated via `utilities/compute_coexpression.py`.
- **Protein & 3D Structural Enhancements**: Protein FASTA, AlphaFold 3D Geometry, and InterPro Scan results are optional enhancements that unlock amino acid sequence views, interactive 3D structures, and functional domain overlays.

---

## Quick start: run the dashboard locally

1. **Clone and enter the repo**
    ```bash
    git clone https://github.com/borbalamuszka/isoform-analysis.git
    cd isoform-analysis
    ```

    > **Alternatives to cloning the repo?**
    > The simplest supported workflow is cloning this GitHub repo. You could instead download a ZIP from GitHub ("Code → Download ZIP") and unpack it, but you still need the same Python environment, data files, and commands as described below.

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

   Choose one of two options to launch and configure the dashboard:

   ### Option A: Start empty & load data using the UI (Recommended)
   Run the dashboard without any command-line arguments:
   ```bash
   python -m isoform_dashboard.dashboard_app
   ```
   Then open `http://127.0.0.1:8050` in your web browser. A welcome modal will prompt you to select your data files using a built-in visual file browser. Once selected, click **Apply Changes & Reload** to load the dataset dynamically.
   
   > 💡 **Help Guide:** Once the dashboard is running, click the **💡 Help Guide** button in the top header to view a comprehensive tabbed manual covering data structures, widget capabilities, and precomputations.

   > 📊 **Web Preprocessing Suite:** You can run raw matrix processing, bootstrap CI generation, co-expression calculations, InterPro API queries, and network drive mapping directly from the UI via the **📊 Preprocess Data** button.

   > 🔄 **Dynamic Dataset Switching:** Re-open the **Configure Data Sources** dialog at any time to dynamically load new data files or switch datasets without restarting the server.

   ### Option B: Run with command-line arguments (Classic)
   Specify input files as command-line arguments:
   ```bash
   python -m isoform_dashboard.dashboard_app \
     --input-mean   data/project/output/isoform_distributions/distributions_mean.tsv \
     --input-sum    data/project/output/isoform_distributions/distributions_sum.tsv \
     --ci-file      data/project/output/isoform_distributions/confidence_intervals.tsv \
     --exons        data/project/expressed_isoforms.gtf \
     --proteins     data/project/proteins.fasta \
     --geometry-dir data/project/output/alphafold_geometry \
     --interpro-dir data/project/output/interpro_results \
     --gene-coexpression-dir data/project/output/coexpression
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
     --interpro-dir data/project/output/interpro_results `
     --gene-coexpression-dir data/project/output/coexpression
   ```

   Then open `http://127.0.0.1:8050` in your browser.

---

## Troubleshooting & FAQs

- **Python Version Requirements:** Ensure Python 3.8+ (Python 3.9 - 3.12 recommended) is installed. You can check your version by running `python --version` or `python3 --version`.
- **Python Command Not Recognized:** On Windows, if `python` is not found, try using `py -m venv isoform_dashboard_env` or ensure "Add Python to PATH" was selected during Python installation. On macOS/Linux, use `python3`.
- **Dependency Installation Issues:** If `pip install -r requirements.txt` fails, upgrade pip first inside your virtual environment using `python -m pip install --upgrade pip`.
- **Port in use:** If port 8050 is busy, pass `--port 8051` when launching.
- **Stopping the app:** Press `Ctrl+C` in your terminal to stop the server.
- **Windows Paths & Multi-line Commands:** Paths are relative to the repo root. Forward slashes (`/`) work on Windows in Python. If copying multi-line commands in PowerShell, ensure backticks (`` ` ``) are used. For Windows cmd.exe, use carets (`^`) for line continuation.
- **Alternatives to Cloning:** You can download a ZIP file from GitHub ("Code → Download ZIP") and extract it if git is unavailable.
- **InterPro REST API Limits:** When querying EBI InterPro Scan services, always provide a contact email to prevent client IP rate-limiting.

---

## Input Data Specification

- **Expression Data File**
  - Tab-separated or CSV matrix of raw count / CPM expression (rows = `transcript_id`, columns = sample IDs).
  - Used by distribution scripts and dashboard.

- **Positions / GTF File**
  - Standard GTF genomic annotation containing exon and CDS features (`gene_id`, `transcript_id`, coordinates).
  - Used by exon structure viewers, coordinate mappers, and domain renderers.

- **Sample Metadata File**
  - Table mapping sample IDs to experimental biological variables (e.g. `condition`, `cell_type`, `tissue`, `disease_status`).
  - Required for sample aggregation and bootstrap CI calculation.

- **Optional Confidence Intervals File**
  - `confidence_intervals.tsv` containing precomputed 95% bootstrap limits (`ci_lower`, `ci_upper`) per isoform and per group.
  - Used by the dashboard to draw error bars on distribution charts.

---

## Main Components & Precomputations

### 1. Isoform Distributions: `isoform_distribution/distributions.py`

Filters low-contribution isoforms by percentage cutoff and aggregates sample columns into Mean/Sum distribution tables.

- **Key Arguments**
  - `--matrix`: Isoform expression matrix (rows = `transcript_id`).
  - `--gtf`: GTF annotation file.
  - `--meta-file`: Metadata table mapping samples to biological groups (**required**).
  - `--output-dir`: Output folder for generated tables (`distributions_mean.tsv`, `distributions_sum.tsv`).
  - `--meta-sample-col`: Column with sample IDs matching matrix headers (default: `sample_id`).
  - `--meta-group-col`: Metadata column to aggregate by (e.g., `condition`, `cell_type`, `tissue`; default: `group`).
  - `--cutoff-pct`: Minimum percentage contribution of an isoform to its gene to be retained (default: 1.5%).

**Biological Metadata Examples (`--meta-file`)**

Example 1 (Neuroscience Project - Group by `condition` or `tissue`):
```tsv
sample_id	condition	tissue	donor_id
Sample_01_Caudate	Control	Caudate	Donor_A
Sample_02_DLPFC	Control	DLPFC	Donor_A
Sample_03_Caudate	MDD	Caudate	Donor_B
Sample_04_DLPFC	MDD	DLPFC	Donor_B
```

Example 2 (Single-Cell / Multi-tissue Project - Group by `cell_type`):
```tsv
sample_id	cell_type	tissue
185VYD_cardiomyocyte	Cardiomyocyte	Heart
206TQZ_astrocyte	Astrocyte	Brain
```

**CLI Execution Example (Group by condition):**
```bash
python -m isoform_distribution.distributions \
  --matrix data/neuro_project/expressed_isoforms_matrix.tsv \
  --gtf data/neuro_project/expressed_isoforms.gtf \
  --meta-file data/neuro_project/metadata.tsv \
  --output-dir data/neuro_project/output/isoform_distributions \
  --table-type aggregated \
  --meta-sample-col sample_id \
  --meta-group-col condition
```

---

### 2. Bootstrapping Confidence Intervals: `isoform_distribution/bootstrap_isoform_means.py`

Resamples sample columns with replacement across iterations (default: 1000) to calculate robust 95% confidence intervals (2.5th and 97.5th percentiles) for isoform expression levels.

- CIs are saved to `confidence_intervals.tsv` and rendered as error bars on distribution plots.

**CLI Execution Example:**
```bash
python3 -m isoform_distribution.bootstrap_isoform_means \
  --input data/neuro_project/expressed_isoforms_matrix.tsv \
  --output-dir data/neuro_project/output/isoform_distributions \
  --iterations 1000 \
  --seed 42
```

---

### 3. Co-expression Networks: `utilities/compute_coexpression.py`

Computes pairwise gene-to-gene and isoform-to-isoform correlation matrices across expression samples.

- **Correlation Method**: Supports `--method spearman` (rank-based, robust against outliers, default) or `--method pearson` (linear).
- **RAM Optimization**:
  - Sets `--threshold 0.3` (or higher) to filter weak edges and maintain a low memory footprint in sparse `.npz` format.
  - Processes correlations in computational blocks via `--chunk_size 1000`.

**CLI Execution Example:**
```bash
python3 -m utilities.compute_coexpression \
  --input data/project/expressed_isoforms_matrix.tsv \
  --output-dir data/project/output/coexpression \
  --threshold 0.3 \
  --chunk_size 1000 \
  --method spearman
```

---

### 4. InterPro Protein Domains: `interpro/interpro_scan.py`

Queries the EBI InterPro Scan REST API for domain predictions based on translated peptide sequences.

- Requires `--email` to identify client requests and prevent API throttling.
- Saves results as per-transcript JSON files (`TRANSCRIPT_ID.json`).

**CLI Execution Example:**
```bash
python3 -m interpro.interpro_scan \
  TRANSCRIPT_ID \
  --fasta data/project/proteins.fasta \
  --email user@institution.edu \
  --output data/project/output/interpro_results/TRANSCRIPT_ID.json
```

---

### 5. 3D Structure Geometry: `utilities/extract_3d_geometry.py`

Converts AlphaFold structure job outputs into per-residue geometry summaries and interactive 3D HTML models with exon highlighting.

**CLI Execution Example:**
```bash
python -m utilities.extract_3d_geometry \
  --alphafold-dir data/project/alphafold_jobs \
  --output-dir data/project/output/alphafold_geometry \
  --gtf data/project/expressed_isoforms.gtf
```