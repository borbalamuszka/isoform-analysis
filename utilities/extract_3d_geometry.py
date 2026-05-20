"""
extract_3d_geometry.py
======================
Extracts 3D geometry data from AlphaFold Server job outputs (.cif models +
full_data JSON) and produces per job:

  1. A per-residue geometry CSV with:
       - residue index, chain, residue name
       - Cα coordinates (x, y, z)
       - per-residue mean pLDDT  (from full_data JSON atom_plddts)
       - per-residue PAE self-score (diagonal of PAE matrix)

  2. An interactive HTML viewer coloured by pLDDT (AlphaFold DB
     style: dark-blue ≥ 90, light-blue 70–90, yellow 50–70, orange < 50)

  3. Per-exon highlighted HTML viewers (one per CDS exon) when --gtf is
     supplied.  Each file shows the full protein greyed out with only that
     exon's residues coloured by pLDDT – ready to be served by the dashboard
     when the user clicks an exon in the left-hand structure panel.

Usage
-----
Run on the entire folds directory (processes every job sub-folder):

    conda run -n minimal_env python extract_3d_geometry.py \\
        --alphafold-dir ./Alphafold/folds_2026_02_25_15_56 \\
        --output-dir   ./data/neuro_project/output/alphafold_geometry \\
        --model-index  0          # which ranked model to use (0 = best)

Include per-exon viewers by pointing at the GTF:

    conda run -n minimal_env python extract_3d_geometry.py \\
        --alphafold-dir ./Alphafold/folds_2026_02_25_15_56 \\
        --output-dir   ./data/neuro_project/output/alphafold_geometry \\
        --gtf          ./data/neuro_project/expressed_isoforms.gtf

Or point at a single job folder:

    conda run -n minimal_env python extract_3d_geometry.py \\
        --job-dir ./Alphafold/folds_2026_02_25_15_56/g186356nnc_ensg000001028377 \\
        --output-dir ./output_geometry \\
        --gtf        ./expressed_isoforms.gtf
"""

import argparse
import json
import logging
import pathlib
import re
import sys
import warnings

import numpy as np
import pandas as pd
try:
    import py3Dmol
    HAS_PY3DMOL = True
except ImportError:
    HAS_PY3DMOL = False
    warnings.warn("py3Dmol not found – HTML viewers will be skipped. "
                  "Install with: pip install py3Dmol")

# ── optional: BioPython for CIF parsing ────────────────────────────────────
try:
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB import PDBIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    warnings.warn("BioPython not found – falling back to manual CIF parser. "
                  "Install with: pip install biopython")

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s | %(message)s")
log = logging.getLogger(__name__)

# ───────────────────────────── pLDDT colour scheme ─────────────────────────
# Matches the AlphaFold DB colour scheme exactly.
PLDDT_COLOURS = [
    (90, "royalblue"),    # very high
    (70, "cornflowerblue"),  # confident
    (50, "yellow"),        # low
    (0,  "orange"),        # very low
]

def plddt_to_colour(value: float) -> str:
    """Return AlphaFold-style hex/name colour for a pLDDT value."""
    for threshold, colour in PLDDT_COLOURS:
        if value >= threshold:
            return colour
    return "orange"


# ───────────────────────────── GTF helpers ─────────────────────────────────

def _parse_gtf_attributes(attributes: str):
    """Return (transcript_id, gene_id) from a GTF attributes string."""
    t_match = re.search(r'transcript_id "([^"]+)"', attributes)
    g_match = re.search(r'gene_id "([^"]+)"', attributes)
    return (
        t_match.group(1) if t_match else None,
        g_match.group(1) if g_match else None,
    )


def build_transcript_id_lookup(isoforms_by_gene: dict) -> dict:
    """Build a reverse map: encoded_id → original_gtf_transcript_id.

    AlphaFold Server job folders are named by taking the GTF transcript ID,
    removing all dots, and lowercasing (e.g. ``G70076.141.nnc`` →
    ``g70076141nnc``).  This map lets us recover the original ID from the
    folder-encoded version.

    Args:
        isoforms_by_gene: Output of :func:`parse_gtf_isoforms`.

    Returns:
        Dict mapping encoded (lowercase, no-dots) transcript IDs to their
        original GTF transcript IDs.
    """
    lookup: dict = {}
    for isoforms in isoforms_by_gene.values():
        for t_id in isoforms:
            encoded = t_id.replace(".", "").lower()
            lookup[encoded] = t_id
    return lookup


def parse_gtf_isoforms(gtf_path: pathlib.Path) -> dict:
    """Parse a GTF file and return exon structures with CDS overlap info.

    Returns:
        ``isoforms_by_gene``: gene_id → transcript_id → list of exon dicts.
        Each exon dict has keys: exon_start, exon_end, cds_start, cds_end, strand.
    """
    from collections import defaultdict

    exons_by_transcript: dict = defaultdict(list)
    cds_by_transcript: dict   = defaultdict(list)
    gene_for_transcript: dict = {}

    with open(gtf_path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            start        = int(parts[3])
            end          = int(parts[4])
            strand       = parts[6] if len(parts) > 6 else "+"
            t_id, g_id   = _parse_gtf_attributes(parts[8])
            if not t_id or not g_id:
                continue
            gene_for_transcript[t_id] = g_id
            if feature_type == "exon":
                exons_by_transcript[t_id].append({"start": start, "end": end, "strand": strand})
            elif feature_type == "CDS":
                cds_by_transcript[t_id].append({"start": start, "end": end})

    isoforms_by_gene: dict = defaultdict(lambda: defaultdict(list))

    for t_id, exons in exons_by_transcript.items():
        g_id = gene_for_transcript.get(t_id)
        if not g_id:
            continue
        cds_regions = cds_by_transcript.get(t_id, [])
        for exon in exons:
            es, ee = exon["start"], exon["end"]
            cds_start = cds_end = None
            for cds in cds_regions:
                if cds["start"] <= ee and cds["end"] >= es:
                    ov_s = max(cds["start"], es)
                    ov_e = min(cds["end"], ee)
                    if cds_start is None:
                        cds_start, cds_end = ov_s, ov_e
                    else:
                        cds_start = min(cds_start, ov_s)
                        cds_end   = max(cds_end,   ov_e)
            isoforms_by_gene[g_id][t_id].append({
                "exon_start": es,
                "exon_end":   ee,
                "cds_start":  cds_start,
                "cds_end":    cds_end,
                "strand":     exon.get("strand", "+"),
            })

    return isoforms_by_gene


def compute_exon_residue_range(exon_idx: int, exons_sorted: list):
    """Return the (1-based) protein residue range covered by a CDS exon.

    Args:
        exon_idx: 0-based position of the exon in *exons_sorted*.
        exons_sorted: Exon dicts sorted ascending by exon_start.

    Returns:
        ``(res_start, res_end)`` (both inclusive, 1-based), or ``None`` if the
        exon has no CDS overlap.
    """
    target = exons_sorted[exon_idx]
    if (target.get("cds_start") is None
            or target.get("cds_end") is None
            or target["cds_end"] <= target["cds_start"]):
        return None

    strand = target.get("strand", "+")

    # On the minus strand the transcript is read in reverse genomic order.
    if strand == "-":
        ordered          = list(reversed(exons_sorted))
        target_order_idx = len(exons_sorted) - 1 - exon_idx
    else:
        ordered          = exons_sorted
        target_order_idx = exon_idx

    cumulative_bp = 0
    for i, exon in enumerate(ordered):
        if i == target_order_idx:
            break
        cs = exon.get("cds_start")
        ce = exon.get("cds_end")
        if cs is not None and ce is not None and ce > cs:
            cumulative_bp += ce - cs

    target_cds_bp = target["cds_end"] - target["cds_start"]
    res_start     = cumulative_bp // 3 + 1
    res_end       = (cumulative_bp + target_cds_bp - 1) // 3 + 1
    return (res_start, res_end)


# ───────────────────────────── CIF helpers ─────────────────────────────────

def _parse_cif_atom_sites_manual(cif_path: pathlib.Path) -> pd.DataFrame:
    """
    Minimal mmCIF _atom_site parser that works without BioPython.
    Returns a DataFrame with columns matching _atom_site field names.
    """
    lines = cif_path.read_text().splitlines()

    # Collect _atom_site column names in order
    in_atom_loop = False
    columns: list[str] = []
    data_rows: list[list[str]] = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line == "loop_":
            # Peek ahead to see if this is an _atom_site loop
            j = i + 1
            peek_cols: list[str] = []
            while j < len(lines) and lines[j].strip().startswith("_atom_site."):
                peek_cols.append(lines[j].strip())
                j += 1
            if peek_cols:
                in_atom_loop = True
                columns = [c.split(".")[1] for c in peek_cols]
                i = j  # skip to data rows
                continue
            else:
                in_atom_loop = False

        if in_atom_loop:
            if line.startswith("_") or line == "#" or line == "":
                in_atom_loop = False
                i += 1
                continue
            # Split whitespace-delimited record
            data_rows.append(line.split())

        i += 1

    if not columns or not data_rows:
        raise ValueError(f"No _atom_site data found in {cif_path}")

    df = pd.DataFrame(data_rows, columns=columns[:len(data_rows[0])])
    return df


def parse_cif_calpha(cif_path: pathlib.Path) -> pd.DataFrame:
    """
    Parse a CIF model file and return a DataFrame of Cα atoms with columns:
      chain, res_seq, res_name, x, y, z, b_factor (= pLDDT in AlphaFold CIFs)
    """
    if HAS_BIOPYTHON:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", str(cif_path))
        records = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if "CA" in residue:
                        atom = residue["CA"]
                        x, y, z = atom.get_vector().get_array()
                        records.append({
                            "chain":    chain.id,
                            "res_seq":  residue.get_id()[1],
                            "res_name": residue.get_resname(),
                            "x": float(x),
                            "y": float(y),
                            "z": float(z),
                            "b_factor": float(atom.get_bfactor()),
                        })
            break  # only first MODEL block
        return pd.DataFrame(records)

    else:
        # Manual parser fallback
        df = _parse_cif_atom_sites_manual(cif_path)
        # Keep only Cα atoms
        ca = df[df["label_atom_id"] == "CA"].copy()
        ca = ca.rename(columns={
            "label_asym_id": "chain",
            "label_seq_id":  "res_seq",
            "label_comp_id": "res_name",
            "Cartn_x": "x",
            "Cartn_y": "y",
            "Cartn_z": "z",
            "B_iso_or_equiv": "b_factor",
        })
        for col in ["x", "y", "z", "b_factor", "res_seq"]:
            ca[col] = pd.to_numeric(ca[col], errors="coerce")
        return ca[["chain", "res_seq", "res_name", "x", "y", "z", "b_factor"]].reset_index(drop=True)


def cif_to_pdb_string(cif_path: pathlib.Path) -> str:
    """
    Convert a CIF file to a PDB-format string for py3Dmol.
    Uses BioPython if available; otherwise returns the raw CIF text
    (py3Dmol can also accept 'mmcif' format directly).
    """
    if HAS_BIOPYTHON:
        from io import StringIO
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", str(cif_path))
        io = PDBIO()
        io.set_structure(structure)
        buf = StringIO()
        io.save(buf)
        return buf.getvalue()
    else:
        # Return raw CIF; py3Dmol accepts mmcif format
        return cif_path.read_text()


# ───────────────────────────── JSON helpers ────────────────────────────────

def load_full_data(json_path: pathlib.Path) -> dict:
    return json.loads(json_path.read_text())


def compute_residue_plddt(full_data: dict) -> pd.DataFrame:
    """
    Average atom-level pLDDT values per residue token.
    Returns DataFrame with columns: chain, res_seq, mean_plddt
    """
    atom_plddts    = np.array(full_data["atom_plddts"],    dtype=float)
    atom_chains    = full_data["atom_chain_ids"]
    token_chains   = full_data["token_chain_ids"]
    token_res_ids  = full_data["token_res_ids"]

    # full_data atom_plddts are per-atom; average by residue token index
    # token_res_ids gives the residue sequence number for each token
    n_tokens = len(token_res_ids)
    n_atoms  = len(atom_plddts)

    # Distribute atoms evenly across tokens (AlphaFold Server outputs
    # atoms in residue order; use integer division to assign token index)
    atoms_per_token = n_atoms // n_tokens
    remainder       = n_atoms % n_tokens

    records = []
    atom_idx = 0
    for tok in range(n_tokens):
        # Some tokens may have one extra atom (remainder)
        count = atoms_per_token + (1 if tok < remainder else 0)
        chunk = atom_plddts[atom_idx: atom_idx + count]
        records.append({
            "chain":      token_chains[tok],
            "res_seq":    token_res_ids[tok],
            "mean_plddt": float(np.mean(chunk)) if len(chunk) > 0 else float("nan"),
        })
        atom_idx += count

    return pd.DataFrame(records)


def compute_pae_diagonal(full_data: dict) -> np.ndarray:
    """Return the diagonal of the PAE matrix (self-PAE per residue)."""
    pae = np.array(full_data["pae"], dtype=float)
    return np.diag(pae)


# ───────────────────────────── HTML viewer ─────────────────────────────────

def make_html_viewer(
    cif_path: pathlib.Path,
    residue_df: pd.DataFrame,
    output_html: pathlib.Path,
    width: int = 700,
    height: int = 600,
) -> None:
    """
    Generate a standalone HTML file with an interactive py3Dmol viewer
    coloured by pLDDT (AlphaFold DB colour scheme).
    """
    if not HAS_PY3DMOL:
        log.warning("py3Dmol not available – skipping HTML viewer for %s", cif_path.name)
        return

    fmt = "mmcif" if not HAS_BIOPYTHON else "pdb"
    mol_str = cif_to_pdb_string(cif_path) if HAS_BIOPYTHON else cif_path.read_text()

    view = py3Dmol.view(width=width, height=height)
    view.addModel(mol_str, fmt)

    # Cartoon backbone, white base
    view.setStyle({"cartoon": {"color": "white"}})

    # Colour each residue by its mean pLDDT
    for _, row in residue_df.iterrows():
        colour = plddt_to_colour(row["mean_plddt"])
        view.setStyle(
            {"chain": row["chain"], "resi": int(row["res_seq"])},
            {"cartoon": {"color": colour}},
        )

    view.zoomTo()
    view.spin(False)

    # Embed into a standalone HTML with a legend
    html_body = view._make_html()

    legend_html = """
<div style="font-family:Arial,sans-serif; padding:8px; background:#1a1a2e;
            color:white; display:inline-block; border-radius:6px; margin-top:6px;">
  <b>pLDDT colour key</b><br>
  <span style="color:royalblue">&#9632;</span> &ge;90 Very high &nbsp;
  <span style="color:cornflowerblue">&#9632;</span> 70&ndash;89 Confident &nbsp;
  <span style="color:yellow">&#9632;</span> 50&ndash;69 Low &nbsp;
  <span style="color:orange">&#9632;</span> &lt;50 Very low
</div>
"""
    full_html = f"""<!DOCTYPE html>
<html>
<head><meta charset="utf-8"><title>AlphaFold 3D – {cif_path.stem}</title></head>
<body style="background:#1a1a2e; margin:0; padding:10px;">
<h2 style="color:white;font-family:Arial;">{cif_path.stem}</h2>
{html_body}
{legend_html}
</body>
</html>"""

    output_html.write_text(full_html)
    log.info("  → HTML viewer saved: %s", output_html)


# ───────────────────────────── Exon HTML viewers ───────────────────────────

def make_exon_html_viewers(
    transcript_id: str,
    isoforms_by_gene: dict,
    plddt_df: pd.DataFrame,
    pdb_str: str,
    job_out: pathlib.Path,
    model_index: int = 0,
    width: int = 700,
    height: int = 600,
) -> list:
    """Generate per-exon highlighted HTML viewers next to the full viewer.

    For each CDS exon a new HTML is written that shows the full protein with:
      - all residues greyed out (lightgrey cartoon, low opacity)
      - that exon's residues coloured by pLDDT at full opacity

    Non-CDS exons get ``None`` in the returned list.

    Files are written as ``viewer_model{N}_exon{i}.html`` (0-based exon index).
    Existing files are reused without regeneration.

    Args:
        transcript_id: E.g. ``"G16681.47.nnc"``.
        isoforms_by_gene: Output of :func:`parse_gtf_isoforms`.
        plddt_df: DataFrame with columns ``res_seq`` and ``mean_plddt``.
        pdb_str: PDB-format string of the model (already converted from CIF).
        job_out: Directory to write HTML files into.
        model_index: Model index used for file naming.
        width: Viewer width in pixels.
        height: Viewer height in pixels.

    Returns:
        List of file paths (``str``) or ``None``, one entry per exon (sorted
        ascending by genomic position). Returns ``[]`` when py3Dmol is
        unavailable or the transcript is not found in the GTF.
    """
    if not HAS_PY3DMOL:
        log.warning("make_exon_html_viewers: py3Dmol not installed – skipping")
        return []

    # Locate this transcript in the GTF-derived isoforms dict
    transcript_exons = None
    for isoforms in isoforms_by_gene.values():
        if transcript_id in isoforms:
            transcript_exons = isoforms[transcript_id]
            break

    if transcript_exons is None:
        log.warning("make_exon_html_viewers: transcript %s not found in GTF", transcript_id)
        return []

    exons_sorted = sorted(transcript_exons, key=lambda e: e["exon_start"])

    plddt_by_res = dict(zip(
        plddt_df["res_seq"].astype(int),
        plddt_df["mean_plddt"].astype(float),
    ))

    results = []
    for exon_idx, exon in enumerate(exons_sorted):
        res_range = compute_exon_residue_range(exon_idx, exons_sorted)

        if res_range is None:
            results.append(None)
            continue

        res_start, res_end = res_range
        out_path = job_out / f"viewer_model{model_index}_exon{exon_idx}.html"

        if out_path.is_file():
            results.append(str(out_path))
            log.debug("make_exon_html_viewers: reusing existing %s", out_path)
            continue

        try:
            view = py3Dmol.view(width=width, height=height)
            view.addModel(pdb_str, "pdb")
            view.setStyle({"cartoon": {"color": "lightgrey", "opacity": 0.25}})

            for res_seq in range(res_start, res_end + 1):
                plddt  = plddt_by_res.get(res_seq, 50.0)
                colour = plddt_to_colour(plddt)
                view.setStyle(
                    {"chain": "A", "resi": res_seq},
                    {"cartoon": {"color": colour, "opacity": 1.0}},
                )

            view.zoomTo()
            view.spin(False)
            html_body = view._make_html()

            cds_str = (
                f"CDS: {exon['cds_start']:,}–{exon['cds_end']:,}"
                if exon.get("cds_start") else ""
            )
            full_html = (
                "<!DOCTYPE html>\n<html>\n<head>"
                "<meta charset='utf-8'>"
                f"<title>Exon {exon_idx + 1} – {transcript_id}</title>"
                "</head>\n"
                "<body style='background:#ffffff; margin:0; padding:10px;'>\n"
                f"<h3 style='font-family:Arial;margin:0 0 4px;'>"
                f"Exon {exon_idx + 1} / {len(exons_sorted)}"
                f" &nbsp; residues {res_start}–{res_end}"
                + (f" &nbsp; ({cds_str})" if cds_str else "")
                + "</h3>\n"
                + html_body
                + "\n</body>\n</html>"
            )

            out_path.write_text(full_html)
            results.append(str(out_path))
            log.info("  → Exon viewer saved: %s (res %d–%d)", out_path, res_start, res_end)

        except Exception as exc:
            log.warning("make_exon_html_viewers: failed for exon %d of %s: %s",
                        exon_idx, transcript_id, exc)
            results.append(None)

    return results


# ───────────────────────────── Job processing ──────────────────────────────

def find_job_files(job_dir: pathlib.Path, model_index: int):
    """
    Locate the CIF model and full_data JSON for a given model_index inside a
    job directory.  Returns (cif_path, full_data_path) or raises FileNotFoundError.
    """
    cif_glob       = list(job_dir.glob(f"*_model_{model_index}.cif"))
    full_data_glob = list(job_dir.glob(f"*_full_data_{model_index}.json"))

    def _pick(files, label):
        if not files:
            raise FileNotFoundError(
                f"No {label} (model index {model_index}) in {job_dir}"
            )
        return files[0]

    return (
        _pick(cif_glob,       "CIF model"),
        _pick(full_data_glob, "full_data JSON"),
    )


def process_job(
    job_dir: pathlib.Path,
    output_dir: pathlib.Path,
    model_index: int = 0,
    isoforms_by_gene: dict = None,
    transcript_id_lookup: dict = None,
) -> bool:
    """
    Process one AlphaFold Server job folder.
    Returns True on success, False on failure.

    Args:
        job_dir: AlphaFold Server job directory.
        output_dir: Root output directory.
        model_index: Which ranked model to use.
        isoforms_by_gene: Optional exon structure dict from :func:`parse_gtf_isoforms`.
                          When provided, per-exon highlighted HTML viewers are
                          generated alongside the full-protein viewer.
        transcript_id_lookup: Optional dict mapping encoded transcript IDs
                              (lowercase, dots removed) to original GTF
                              transcript IDs.  Built via
                              :func:`build_transcript_id_lookup`.
    """
    job_id = job_dir.name
    log.info("Processing job: %s", job_id)

    # Derive output folder name: transcript part only (strip trailing _ensg... suffix)
    encoded_id = re.sub(r'_ensg\d+$', '', job_id, flags=re.IGNORECASE)

    # Recover the original GTF transcript ID (e.g. G70076.141.nnc) from the
    # folder-encoded version (e.g. g70076141nnc) using the lookup if available.
    if transcript_id_lookup and encoded_id in transcript_id_lookup:
        transcript_id = transcript_id_lookup[encoded_id]
    else:
        transcript_id = encoded_id

    try:
        cif_path, full_data_path = find_job_files(job_dir, model_index)
    except FileNotFoundError as e:
        log.warning("  Skipping: %s", e)
        return False

    # ── Parse CIF for Cα coordinates ──
    try:
        ca_df = parse_cif_calpha(cif_path)
    except Exception as e:
        log.warning("  CIF parse failed for %s: %s", job_id, e)
        return False

    # ── Load full_data JSON ──
    full_data = load_full_data(full_data_path)

    # ── Per-residue pLDDT from JSON ──
    try:
        plddt_df = compute_residue_plddt(full_data)
    except Exception as e:
        log.warning("  pLDDT aggregation failed for %s: %s", job_id, e)
        # Fall back to b_factor from CIF
        plddt_df = ca_df[["chain", "res_seq"]].copy()
        plddt_df["mean_plddt"] = ca_df["b_factor"]

    # ── PAE diagonal ──
    try:
        pae_diag = compute_pae_diagonal(full_data)
    except Exception as e:
        log.warning("  PAE extraction failed for %s: %s", job_id, e)
        pae_diag = np.full(len(plddt_df), float("nan"))

    plddt_df["pae_self"] = pae_diag[:len(plddt_df)]

    # ── Merge coordinates + pLDDT + PAE ──
    geometry_df = pd.merge(
        ca_df,
        plddt_df[["chain", "res_seq", "mean_plddt", "pae_self"]],
        on=["chain", "res_seq"],
        how="left",
    )
    geometry_df.insert(0, "transcript_id", transcript_id)
    geometry_df.insert(1, "model_index", model_index)

    # ── Save per-residue CSV ──
    # Use encoded_id for the output folder name to stay consistent with
    # existing files on disk (AlphaFold Server lowercases the job name).
    job_out = output_dir / encoded_id
    job_out.mkdir(parents=True, exist_ok=True)
    geom_csv = job_out / f"geometry_model{model_index}.csv"
    geometry_df.to_csv(geom_csv, index=False, float_format="%.4f")
    log.info("  → Geometry CSV saved: %s", geom_csv)

    # ── Full-protein HTML viewer ──
    html_path = job_out / f"viewer_model{model_index}.html"
    make_html_viewer(cif_path, plddt_df, html_path)

    # ── Per-exon highlighted HTML viewers (requires --gtf) ──
    if isoforms_by_gene is not None and HAS_PY3DMOL:
        # Re-use the PDB string that make_html_viewer already produced; we need
        # it as a plain string here.  Convert from CIF once and reuse.
        pdb_str = cif_to_pdb_string(cif_path) if HAS_BIOPYTHON else cif_path.read_text()
        make_exon_html_viewers(
            transcript_id=transcript_id,
            isoforms_by_gene=isoforms_by_gene,
            plddt_df=plddt_df[["res_seq", "mean_plddt"]],
            pdb_str=pdb_str,
            job_out=job_out,
            model_index=model_index,
            width=700,
            height=600,
        )

    return True


# ───────────────────────────── CLI ────────────────────────────────────────

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Extract 3D geometry data from AlphaFold Server job outputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--alphafold-dir", metavar="DIR",
        help="Directory containing multiple job sub-folders (batch mode).",
    )
    source.add_argument(
        "--job-dir", metavar="DIR",
        help="Single job folder to process.",
    )
    p.add_argument(
        "--output-dir", metavar="DIR", required=True,
        help="Directory where CSVs and HTML viewers are written.",
    )
    p.add_argument(
        "--model-index", type=int, default=0, metavar="N",
        help="Which ranked model to extract (0 = best, default: 0).",
    )
    p.add_argument(
        "--all-models", action="store_true",
        help="Process all 5 model indices (0–4) per job.",
    )
    p.add_argument(
        "--gtf", metavar="FILE",
        help=(
            "GTF annotation file (e.g. expressed_isoforms.gtf). "
            "When provided, per-exon highlighted HTML viewers are generated "
            "alongside the full-protein viewer for every processed job."
        ),
    )
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine which model indices to process
    model_indices = list(range(5)) if args.all_models else [args.model_index]

    # Parse GTF if provided
    isoforms_by_gene = None
    transcript_id_lookup = None
    if args.gtf:
        gtf_path = pathlib.Path(args.gtf)
        if not gtf_path.is_file():
            log.error("GTF file not found: %s", gtf_path)
            sys.exit(1)
        log.info("Parsing GTF: %s", gtf_path)
        isoforms_by_gene = parse_gtf_isoforms(gtf_path)
        n_transcripts = sum(len(v) for v in isoforms_by_gene.values())
        log.info("GTF parsed: %d genes, %d transcripts", len(isoforms_by_gene), n_transcripts)
        transcript_id_lookup = build_transcript_id_lookup(isoforms_by_gene)
        log.info("Transcript ID lookup built: %d entries", len(transcript_id_lookup))

    # Collect job directories
    if args.alphafold_dir:
        alphafold_dir = pathlib.Path(args.alphafold_dir)
        job_dirs = sorted(
            d for d in alphafold_dir.iterdir()
            if d.is_dir() and not d.name.startswith(".")
        )
        log.info("Found %d job folders in %s", len(job_dirs), alphafold_dir)
    else:
        job_dirs = [pathlib.Path(args.job_dir)]

    # Process
    n_ok = 0
    for job_dir in job_dirs:
        for midx in model_indices:
            if process_job(
                job_dir,
                output_dir,
                model_index=midx,
                isoforms_by_gene=isoforms_by_gene,
                transcript_id_lookup=transcript_id_lookup,
            ):
                n_ok += 1

    if n_ok == 0:
        log.error("No jobs were processed successfully.")
        sys.exit(1)

    log.info("Done – processed %d job/model pair(s).", n_ok)


if __name__ == "__main__":
    main()
