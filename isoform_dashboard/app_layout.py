"""Dashboard layout and callbacks.

This module creates the Dash app and defines all interactive callbacks.
"""
import logging
import os
import functools
import pandas as pd
from datetime import datetime
from dash import Dash, dcc, html, Input, Output, State, ALL, dash_table, callback_context, no_update
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import string

log = logging.getLogger(__name__)

from .data_processing import (
    prepare_table_data,
    compute_min_spearman_per_gene,
    compute_gene_ranking,
    compute_gene_ranking_by_expression,
    gene_has_cds,
    calculate_entropy_and_correlation,
)
from .alphafold_geometry import (
    build_alphafold_geometry_mapping,
    resolve_alphafold_geometry,
    extend_geometry_mapping_by_sequence,
    discover_exon_viewers,
)
from .visualizations import (
    fig_summed_vs_top_entropy_colored_by_min_spearman,
    create_exon_visualization,
    fig_isoform_sample_panels
)
from .coexpression_network import create_coexpression_widget, generate_network_elements
from .config import Colors, Dimensions

# Residue-range helper (pure function, no heavy deps)
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utilities.extract_3d_geometry import compute_exon_residue_range


def get_system_roots():
    if os.name == 'nt':
        drives = []
        for letter in string.ascii_uppercase:
            drive = f"{letter}:\\"
            if os.path.exists(drive):
                drives.append(drive)
        return drives
    else:
        roots = ['/']
        for sub in ['/mnt', '/media', '/home']:
            if os.path.exists(sub):
                roots.append(sub)
        return roots


def get_directory_contents(path, target_type='all'):
    try:
        path = os.path.abspath(path)
        if not os.path.exists(path):
            return [], f"Path does not exist: {path}"
        
        items = []
        parent = os.path.dirname(path)
        if parent and parent != path:
            items.append({'name': '..', 'is_dir': True, 'path': parent})
            
        for name in os.listdir(path):
            full_path = os.path.join(path, name)
            try:
                is_dir = os.path.isdir(full_path)
            except Exception:
                continue
            
            if target_type == 'dir' and not is_dir:
                continue
                
            size = None
            if not is_dir:
                try:
                    size = os.path.getsize(full_path)
                except Exception:
                    pass
            
            items.append({
                'name': name,
                'is_dir': is_dir,
                'path': full_path,
                'size': size
            })
            
        items.sort(key=lambda x: (not x['is_dir'], x['name'].lower()))
        return items, None
    except Exception as e:
        return [], str(e)


def _axis_ranges_for_results(results_df: pd.DataFrame):
    if results_df is None or results_df.empty:
        return None, None
    x_min = results_df["summed_isoform_entropy"].min()
    x_max = results_df["summed_isoform_entropy"].max()
    y_min = results_df["top_isoform_entropy"].min()
    y_max = results_df["top_isoform_entropy"].max()
    return [x_min, x_max], [y_min, y_max]


def create_app(df_mean: pd.DataFrame, df_sum: pd.DataFrame, results_df_mean: pd.DataFrame,
               results_df_sum: pd.DataFrame, sample_cols, ci_df: pd.DataFrame, ci_columns: list,
               global_col_mean: str, global_col_sum: str, isoforms_by_gene, gene_names: dict,
               has_sum: bool,
               geometry_dir: str = None,
               protein_sequences: dict = None,
               domain_mapping: dict = None,
               default_ranking: str = "spearman",
               gene_coexpression=None,
               gene_coexpression_idx=None,
               isoform_coexpression=None,
               isoform_coexpression_idx=None,
               path_mean=None,
               path_sum=None,
               path_gtf=None,
               path_proteins=None,
               path_interpro_dir=None,
               path_coexpression_dir=None,
               path_ci=None):
    """Create and configure the Dash application.

    Args:
        df_mean: Mean expression data
        df_sum: Sum expression data
        results_df_mean: Computed results for mean data
        results_df_sum: Computed results for sum data
        sample_cols: List of sample column names
        ci_df: Confidence interval data
        ci_columns: List of CI column names
        global_col_mean: Global column name for mean data
        global_col_sum: Global column name for sum data
        isoforms_by_gene: Dictionary of exon structures
        gene_names: Dictionary mapping gene_id to gene_name
        has_sum: Whether sum data is available
        geometry_dir: Path to alphafold_geometry output directory (optional).
        protein_sequences: Dict mapping transcript_id -> amino-acid sequence (optional).
        domain_mapping: Dict mapping transcript_id -> list of domain dicts (optional).
        default_ranking: Which ranking to sort the table by at startup.
            'spearman' (default) or 'expression'.
        gene_coexpression: Precomputed sparse gene coexpression matrix (optional).
        gene_coexpression_idx: Index array for gene coexpression matrix (optional).
        isoform_coexpression: Precomputed sparse isoform coexpression matrix (optional).
        isoform_coexpression_idx: Index array for isoform coexpression matrix (optional).

    Returns:
        Configured Dash app
    """
    protein_sequences = protein_sequences or {}
    domain_mapping = domain_mapping or {}
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.title = "Isoform Analysis Dashboard"

    # Build AlphaFold geometry mapping (transcript_nodots_gene_nodots -> file paths)
    af_geometry_mapping = build_alphafold_geometry_mapping(geometry_dir) if geometry_dir else {}
    if geometry_dir and not af_geometry_mapping:
        log.warning("create_app: --geometry-dir=%r was given but no geometry entries were found", geometry_dir)
    elif not geometry_dir:
        log.warning("create_app: --geometry-dir not provided; 3D structure features will be disabled")
    else:
        log.info("create_app: AlphaFold geometry mapping loaded: %d entries from %r",
                 len(af_geometry_mapping), geometry_dir)

    # Discover pre-generated per-exon HTML viewers
    if af_geometry_mapping:
        discover_exon_viewers(af_geometry_mapping)

    # Extend geometry mapping with aliases for transcripts that share an
    # identical amino acid sequence with a transcript that already has geometry.
    if af_geometry_mapping and protein_sequences:
        n_aliases = extend_geometry_mapping_by_sequence(af_geometry_mapping, protein_sequences)
        if n_aliases:
            log.info(
                "create_app: added %d sequence-identity geometry aliases; "
                "total mapping size now %d",
                n_aliases, len(af_geometry_mapping),
            )

    global_style = {
        "fontFamily": "Arial, sans-serif",
        "fontSize": "14px",
        "color": "#333333",
    }

    # Mutable state container
    state = {
        'df_mean': df_mean,
        'df_sum': df_sum,
        'results_df_mean': results_df_mean,
        'results_df_sum': results_df_sum,
        'sample_cols': sample_cols or [],
        'global_col_mean': global_col_mean or "",
        'global_col_sum': global_col_sum or "",
        'has_sum': has_sum,
        'isoforms_by_gene': isoforms_by_gene or {},
        'gene_names': gene_names or {},
        'af_geometry_mapping': af_geometry_mapping or {},
        'domain_mapping': domain_mapping or {},
        'protein_sequences': protein_sequences or {},
        'gene_coexpression': gene_coexpression,
        'gene_coexpression_idx': gene_coexpression_idx,
        'isoform_coexpression': isoform_coexpression,
        'isoform_coexpression_idx': isoform_coexpression_idx,
        
        # Derived tables
        'table_df_mean': None,
        'table_df_sum': None,
        'table_mean_records': [],
        'table_sum_records': [],
        'transcript_to_exons': {},
        'genes_with_cds': set(),
        'genes_with_3d': set(),
        'genes_with_domains': set(),
        
        # Paths
        'path_mean': path_mean or "",
        'path_sum': path_sum or "",
        'path_gtf': path_gtf or "",
        'path_geom': geometry_dir or "",
        'path_coexp': path_coexpression_dir or "",
        'path_fasta': path_proteins or "",
        'path_interpro': path_interpro_dir or "",
        'path_ci': path_ci or "",
        'ci_df': ci_df,
        'ci_columns': ci_columns or [],
        'ci_dict': {},
        
        # Scatter ranges
        'scatter_axis_ranges': {'mean': (None, None), 'sum': (None, None)},
        
        # Loading progress state
        'loading_progress': {'step': 0, 'total': 7, 'msg': '', 'done': False, 'error': None, 'updated_time': 0}
    }

    def recompute_derived_state():
        state['genes_with_cds'] = set()
        if state['isoforms_by_gene']:
            state['genes_with_cds'] = {
                gene_id for gene_id in state['isoforms_by_gene'].keys()
                if gene_has_cds(gene_id, state['isoforms_by_gene'])
            }

        state['genes_with_3d'] = set()
        if state['isoforms_by_gene'] and state['af_geometry_mapping']:
            for gene_id, transcripts in state['isoforms_by_gene'].items():
                for transcript_id in transcripts.keys():
                    key = transcript_id.replace(".", "").lower()
                    if key in state['af_geometry_mapping']:
                        state['genes_with_3d'].add(gene_id)
                        break

        state['genes_with_domains'] = set()
        if state['isoforms_by_gene'] and state['domain_mapping']:
            for gene_id, transcripts in state['isoforms_by_gene'].items():
                for transcript_id in transcripts.keys():
                    if state['domain_mapping'].get(transcript_id):
                        state['genes_with_domains'].add(gene_id)
                        break

        state['transcript_to_exons'] = {}
        if state['isoforms_by_gene']:
            for _gene_exons in state['isoforms_by_gene'].values():
                for _tid, _exons in _gene_exons.items():
                    state['transcript_to_exons'][_tid] = sorted(_exons, key=lambda e: e["exon_start"])

        state['ci_dict'] = state['ci_df'].to_dict('index') if state['ci_df'] is not None else {}

        if state['results_df_mean'] is not None and not state['results_df_mean'].empty:
            state['results_df_mean']["min_spearman"] = compute_min_spearman_per_gene(state['results_df_mean'])
            state['results_df_mean']["rank"] = compute_gene_ranking(state['results_df_mean'])
            state['results_df_mean']["rank_by_expression"] = compute_gene_ranking_by_expression(state['results_df_mean'])

        if state['has_sum'] and state['results_df_sum'] is not None and not state['results_df_sum'].empty:
            state['results_df_sum']["min_spearman"] = compute_min_spearman_per_gene(state['results_df_sum'])
            state['results_df_sum']["rank"] = compute_gene_ranking(state['results_df_sum'])
            state['results_df_sum']["rank_by_expression"] = compute_gene_ranking_by_expression(state['results_df_sum'])

        if state['results_df_mean'] is not None and not state['results_df_mean'].empty:
            state['table_df_mean'] = prepare_table_data(
                state['results_df_mean'], state['isoforms_by_gene'], state['gene_names'],
                af_geometry_mapping=state['af_geometry_mapping'],
                default_ranking=default_ranking,
                protein_sequences=state['protein_sequences'],
                domain_mapping=state['domain_mapping']
            )
            state['table_mean_records'] = state['table_df_mean'].to_dict('records')
        else:
            state['table_df_mean'] = None
            state['table_mean_records'] = []

        if state['has_sum'] and state['results_df_sum'] is not None and not state['results_df_sum'].empty:
            state['table_df_sum'] = prepare_table_data(
                state['results_df_sum'], state['isoforms_by_gene'], state['gene_names'],
                af_geometry_mapping=state['af_geometry_mapping'],
                default_ranking=default_ranking,
                protein_sequences=state['protein_sequences'],
                domain_mapping=state['domain_mapping']
            )
            state['table_sum_records'] = state['table_df_sum'].to_dict('records')
        else:
            state['table_df_sum'] = None
            state['table_sum_records'] = []

        state['scatter_axis_ranges'] = {
            'mean': _axis_ranges_for_results(state['results_df_mean']),
            'sum': _axis_ranges_for_results(state['results_df_sum']) if state['has_sum'] else (None, None)
        }

    # Populate initial derived state
    recompute_derived_state()
    state['recompute_derived_state'] = recompute_derived_state

    welcome_banner = html.Div(
        id="welcome-banner",
        style={
            "display": "block" if df_mean is None or df_mean.empty else "none",
            "backgroundColor": "#ebf5fb",
            "borderLeft": "6px solid #3498db",
            "padding": "16px 20px",
            "borderRadius": "4px",
            "marginBottom": "20px",
            "fontFamily": "Segoe UI, Tahoma, Geneva, Verdana, sans-serif"
        },
        children=[
            html.H4("Welcome to the Isoform Analysis Dashboard!", style={"margin": "0 0 8px 0", "color": "#2980b9"}),
            html.P([
                "No dataset is currently loaded. To visualize your data, please click the ",
                html.Strong("⚙ Configure Data Sources"),
                " button in the top-right corner to load your Mean Expression TSV file and other optional data files (such as Sum Expression TSV, Exons GTF, AlphaFold 3D structures, etc.)."
            ], style={"margin": "0", "color": "#2c3e50", "fontSize": "14px"})
        ]
    )

    settings_popup_layout = html.Div(
        id="settings-popup-window",
        style={
            "display": "none",
            "position": "fixed",
            "zIndex": "2000",
            "left": "0",
            "top": "0",
            "width": "100%",
            "height": "100%",
            "overflow": "auto",
            "backgroundColor": "rgba(0,0,0,0.5)",
            "backdropFilter": "blur(4px)",
        },
        children=[
            html.Div(
                style={
                    "backgroundColor": "#ffffff",
                    "margin": "5% auto",
                    "padding": "24px",
                    "border": "1px solid #888",
                    "width": "65%",
                    "maxWidth": "750px",
                    "borderRadius": "8px",
                    "boxShadow": "0 4px 20px rgba(0,0,0,0.2)",
                    "fontFamily": "Segoe UI, Tahoma, Geneva, Verdana, sans-serif",
                },
                children=[
                    html.Div([
                        html.H3("Configure Data Sources", style={"margin": "0 0 20px 0", "color": "#2c3e50", "display": "inline-block"}),
                        html.Button("×", id="btn-close-settings", n_clicks=0, style={
                            "float": "right",
                            "fontSize": "28px",
                            "fontWeight": "bold",
                            "border": "none",
                            "background": "none",
                            "cursor": "pointer",
                            "color": "#aaa",
                        })
                    ]),
                    
                    # Field 1: Mean Expression (TSV)
                    html.Div([
                        html.Label("Mean Expression TSV File (Required for visualization)", className="settings-label"),
                        html.Div([
                            dcc.Input(id="input-path-mean", value=state['path_mean'], placeholder="e.g. Z:\\distributions_condition_mean.tsv", className="settings-input"),
                            html.Button("Browse...", id="btn-browse-mean", n_clicks=0, className="btn-browse")
                        ], className="settings-input-container")
                    ], className="settings-row"),

                    # Field 2: Sum Expression (TSV)
                    html.Div([
                        html.Label("Sum Expression TSV File (Optional)", className="settings-label"),
                        html.Div([
                            dcc.Input(id="input-path-sum", value=state['path_sum'], placeholder="e.g. Z:\\distributions_condition_sum.tsv", className="settings-input"),
                            html.Button("Browse...", id="btn-browse-sum", n_clicks=0, className="btn-browse")
                        ], className="settings-input-container")
                    ], className="settings-row"),

                    # Field 3: Exons GTF
                    html.Div([
                        html.Label("Exons GTF File (Optional, required for structure visualization)", className="settings-label"),
                        html.Div([
                            dcc.Input(id="input-path-gtf", value=state['path_gtf'], placeholder="e.g. Z:\\expressed_isoforms.gtf", className="settings-input"),
                            html.Button("Browse...", id="btn-browse-gtf", n_clicks=0, className="btn-browse")
                        ], className="settings-input-container")
                    ], className="settings-row"),

                    # Field 4: AlphaFold Geometry Directory
                    html.Div([
                        html.Label("AlphaFold Geometry Directory (Optional, required for 3D structures)", className="settings-label"),
                        html.Div([
                            dcc.Input(id="input-path-geom", value=state['path_geom'], placeholder="e.g. Z:\\alphafold_geometry", className="settings-input"),
                            html.Button("Browse...", id="btn-browse-geom", n_clicks=0, className="btn-browse")
                        ], className="settings-input-container")
                    ], className="settings-row"),

                    # Field 5: Co-expression Directory
                    html.Div([
                        html.Label("Co-expression Matrices Directory (Optional, containing gene_coexpression.npz, etc.)", className="settings-label"),
                        html.Div([
                            dcc.Input(id="input-path-coexp", value=state['path_coexp'], placeholder="e.g. Z:\\gene_coexpression_dir", className="settings-input"),
                            html.Button("Browse...", id="btn-browse-coexp", n_clicks=0, className="btn-browse")
                        ], className="settings-input-container")
                    ], className="settings-row"),

                    # Field 6: Protein FASTA
                    html.Div([
                        html.Label("Protein Sequences FASTA (Optional)", className="settings-label"),
                        html.Div([
                            dcc.Input(id="input-path-fasta", value=state['path_fasta'], placeholder="e.g. Z:\\proteins.fasta", className="settings-input"),
                            html.Button("Browse...", id="btn-browse-fasta", n_clicks=0, className="btn-browse")
                        ], className="settings-input-container")
                    ], className="settings-row"),

                    # Field 7: InterPro Results Directory (Optional)
                    html.Div([
                        html.Label("InterPro Results Directory (Optional)", className="settings-label"),
                        html.Div([
                            dcc.Input(id="input-path-interpro", value=state['path_interpro'], placeholder="e.g. Z:\\interpro_dir", className="settings-input"),
                            html.Button("Browse...", id="btn-browse-interpro", n_clicks=0, className="btn-browse")
                        ], className="settings-input-container")
                    ], className="settings-row"),

                    # Field 8: Confidence Intervals File
                    html.Div([
                        html.Label("Bootstrap Confidence Intervals TSV File (Optional)", className="settings-label"),
                        html.Div([
                            dcc.Input(id="input-path-ci", value=state['path_ci'], placeholder="e.g. Z:\\confidence_intervals.tsv", className="settings-input"),
                            html.Button("Browse...", id="btn-browse-ci", n_clicks=0, className="btn-browse")
                        ], className="settings-input-container")
                    ], className="settings-row"),

                    # Progress Container
                    html.Div(
                        id="progress-container",
                        style={"display": "none", "marginTop": "20px", "marginBottom": "20px"},
                        children=[
                            html.Div(id="progress-status-msg", style={"fontWeight": "bold", "marginBottom": "8px", "color": "#2c3e50"}),
                            html.Div(
                                style={"width": "100%", "backgroundColor": "#e0e0e0", "borderRadius": "4px", "height": "16px", "overflow": "hidden"},
                                children=[
                                    html.Div(id="progress-bar-fill", style={"width": "0%", "height": "100%", "backgroundColor": "#2ecc71", "transition": "width 0.2s"})
                                ]
                            )
                        ]
                    ),

                    # Feedback message
                    html.Div(id="settings-feedback", style={
                        "marginTop": "20px",
                        "padding": "12px",
                        "borderRadius": "4px",
                        "display": "none",
                        "fontSize": "13px"
                    }),

                    # Footer
                    html.Div([
                        html.Button("Apply Changes & Reload", id="btn-apply-settings", n_clicks=0, style={
                            "backgroundColor": "#2ecc71",
                            "color": "white",
                            "border": "none",
                            "padding": "10px 20px",
                            "borderRadius": "4px",
                            "cursor": "pointer",
                            "fontWeight": "bold",
                            "fontSize": "14px",
                            "marginRight": "10px",
                            "transition": "background-color 0.15s"
                        }),
                        html.Button("Cancel", id="btn-cancel-settings", n_clicks=0, style={
                            "backgroundColor": "#95a5a6",
                            "color": "white",
                            "border": "none",
                            "padding": "10px 20px",
                            "borderRadius": "4px",
                            "cursor": "pointer",
                            "fontWeight": "bold",
                            "fontSize": "14px",
                            "transition": "background-color 0.15s"
                        })
                    ], style={"marginTop": "25px", "textAlign": "right"})
                ]
            )
        ]
    )

    preprocess_popup_layout = html.Div(
        id="preprocess-popup-window",
        style={
            "display": "none",
            "position": "fixed",
            "zIndex": "2000",
            "left": "0",
            "top": "0",
            "width": "100%",
            "height": "100%",
            "overflow": "auto",
            "backgroundColor": "rgba(0,0,0,0.5)",
            "backdropFilter": "blur(4px)",
        },
        children=[
            html.Div(
                style={
                    "backgroundColor": "#ffffff",
                    "margin": "5% auto",
                    "padding": "24px",
                    "border": "1px solid #888",
                    "width": "65%",
                    "maxWidth": "750px",
                    "borderRadius": "8px",
                    "boxShadow": "0 4px 20px rgba(0,0,0,0.2)",
                    "fontFamily": "Segoe UI, Tahoma, Geneva, Verdana, sans-serif",
                },
                children=[
                    html.Div([
                        html.H3("Preprocess Data & Calculate Statistics", style={"margin": "0 0 20px 0", "color": "#2c3e50", "display": "inline-block"}),
                        html.Button("×", id="btn-close-preprocess", n_clicks=0, style={
                            "float": "right",
                            "fontSize": "28px",
                            "fontWeight": "bold",
                            "border": "none",
                            "background": "none",
                            "cursor": "pointer",
                            "color": "#aaa",
                        })
                    ]),

                    html.Div([
                        html.Label("Select Precomputation Tool to Run", className="settings-label"),
                        dcc.Dropdown(
                            id="dropdown-select-tool",
                            options=[
                                {'label': 'Calculate Distribution & Confidence Intervals (CIs)', 'value': 'dist'},
                                {'label': 'Calculate Co-expression Matrix', 'value': 'coexp'},
                                {'label': 'Query InterPro Protein Domains (EBI Scan API)', 'value': 'interpro'},
                                {'label': '🔗 Map Remote Windows Network Share (net use)', 'value': 'map-drive'}
                            ],
                            value='dist',
                            clearable=False,
                            className="settings-dropdown",
                            style={"color": "#2c3e50"}
                        )
                    ], className="settings-row"),

                    # Tool 1: Calculate Distribution
                    html.Div(
                        id="tool-config-dist",
                        style={"display": "block"},
                        children=[
                            html.Div([
                                html.Label("Raw Expression Matrix (TSV/CSV)", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-dist-matrix", placeholder="e.g. Z:\\expressed_isoforms_matrix.txt", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-dist-matrix", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Exons GTF File", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-dist-gtf", placeholder="e.g. Z:\\expressed_isoforms.gtf", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-dist-gtf", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Metadata File (CSV/TSV)", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-dist-meta", placeholder="e.g. Z:\\sample_metadata.txt", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-dist-meta", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Output Directory", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-dist-outdir", placeholder="e.g. Z:\\data\\isoform_distributions", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-dist-outdir", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Metadata Sample ID Column", className="settings-label"),
                                dcc.Input(id="tool-dist-meta-sample", value="sample_id", className="settings-input")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Metadata Grouping Column", className="settings-label"),
                                dcc.Input(id="tool-dist-meta-group", value="region", className="settings-input")
                            ], className="settings-row"),
                            html.Div([
                                html.Div([
                                    html.Label("Filtering Percentage Cutoff", className="settings-label"),
                                    dcc.Input(id="tool-dist-cutoff", type="number", value=1.5, step=0.1, className="settings-input")
                                ], style={"width": "48%", "display": "inline-block"}),
                                html.Div([
                                    dcc.Checklist(
                                        id="tool-dist-run-bootstrap",
                                        options=[{'label': ' Run Bootstrap CIs', 'value': 'run_ci'}],
                                        value=['run_ci'],
                                        style={"marginTop": "28px"}
                                    )
                                ], style={"width": "48%", "display": "inline-block", "float": "right"})
                            ], className="settings-row"),
                            
                            # Bootstrap parameters (expandable/indented)
                            html.Div(
                                id="tool-dist-bootstrap-params",
                                style={"display": "block", "paddingLeft": "20px", "borderLeft": "3px solid #eee", "marginTop": "10px"},
                                children=[
                                    html.Div([
                                        html.Div([
                                            html.Label("Bootstrap Iterations", className="settings-label"),
                                            dcc.Input(id="tool-dist-bootstrap-iter", type="number", value=1000, step=100, className="settings-input")
                                        ], style={"width": "31%", "display": "inline-block"}),
                                        html.Div([
                                            html.Label("Bootstrap Include Key", className="settings-label"),
                                            dcc.Input(id="tool-dist-bootstrap-key", value="adult", className="settings-input")
                                        ], style={"width": "31%", "display": "inline-block", "marginLeft": "3%"}),
                                        html.Div([
                                            html.Label("Random Seed", className="settings-label"),
                                            dcc.Input(id="tool-dist-bootstrap-seed", type="number", value=42, className="settings-input")
                                        ], style={"width": "31%", "display": "inline-block", "float": "right"})
                                    ], className="settings-row")
                                ]
                            )
                        ]
                    ),

                    # Tool 2: Calculate Co-expression
                    html.Div(
                        id="tool-config-coexp",
                        style={"display": "none"},
                        children=[
                            html.Div([
                                html.Label("Raw Expression Matrix (TSV/CSV)", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-coexp-expression", placeholder="e.g. Z:\\expressed_isoforms_matrix.tsv", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-coexp-expr", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Exons GTF File", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-coexp-gtf", placeholder="e.g. Z:\\expressed_isoforms.gtf", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-coexp-gtf", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Output Directory", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-coexp-outdir", placeholder="e.g. Z:\\data\\co-expression", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-coexp-outdir", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Correlation Method", className="settings-label"),
                                dcc.Dropdown(
                                    id="tool-coexp-method",
                                    options=[
                                        {"label": "Spearman (Rank-based)", "value": "spearman"},
                                        {"label": "Pearson (Value-based)", "value": "pearson"}
                                    ],
                                    value="spearman",
                                    clearable=False,
                                    className="settings-dropdown",
                                    style={"color": "#2c3e50"}
                                )
                            ], className="settings-row"),
                            html.Div([
                                html.Div([
                                    html.Label("Correlation Threshold", className="settings-label"),
                                    dcc.Input(id="tool-coexp-threshold", type="number", value=0.3, step=0.01, className="settings-input"),
                                    html.Small("Note: Values < 0.2 can cause extremely high memory usage.", style={"color": "#7f8c8d", "fontSize": "11px", "display": "block", "marginTop": "4px"})
                                ], style={"width": "48%", "display": "inline-block", "verticalAlign": "top"}),
                                html.Div([
                                    html.Label("Chunk Size", className="settings-label"),
                                    dcc.Input(id="tool-coexp-chunk", type="number", value=250, step=50, className="settings-input")
                                ], style={"width": "48%", "display": "inline-block", "float": "right", "verticalAlign": "top"})
                            ], className="settings-row", style={"marginBottom": "20px"})
                        ]
                    ),

                    # Tool 3: InterPro
                    html.Div(
                        id="tool-config-interpro",
                        style={"display": "none"},
                        children=[
                            html.Div([
                                html.Label("FASTA Sequences File", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-interpro-fasta", placeholder="e.g. Z:\\proteins.fasta", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-interpro-fasta", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Exons GTF File (Optional, for filtering)", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-interpro-gtf", placeholder="e.g. Z:\\expressed_isoforms.gtf", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-interpro-gtf", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Output Directory", className="settings-label"),
                                html.Div([
                                    dcc.Input(id="tool-interpro-outdir", placeholder="e.g. Z:\\data\\interpro_results", className="settings-input"),
                                    html.Button("Browse...", id="btn-browse-tool-interpro-outdir", n_clicks=0, className="btn-browse")
                                ], className="settings-input-container")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Contact Email (required for EBI services API limit tracking)", className="settings-label"),
                                dcc.Input(id="tool-interpro-email", value="bm708@cam.ac.uk", className="settings-input")
                            ], className="settings-row"),
                            html.Div([
                                dcc.Checklist(
                                    id="tool-interpro-skip",
                                    options=[{'label': ' Skip sequences that already have result files', 'value': 'skip'}],
                                    value=['skip'],
                                    style={"marginTop": "10px"}
                                )
                            ], className="settings-row")
                        ]
                    ),

                    # Tool 4: Map Network Drive
                    html.Div(
                        id="tool-config-map-drive",
                        style={"display": "none"},
                        children=[
                            html.Div([
                                html.Label("Remote Share Path (e.g. \\\\192.168.1.15\\IsoformData)", className="settings-label"),
                                dcc.Input(id="tool-map-share", placeholder="e.g. \\\\192.168.1.15\\IsoformData", className="settings-input")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Drive Letter to Map (e.g. Z:)", className="settings-label"),
                                dcc.Input(id="tool-map-letter", value="Z:", className="settings-input")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Username (Optional)", className="settings-label"),
                                dcc.Input(id="tool-map-user", placeholder="Windows username on the client machine", className="settings-input")
                            ], className="settings-row"),
                            html.Div([
                                html.Label("Password (Optional)", className="settings-label"),
                                dcc.Input(id="tool-map-pass", type="password", placeholder="Windows password on the client machine", className="settings-input")
                            ], className="settings-row")
                        ]
                    ),

                    # Execution Button
                    html.Div([
                        html.Button("▶ Run Precomputation Tool", id="btn-run-tool", n_clicks=0, style={
                            "backgroundColor": "#2ecc71",
                            "color": "white",
                            "border": "none",
                            "padding": "10px 20px",
                            "borderRadius": "4px",
                            "cursor": "pointer",
                            "fontWeight": "bold",
                            "marginRight": "10px",
                            "fontSize": "14px"
                        }),
                        html.Button("⏹ Terminate Running Tool", id="btn-kill-tool", n_clicks=0, style={
                            "backgroundColor": "#e74c3c",
                            "color": "white",
                            "border": "none",
                            "padding": "10px 20px",
                            "borderRadius": "4px",
                            "cursor": "pointer",
                            "fontWeight": "bold",
                            "fontSize": "14px",
                            "display": "none"
                        }),
                        html.Button("Cancel", id="btn-cancel-preprocess", n_clicks=0, style={
                            "backgroundColor": "#95a5a6",
                            "color": "white",
                            "border": "none",
                            "padding": "10px 20px",
                            "borderRadius": "4px",
                            "cursor": "pointer",
                            "fontWeight": "bold",
                            "fontSize": "14px",
                            "float": "right",
                            "transition": "background-color 0.15s"
                        })
                    ], style={"marginTop": "20px"}),

                    # Preprocessing Progress Bar
                    html.Div(
                        id="tool-progress-container",
                        style={"marginTop": "20px", "display": "none"},
                        children=[
                            html.Div(
                                id="tool-progress-status-msg",
                                style={"fontSize": "13px", "fontWeight": "bold", "marginBottom": "5px", "color": "#2c3e50"}
                            ),
                            html.Div(
                                style={
                                    "width": "100%",
                                    "height": "18px",
                                    "backgroundColor": "#f0f0f0",
                                    "borderRadius": "9px",
                                    "overflow": "hidden"
                                },
                                children=[
                                    html.Div(
                                        id="tool-progress-bar-fill",
                                        style={
                                            "width": "0%",
                                            "height": "100%",
                                            "backgroundColor": "#2ecc71",
                                            "transition": "width 0.2s"
                                        }
                                    )
                                ]
                            )
                        ]
                    ),

                    # Console log output
                    html.Div(
                        id="tool-console-container",
                        style={"marginTop": "20px", "display": "none"},
                        children=[
                            html.Div("Live Execution Console Log", style={"fontWeight": "bold", "marginBottom": "5px", "color": "#2c3e50"}),
                            html.Pre(
                                id="tool-console-output",
                                style={
                                    "backgroundColor": "#1e1e1e",
                                    "color": "#d4d4d4",
                                    "padding": "12px",
                                    "borderRadius": "4px",
                                    "height": "200px",
                                    "overflowY": "scroll",
                                    "fontFamily": "Consolas, Courier New, monospace",
                                    "fontSize": "12px",
                                    "whiteSpace": "pre-wrap"
                                }
                            )
                        ]
                    ),
                ]
            )
        ]
    )

    file_selector_layout = html.Div(
        id="file-selector-popup",
        style={
            "display": "none",
            "position": "fixed",
            "zIndex": "2001",
            "left": "0",
            "top": "0",
            "width": "100%",
            "height": "100%",
            "backgroundColor": "rgba(0,0,0,0.6)",
            "backdropFilter": "blur(2px)",
        },
        children=[
            html.Div(
                style={
                    "backgroundColor": "#ffffff",
                    "margin": "8% auto",
                    "padding": "20px",
                    "border": "1px solid #888",
                    "width": "55%",
                    "maxWidth": "650px",
                    "borderRadius": "8px",
                    "boxShadow": "0 5px 25px rgba(0,0,0,0.3)",
                    "fontFamily": "Segoe UI, Tahoma, Geneva, Verdana, sans-serif",
                },
                children=[
                    html.Div([
                        html.H4("Select File or Folder", id="file-selector-title", style={"margin": "0 0 15px 0", "color": "#2c3e50", "display": "inline-block"}),
                        html.Button("×", id="btn-close-file-selector", n_clicks=0, style={
                            "float": "right",
                            "fontSize": "24px",
                            "fontWeight": "bold",
                            "border": "none",
                            "background": "none",
                            "cursor": "pointer",
                            "color": "#aaa",
                        })
                    ]),
                    
                    # Quick Jump System Roots
                    html.Div([
                        html.Label("Quick Jump: ", style={"fontWeight": "bold", "marginRight": "8px", "fontSize": "12px"}),
                        html.Div(id="system-roots-container", style={"display": "inline-block"})
                    ], style={"marginBottom": "10px"}),

                    # Path breadcrumb & manual input row
                    html.Div([
                        html.Label("Current Directory:", style={"fontWeight": "bold", "display": "block", "marginBottom": "3px", "fontSize": "12px"}),
                        html.Div(id="file-browser-breadcrumbs", style={"padding": "6px", "backgroundColor": "#f8f9fa", "border": "1px solid #ddd", "borderRadius": "4px", "marginBottom": "8px", "fontSize": "13px", "wordBreak": "break-all"})
                    ]),
                    
                    # File/folder list box
                    html.Div(id="file-browser-list-container", className="file-browser-list"),
                    
                    # Selected path preview
                    html.Div([
                        html.Label("Selected Path:", style={"fontWeight": "bold", "display": "block", "marginTop": "12px", "marginBottom": "4px", "fontSize": "12px"}),
                        dcc.Input(id="file-selector-selected-input", style={"width": "100%", "padding": "6px", "border": "1px solid #ccc", "borderRadius": "4px", "fontSize": "13px"})
                    ]),
                    
                    # Footer actions
                    html.Div([
                        html.Button("Select Current Folder", id="btn-select-current-folder", n_clicks=0, style={
                            "backgroundColor": "#e67e22", "color": "white", "border": "none", "padding": "8px 14px", "borderRadius": "4px", "cursor": "pointer", "fontWeight": "bold", "marginRight": "8px"
                        }),
                        html.Button("Confirm Selection", id="btn-confirm-file-selection", n_clicks=0, style={
                            "backgroundColor": "#2ecc71", "color": "white", "border": "none", "padding": "8px 14px", "borderRadius": "4px", "cursor": "pointer", "fontWeight": "bold", "marginRight": "8px"
                        }),
                        html.Button("Cancel", id="btn-cancel-file-selection", n_clicks=0, style={
                            "backgroundColor": "#95a5a6", "color": "white", "border": "none", "padding": "8px 14px", "borderRadius": "4px", "cursor": "pointer", "fontWeight": "bold"
                        })
                    ], style={"marginTop": "20px", "textAlign": "right"})
                ]
            )
        ]
    )

    app.layout = html.Div([

        # Settings and File selector stores
        dcc.Store(id="settings-popup-open", data=False),
        dcc.Store(id="file-selector-popup-open", data=False),
        dcc.Store(id="file-browser-current-dir", data=os.getcwd()),
        dcc.Store(id="file-browser-target-field", data=""),
        dcc.Store(id="data-sources-updated", data=0),
        dcc.Interval(id="progress-interval", interval=500, disabled=True, n_intervals=0),
        dcc.Interval(id="tool-interval", interval=1000, disabled=True, n_intervals=0),

        # Dialog Box Popups
        settings_popup_layout,
        preprocess_popup_layout,
        file_selector_layout,

        dcc.Store(id="selected-domain"),
        # Title and header row with settings button
        html.Div([
            html.Div([
                html.H3("Isoform Entropy Dashboard", style={**global_style, "margin": "0 0 10px 0", "display": "inline-block"}),
                html.Button(
                    "⚙ Configure Data Sources",
                    id="btn-open-settings",
                    n_clicks=0,
                    style={
                        "float": "right",
                        "backgroundColor": "#3498db",
                        "color": "white",
                        "border": "none",
                        "padding": "8px 16px",
                        "borderRadius": "4px",
                        "cursor": "pointer",
                        "fontWeight": "bold",
                        "fontSize": "13px",
                        "transition": "background-color 0.2s",
                        "boxShadow": "0 2px 4px rgba(0,0,0,0.1)",
                        "marginLeft": "10px"
                    }
                ),
                html.Button(
                    "📊 Preprocess Data",
                    id="btn-open-preprocess",
                    n_clicks=0,
                    style={
                        "float": "right",
                        "backgroundColor": "#2ecc71",
                        "color": "white",
                        "border": "none",
                        "padding": "8px 16px",
                        "borderRadius": "4px",
                        "cursor": "pointer",
                        "fontWeight": "bold",
                        "fontSize": "13px",
                        "transition": "background-color 0.2s",
                        "boxShadow": "0 2px 4px rgba(0,0,0,0.1)"
                    }
                )
            ], style={"width": "100%"}),
            html.Div([
                html.Label("Dataset: ", style={**global_style, "fontWeight": "bold", "marginRight": "10px"}),
                dcc.RadioItems(
                    id='dataset-toggle',
                    options=[
                        {'label': ' Mean', 'value': 'mean'},
                        {'label': ' Sum', 'value': 'sum', 'disabled': not has_sum}
                    ],
                    value='mean',
                    inline=True,
                    style={'display': 'inline-block'}
                )
            ], style={"marginBottom": "20px", "paddingBottom": "10px", "borderBottom": "2px solid #ddd"}),
        ]),

        welcome_banner,
        
        # Split view container
        html.Div(
            id="main-dashboard-content",
            style={"display": "flex" if df_mean is not None and not df_mean.empty else "none", "width": "100%", "marginBottom": "20px"},
            children=[
            # Left panel - Visualizations
            html.Div([
                # 1. Isoform Distribution (first)
                html.Div([
                    html.Div([
                        html.Div(id="selected-gene-label", style={**global_style, "fontWeight": "bold", "marginBottom": "8px"}),
                        # Add CI toggle if CIs are available
                        html.Div([
                            dcc.Checklist(
                                id='ci-toggle',
                                options=[{'label': ' Show Confidence Intervals', 'value': 'show_ci'}],
                                value=[],
                                style={'display': 'inline-block' if ci_df is not None else 'none'}
                            )
                        ], style={"marginBottom": "8px"}),
                    ]),
                    dcc.Graph(id="isoform-distribution"),
                ], style={"width": "100%", "marginBottom": "10px"}),

                html.Hr(style={"margin": "8px 0"}),
                
                # 2. Exon Structure (second)
                html.Div([
                    html.Div(
                        id="exon-structure-label",
                        style={**global_style, "fontWeight": "bold", "marginBottom": "2px"}
                    ),
                    dcc.Graph(id="exon-visualization"),
                    html.Div(
                        id="selection-details-wrapper",
                        children=[
                            html.Div(
                                id="transcript-details-wrapper",
                                children=[
                                    html.Div(
                                        [
                                            html.Span("Transcript", style={"fontWeight": "bold"}),
                                            dcc.Clipboard(
                                                id="transcript-details-clipboard",
                                                target_id="transcript-details-text",
                                                style={"float": "right"},
                                                title="Copy transcript details",
                                            ),
                                        ],
                                        style={"marginBottom": "6px"},
                                    ),
                                    html.Pre(
                                        id="transcript-details-text",
                                        style={
                                            "whiteSpace": "pre-wrap",
                                            "margin": 0,
                                            "fontFamily": "monospace",
                                            "fontSize": "12px",
                                        },
                                    ),
                                ],
                                style={
                                    "flex": "1",
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "borderRadius": "4px",
                                    "backgroundColor": "#fafafa",
                                    "display": "none",
                                },
                            ),
                            html.Div(
                                id="domain-details-wrapper",
                                children=[
                                    html.Div(
                                        [
                                            html.Span("Domain", style={"fontWeight": "bold"}),
                                            dcc.Clipboard(
                                                id="domain-details-clipboard",
                                                target_id="domain-details-text",
                                                style={"float": "right"},
                                                title="Copy domain details",
                                            ),
                                        ],
                                        style={"marginBottom": "6px"},
                                    ),
                                    html.Pre(
                                        id="domain-details-text",
                                        style={
                                            "whiteSpace": "pre-wrap",
                                            "margin": 0,
                                            "fontFamily": "monospace",
                                            "fontSize": "12px",
                                        },
                                    ),
                                    html.Div(
                                        id="domain-details-link",
                                        style={"marginTop": "2px"},
                                    ),
                                ],
                                style={
                                    "flex": "1",
                                    "padding": "8px",
                                    "border": "1px solid #ddd",
                                    "borderRadius": "4px",
                                    "backgroundColor": "#fafafa",
                                    "display": "none",
                                },
                            ),
                        ],
                        style={
                            "marginTop": "8px",
                            "display": "none",
                            "gap": "10px",
                        },
                    ),
                ], style={"width": "100%", "marginBottom": "18px"}),
                
                html.Hr(),
                
                # 3. Scatter Plot (third)
                html.Div([
                    html.H4(
                        "Summed vs Top Isoform Entropy (colored by min Spearman)",
                        style={**global_style, "marginBottom": "6px", "marginTop": "0px"}
                    ),
                    html.Div([
                        dcc.Checklist(
                            id='negative-spearman-toggle',
                            options=[{'label': ' Show only negative Spearman', 'value': 'negative_only'}],
                            value=[],
                            style={'display': 'inline-block', 'marginRight': '20px'}
                        ),
                    ], style={"marginBottom": "2px", "marginTop": "2px"}),
                    dcc.Graph(id="plot-entropy-scatter"),
                ], style={"width": "100%"}),

                # 4. Co-expression Network Section
                html.Div(
                    [
                        html.Hr(style={"margin": "20px 0"}),
                        html.Div([
                            html.H4("Co-expression Network", style={**global_style, "marginBottom": "10px"}),
                            html.Div(
                                "Click on a gene node to expand it and reveal the internal isoform-to-isoform correlation structure. The network centers around the currently selected gene.",
                                style={"fontSize": "13px", "color": "#666", "marginBottom": "10px"}
                            ),
                            dcc.Loading(
                                id="loading-network",
                                type="circle",
                                children=[
                                    create_coexpression_widget(),
                                    html.Div([
                                        html.Label("Correlation Threshold (|r|): ", style={"fontSize": "13px", "fontWeight": "bold", "marginRight": "10px"}),
                                        dcc.Input(
                                            id="network-threshold-input",
                                            type="number",
                                            min=0.0,
                                            max=1.0,
                                            step=0.01,
                                            value=0.3,
                                            debounce=True,
                                            style={"width": "80px", "padding": "4px", "borderRadius": "4px", "border": "1px solid #ccc"}
                                        ),
                                         html.Label("Max Gene Neighbors: ", style={"fontSize": "13px", "fontWeight": "bold", "marginLeft": "20px", "marginRight": "10px"}),
                                        dcc.Input(
                                            id="network-neighbors-input",
                                            type="number",
                                            min=1,
                                            max=50,
                                            step=1,
                                            value=10,
                                            debounce=True,
                                            style={"width": "60px", "padding": "4px", "borderRadius": "4px", "border": "1px solid #ccc"}
                                        ),
                                        html.Span(" (Press Enter to apply)", style={"fontSize": "11px", "color": "#888", "marginLeft": "10px"})
                                    ], style={"marginTop": "10px", "display": "flex", "alignItems": "center"})
                                ]
                            )
                        ], style={"width": "100%", "marginTop": "20px"}),
                    ],
                    id="coexpression-network-container",
                    style={"display": "block" if gene_coexpression is not None else "none"}
                ),
            ], style={"width": "60%", "display": "inline-block", "verticalAlign": "top", "paddingRight": "20px"}),
            
            # Right panel - Protein Sequence, 3D Structure and Data Table
            html.Div([
                # Protein Sequence Display (at top)
                html.Div([
                    html.H4("Protein Sequence", style={**global_style, "marginTop": "0", "marginBottom": "10px"}),
                    html.Div(id="protein-sequence-container", children=[
                        html.P("Select a transcript to view its protein sequence",
                               style={"textAlign": "center", "color": "#999", "padding": "20px"})
                    ], style={
                        "border": "1px solid #ddd",
                        "borderRadius": "4px",
                        "padding": "10px",
                        "maxHeight": "200px",
                        "overflowY": "auto",
                        "backgroundColor": "#f9f9f9",
                        "fontFamily": "monospace",
                        "fontSize": "12px",
                    })
                ], style={"marginBottom": "30px"}),

                # 3D Structure Viewer (below protein sequence)
                html.Div([
                    html.H4("3D Structure", style={**global_style, "marginBottom": "8px"}),
                    # Exon selector — shown once a transcript with AlphaFold geometry is selected
                    html.Div(
                        id="exon-selector-bar",
                        children=[],
                        style={"marginBottom": "8px", "display": "none"},
                    ),
                    html.Div(id="glb-viewer-container", children=[
                        html.P("Select a transcript with a 3D model to view it here",
                               style={"textAlign": "center", "color": "#999", "padding": "20px"})
                    ], style={
                        "border": "1px solid #ddd",
                        "borderRadius": "4px",
                        "padding": "10px",
                        "minHeight": "400px",
                        "backgroundColor": "#f9f9f9"
                    })
                ], style={"marginBottom": "30px"}),

                # Gene Rankings & Correlations Table (below)
                html.H4("Gene Rankings & Correlations", style={**global_style, "marginBottom": "15px"}),
                
                # Add CDS, 3D structure, and InterPro domain filter toggles
                html.Div([
                    dcc.Checklist(
                        id='cds-filter-toggle',
                        options=[{'label': ' Genes with CDS', 'value': 'filter_cds'}],
                        value=[],
                        style={'display': 'inline-block', 'marginRight': '20px'}
                    ),
                    dcc.Checklist(
                        id='3d-filter-toggle',
                        options=[{'label': ' Genes with 3D structure', 'value': 'filter_3d'}],
                        value=[],
                        style={'display': 'inline-block', 'marginRight': '20px'}
                    ),
                    dcc.Checklist(
                        id='domain-filter-toggle',
                        options=[{'label': ' Genes with domains', 'value': 'filter_domains'}],
                        value=[],
                        style={'display': 'inline-block'}
                    )
                ], style={"marginBottom": "15px"}),

                dash_table.DataTable(
                    id='gene-table',
                    columns=[],  # Will be populated by callback
                    data=[],  # Will be populated by callback
                    style_table={
                        'overflowY': 'auto',
                        'border': '1px solid #ddd'
                    },
                    style_cell={
                        'textAlign': 'left',
                        'padding': '10px',
                        'fontFamily': 'Arial, sans-serif',
                        'fontSize': '13px',
                        'minWidth': '70px',
                        'maxWidth': '160px',
                        'whiteSpace': 'normal',
                        'height': 'auto'
                    },
                    style_cell_conditional=[
                        {
                            'if': {'column_id': 'Gene ID'},
                            'cursor': 'pointer'
                        },
                        {
                            'if': {'column_id': 'Gene Name'},
                            'minWidth': '80px',
                            'maxWidth': '140px',
                            'width': '120px'
                        },
                        {
                            'if': {'column_id': ['Has CDS', 'Has 3D', 'Has Domains']},
                            'minWidth': '60px',
                            'maxWidth': '70px',
                            'width': '70px',
                            'textAlign': 'center'
                        },
                        {
                            'if': {'column_id': ['Rank', '# Isoforms', 'Min Spearman', 'Top Entropy', 'Summed Entropy']},
                            'minWidth': '30px',
                            'maxWidth': '60px',
                            'width': '60px'
                        }
                    ],
                    style_header={
                        'backgroundColor': '#f0f0f0',
                        'fontWeight': 'bold',
                        'borderBottom': '2px solid #333'
                    },
                    style_data_conditional=[
                        {
                            'if': {'row_index': 'odd'},
                            'backgroundColor': '#f9f9f9'
                        },
                        {
                            'if': {'state': 'selected'},
                            'backgroundColor': '#FFD700',
                            'border': '1px solid gold'
                        }
                    ],
                    row_selectable='single',
                    selected_rows=[],
                    page_action='native',
                    page_current=0,
                    page_size=20,
                    filter_action="native",
                    sort_action="native",
                    sort_mode="multi",
                    editable=False,
                    cell_selectable=True,
                    markdown_options={"link_target": "_blank"}
                ),
            ], style={"width": "38%", "display": "inline-block", "verticalAlign": "top", "paddingLeft": "20px", "borderLeft": "2px solid #ddd"}),
        ]),
        
        dcc.Store(id="selected-gene"),
        dcc.Store(id="expanded-network-genes", data=[]),
        dcc.Store(id="expanded-gene-position", data=None),
        dcc.Store(id="selected-transcript"),
        dcc.Store(id="selected-exon"),
    ], style={"margin": "20px", **global_style})


    # Register all callbacks
    _register_callbacks(app, state)

    return app


def _register_callbacks(app, state):
    """Register all dashboard callbacks.

    Args:
        app: Dash app instance
        isoforms_by_gene: Dictionary of exon structures
        df_mean: Mean expression DataFrame (server-side)
        df_sum: Sum expression DataFrame (server-side)
        results_df_mean: Results DataFrame for mean
        results_df_sum: Results DataFrame for sum
        table_df_mean: Pre-computed table data for mean dataset (server-side)
        table_df_sum: Pre-computed table data for sum dataset (server-side)
        ci_df: Confidence intervals DataFrame
        ci_columns: List of CI column names
        global_col_mean: Global column for mean data
        global_col_sum: Global column for sum data
        sample_cols: List of sample columns
        af_geometry_mapping: Output of build_alphafold_geometry_mapping() (optional).
        protein_sequences: Dict mapping transcript_id -> amino-acid sequence (optional).
        domain_mapping: Dict mapping transcript_id -> list of domain dicts (optional).
    """
    # ----------------------------------------------------
    # Dialog Box & File Selector Navigation Callbacks
    # ----------------------------------------------------
    
    def bg_load_data(path_mean, path_sum, path_gtf, path_geom, path_coexp, path_fasta, path_interpro, path_ci):
        def get_path_mtime(p_target):
            if not p_target or not os.path.exists(p_target):
                return None
            try:
                if os.path.isdir(p_target):
                    # Check top-level elements modification times to prevent recursive os.walk freezes
                    mtime = os.path.getmtime(p_target)
                    try:
                        mtime += len(os.listdir(p_target))
                    except:
                        pass
                    return mtime
                else:
                    return os.path.getmtime(p_target)
            except:
                return None

        try:
            # Step 1: Mean Expression TSV loading
            mtime_mean = get_path_mtime(path_mean)
            if (path_mean == state.get('path_mean') and 
                mtime_mean is not None and mtime_mean == state.get('mtime_mean') and 
                state.get('df_mean') is not None):
                state['loading_progress'].update({'step': 2, 'msg': 'Step 2/8: Reusing cached Mean Expression data...'})
                df_mean = state['df_mean']
                results_df_mean = state['results_df_mean']
                sample_cols = state['sample_cols']
                global_col_mean = state['global_col_mean']
            else:
                state['loading_progress'].update({'step': 1, 'msg': 'Step 1/8: Reading Mean Expression TSV...'})
                if not path_mean:
                    raise ValueError("Mean Expression TSV file is required.")
                if not os.path.exists(path_mean):
                    raise FileNotFoundError(f"Mean Expression file not found: {path_mean}")
                if os.path.isdir(path_mean):
                    raise IsADirectoryError(f"Mean Expression path is a directory: {path_mean}")
                    
                df_mean = pd.read_csv(path_mean, sep="\t")
                if "gene_id" not in df_mean.columns or "transcript_id" not in df_mean.columns:
                    raise ValueError("Mean Expression input must contain 'gene_id' and 'transcript_id' columns.")
                if len(df_mean.columns) < 3:
                    raise ValueError("Mean Expression input must have at least 3 columns.")
                    
                global_col_mean = df_mean.columns[2]
                meta_cols_mean = {"gene_id", "transcript_id", global_col_mean}
                sample_cols = [c for c in df_mean.columns if c not in meta_cols_mean and pd.api.types.is_numeric_dtype(df_mean[c])]
                if len(sample_cols) < 2:
                    raise ValueError("Need at least two numeric sample columns to compute correlations.")
                    
                # Step 2: Mean Expression Calculations
                state['loading_progress'].update({'step': 2, 'msg': 'Step 2/8: Calculating Mean Expression Entropy & Correlations...'})
                results_mean = calculate_entropy_and_correlation(df_mean, sample_cols, global_col_mean)
                results_df_mean = pd.DataFrame(results_mean)
                
            # Step 3: Sum Expression
            mtime_sum = get_path_mtime(path_sum) if path_sum else None
            if (path_sum == state.get('path_sum') and 
                (not path_sum or (mtime_sum is not None and mtime_sum == state.get('mtime_sum'))) and 
                'df_sum' in state):
                df_sum = state['df_sum']
                results_df_sum = state['results_df_sum']
                global_col_sum = state['global_col_sum']
                has_sum = state['has_sum']
            else:
                state['loading_progress'].update({'step': 3, 'msg': 'Step 3/8: Processing Sum Expression TSV...'})
                df_sum = None
                results_df_sum = pd.DataFrame()
                global_col_sum = ""
                has_sum = False
                
                if path_sum:
                    if not os.path.exists(path_sum):
                        raise FileNotFoundError(f"Sum Expression file not found: {path_sum}")
                    if os.path.isdir(path_sum):
                        raise IsADirectoryError(f"Sum Expression path is a directory: {path_sum}")
                    df_sum = pd.read_csv(path_sum, sep="\t")
                    if "gene_id" not in df_sum.columns or "transcript_id" not in df_sum.columns:
                        raise ValueError("Sum Expression input must contain 'gene_id' and 'transcript_id' columns.")
                    has_sum = True
                    global_col_sum = df_sum.columns[2]
                    results_sum = calculate_entropy_and_correlation(df_sum, sample_cols, global_col_sum)
                    results_df_sum = pd.DataFrame(results_sum)
                else:
                    global_col_sum = global_col_mean
                    
            # Step 4: Exons GTF
            mtime_gtf = get_path_mtime(path_gtf) if path_gtf else None
            if (path_gtf == state.get('path_gtf') and 
                (not path_gtf or (mtime_gtf is not None and mtime_gtf == state.get('mtime_gtf'))) and 
                state.get('isoforms_by_gene') is not None):
                isoforms_by_gene = state['isoforms_by_gene']
                gene_names = state['gene_names']
            else:
                state['loading_progress'].update({'step': 4, 'msg': 'Step 4/8: Parsing GTF Exon Positions & Gene Names...'})
                isoforms_by_gene = {}
                gene_names = {}
                if path_gtf:
                    if not os.path.exists(path_gtf):
                        raise FileNotFoundError(f"GTF file not found: {path_gtf}")
                    if os.path.isdir(path_gtf):
                        raise IsADirectoryError(f"GTF path is a directory: {path_gtf}")
                    from .gtf_parser import parse_isoform_file, parse_gene_names
                    isoforms_by_gene = parse_isoform_file(path_gtf)
                    gene_names = parse_gene_names(path_gtf)
                    
            # Step 5: AlphaFold Geometry
            mtime_geom = get_path_mtime(path_geom) if path_geom else None
            if (path_geom == state.get('path_geom') and 
                (not path_geom or (mtime_geom is not None and mtime_geom == state.get('mtime_geom'))) and 
                state.get('af_geometry_mapping') is not None):
                af_geometry_mapping = state['af_geometry_mapping']
            else:
                state['loading_progress'].update({'step': 5, 'msg': 'Step 5/8: Loading AlphaFold 3D Structures...'})
                af_geometry_mapping = {}
                if path_geom:
                    if not os.path.exists(path_geom):
                        raise FileNotFoundError(f"Geometry directory not found: {path_geom}")
                    if os.path.isfile(path_geom):
                        raise ValueError(f"Geometry path is a file: {path_geom}")
                    from .alphafold_geometry import build_alphafold_geometry_mapping, discover_exon_viewers
                    af_geometry_mapping = build_alphafold_geometry_mapping(path_geom)
                    if af_geometry_mapping:
                        discover_exon_viewers(af_geometry_mapping)
                        
            # Step 6: Coexpression Matrices
            mtime_coexp = get_path_mtime(path_coexp) if path_coexp else None
            if (path_coexp == state.get('path_coexp') and 
                (not path_coexp or (mtime_coexp is not None and mtime_coexp == state.get('mtime_coexp'))) and 
                'gene_coexpression' in state):
                gene_coexpression = state['gene_coexpression']
                gene_coexpression_idx = state['gene_coexpression_idx']
                isoform_coexpression = state['isoform_coexpression']
                isoform_coexpression_idx = state['isoform_coexpression_idx']
            else:
                state['loading_progress'].update({'step': 6, 'msg': 'Step 6/8: Loading Coexpression Sparse Matrices...'})
                gene_coexpression = None
                gene_coexpression_idx = None
                isoform_coexpression = None
                isoform_coexpression_idx = None
                if path_coexp:
                    if not os.path.exists(path_coexp):
                        raise FileNotFoundError(f"Coexpression directory not found: {path_coexp}")
                    if os.path.isfile(path_coexp):
                        raise ValueError(f"Coexpression path is a file: {path_coexp}")
                    import pickle
                    from scipy.sparse import load_npz
                    gene_mat_path = os.path.join(path_coexp, "gene_coexpression.npz")
                    gene_idx_path = os.path.join(path_coexp, "gene_index.pkl")
                    if os.path.exists(gene_mat_path) and os.path.exists(gene_idx_path):
                        gene_coexpression = load_npz(gene_mat_path)
                        with open(gene_idx_path, 'rb') as f:
                           gene_coexpression_idx = pickle.load(f)
                           
                    iso_mat_path = os.path.join(path_coexp, "isoform_coexpression.npz")
                    iso_idx_path = os.path.join(path_coexp, "isoform_index.pkl")
                    if os.path.exists(iso_mat_path) and os.path.exists(iso_idx_path):
                        isoform_coexpression = load_npz(iso_mat_path)
                        with open(iso_idx_path, 'rb') as f:
                           isoform_coexpression_idx = pickle.load(f)
                           
            # Step 7: FASTA and InterPro Domains
            mtime_fasta = get_path_mtime(path_fasta) if path_fasta else None
            if (path_fasta == state.get('path_fasta') and 
                (not path_fasta or (mtime_fasta is not None and mtime_fasta == state.get('mtime_fasta'))) and 
                state.get('protein_sequences') is not None):
                protein_sequences = state['protein_sequences']
            else:
                state['loading_progress'].update({'step': 7, 'msg': 'Step 7/8: Loading FASTA Sequences & InterPro Domains...'})
                protein_sequences = {}
                if path_fasta:
                    if not os.path.exists(path_fasta):
                        raise FileNotFoundError(f"FASTA file not found: {path_fasta}")
                    if os.path.isdir(path_fasta):
                        raise IsADirectoryError(f"FASTA path is a directory: {path_fasta}")
                    from .dashboard_app import load_protein_sequences
                    protein_sequences = load_protein_sequences(path_fasta)
                    
            mtime_interpro = get_path_mtime(path_interpro) if path_interpro else None
            if (path_interpro == state.get('path_interpro') and 
                path_gtf == state.get('path_gtf') and 
                (not path_interpro or (mtime_interpro is not None and mtime_interpro == state.get('mtime_interpro'))) and 
                (not path_gtf or (mtime_gtf is not None and mtime_gtf == state.get('mtime_gtf'))) and 
                state.get('domain_mapping') is not None):
                domain_mapping = state['domain_mapping']
            else:
                domain_mapping = {}
                if path_interpro:
                    if not os.path.exists(path_interpro):
                        raise FileNotFoundError(f"InterPro directory not found: {path_interpro}")
                    if os.path.isfile(path_interpro):
                        raise ValueError(f"InterPro path is a file: {path_interpro}")
                    from .interpro_parser import build_domain_mapping
                    domain_mapping = build_domain_mapping(
                        path_interpro,
                        isoforms_by_gene,
                        min_evalue=1e-5
                    )
                    
            reused_geom = (path_geom == state.get('path_geom') and (not path_geom or (mtime_geom is not None and mtime_geom == state.get('mtime_geom'))))
            reused_fasta = (path_fasta == state.get('path_fasta') and (not path_fasta or (mtime_fasta is not None and mtime_fasta == state.get('mtime_fasta'))))
            if not (reused_geom and reused_fasta) and af_geometry_mapping and protein_sequences:
                from .alphafold_geometry import extend_geometry_mapping_by_sequence
                extend_geometry_mapping_by_sequence(af_geometry_mapping, protein_sequences)
                
            # Step 8: Bootstrap Confidence Intervals
            mtime_ci = get_path_mtime(path_ci) if path_ci else None
            if (path_ci == state.get('path_ci') and 
                (not path_ci or (mtime_ci is not None and mtime_ci == state.get('mtime_ci'))) and 
                'ci_df' in state):
                ci_df = state['ci_df']
                ci_columns = state['ci_columns']
            else:
                state['loading_progress'].update({'step': 8, 'msg': 'Step 8/8: Loading Confidence Intervals TSV...'})
                ci_df = None
                ci_columns = []
                if path_ci:
                    if not os.path.exists(path_ci):
                        raise FileNotFoundError(f"Confidence Intervals file not found: {path_ci}")
                    if os.path.isdir(path_ci):
                        raise IsADirectoryError(f"Confidence Intervals path is a directory: {path_ci}")
                    ci_df = pd.read_csv(path_ci, sep="\t")
                    if 'isoform' in ci_df.columns:
                        ci_df = ci_df.set_index('isoform')
                    elif ci_df.index.name != 'isoform':
                        ci_df.index.name = 'isoform'
                    ci_columns = [col for col in ci_df.columns if col.startswith('ci_')]
                    
            # Apply loaded variables to state
            state['df_mean'] = df_mean
            state['df_sum'] = df_sum
            state['results_df_mean'] = results_df_mean
            state['results_df_sum'] = results_df_sum
            state['sample_cols'] = sample_cols
            state['global_col_mean'] = global_col_mean
            state['global_col_sum'] = global_col_sum
            state['has_sum'] = has_sum
            state['isoforms_by_gene'] = isoforms_by_gene
            state['gene_names'] = gene_names
            state['af_geometry_mapping'] = af_geometry_mapping
            state['domain_mapping'] = domain_mapping
            state['protein_sequences'] = protein_sequences
            state['gene_coexpression'] = gene_coexpression
            state['gene_coexpression_idx'] = gene_coexpression_idx
            state['isoform_coexpression'] = isoform_coexpression
            state['isoform_coexpression_idx'] = isoform_coexpression_idx
            state['ci_df'] = ci_df
            state['ci_columns'] = ci_columns or []
            
            state['path_mean'] = path_mean
            state['path_sum'] = path_sum
            state['path_gtf'] = path_gtf
            state['path_geom'] = path_geom
            state['path_coexp'] = path_coexp
            state['path_fasta'] = path_fasta
            state['path_interpro'] = path_interpro
            state['path_ci'] = path_ci
            
            state['mtime_mean'] = mtime_mean
            state['mtime_sum'] = mtime_sum
            state['mtime_gtf'] = mtime_gtf
            state['mtime_geom'] = mtime_geom
            state['mtime_coexp'] = mtime_coexp
            state['mtime_fasta'] = mtime_fasta
            state['mtime_interpro'] = mtime_interpro
            state['mtime_ci'] = mtime_ci
            state['path_ci'] = path_ci
            
            _base_scatter.cache_clear()
            _base_exon_fig.cache_clear()
            
            # Perform derived calculation
            state['recompute_derived_state']()
            
            # Complete
            state['loading_progress'].update({'done': True})
        except Exception as e:
            import traceback
            traceback.print_exc()
            state['loading_progress'].update({'error': str(e)})
    
    @app.callback(
        Output("settings-popup-window", "style"),
        [Input("btn-open-settings", "n_clicks"),
         Input("btn-close-settings", "n_clicks"),
         Input("btn-cancel-settings", "n_clicks"),
         Input("data-sources-updated", "data")],
        [State("settings-popup-window", "style")]
    )
    def toggle_settings(open_clicks, close_clicks, cancel_clicks, updated_time, current_style):
        trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]
        if not current_style:
            current_style = {"display": "none"}
        else:
            current_style = current_style.copy()
            
        if trigger_id == "btn-open-settings":
            current_style["display"] = "block"
        elif trigger_id in ["btn-close-settings", "btn-cancel-settings", "data-sources-updated"]:
            current_style["display"] = "none"
        return current_style

    @app.callback(
        Output("preprocess-popup-window", "style"),
        [Input("btn-open-preprocess", "n_clicks"),
         Input("btn-close-preprocess", "n_clicks"),
         Input("btn-cancel-preprocess", "n_clicks")],
        [State("preprocess-popup-window", "style")]
    )
    def toggle_preprocess(open_clicks, close_clicks, cancel_clicks, current_style):
        trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]
        if not current_style:
            current_style = {"display": "none"}
        else:
            current_style = current_style.copy()
            
        if trigger_id == "btn-open-preprocess":
            current_style["display"] = "block"
        elif trigger_id in ["btn-close-preprocess", "btn-cancel-preprocess"]:
            current_style["display"] = "none"
        return current_style

    @app.callback(
        [Output("tool-config-dist", "style"),
         Output("tool-config-coexp", "style"),
         Output("tool-config-interpro", "style"),
         Output("tool-config-map-drive", "style")],
        [Input("dropdown-select-tool", "value")]
    )
    def toggle_preprocess_tool_config(selected_tool):
        if selected_tool == 'dist':
            return {"display": "block"}, {"display": "none"}, {"display": "none"}, {"display": "none"}
        elif selected_tool == 'coexp':
            return {"display": "none"}, {"display": "block"}, {"display": "none"}, {"display": "none"}
        elif selected_tool == 'interpro':
            return {"display": "none"}, {"display": "none"}, {"display": "block"}, {"display": "none"}
        else: # map-drive
            return {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "block"}

    @app.callback(
        Output("tool-dist-bootstrap-params", "style"),
        [Input("tool-dist-run-bootstrap", "value")]
    )
    def toggle_bootstrap_params(run_bootstrap_value):
        if run_bootstrap_value and 'run_ci' in run_bootstrap_value:
            return {"display": "block", "paddingLeft": "20px", "borderLeft": "3px solid #eee", "marginTop": "10px"}
        return {"display": "none"}

    @app.callback(
        [Output("tool-dist-matrix", "value"),
         Output("tool-dist-gtf", "value"),
         Output("tool-coexp-expression", "value"),
         Output("tool-coexp-gtf", "value"),
         Output("tool-interpro-fasta", "value"),
         Output("tool-interpro-gtf", "value"),
         Output("tool-dist-outdir", "value"),
         Output("tool-coexp-outdir", "value"),
         Output("tool-interpro-outdir", "value")],
        [Input("btn-open-preprocess", "n_clicks")],
        [State("input-path-mean", "value"),
         State("input-path-gtf", "value"),
         State("input-path-fasta", "value"),
         State("input-path-interpro", "value")],
        prevent_initial_call=True
    )
    def prepopulate_preprocess_paths(n_clicks, path_mean, path_gtf, path_fasta, path_interpro):
        if not n_clicks:
            raise PreventUpdate
        
        default_dist_out = "data/isoform_distributions"
        default_coexp_out = "data/co-expression"
        default_interpro_out = path_interpro if path_interpro else "data/interpro_results"
        
        return path_mean, path_gtf, path_mean, path_gtf, path_fasta, path_gtf, default_dist_out, default_coexp_out, default_interpro_out

    @app.callback(
        [Output("file-selector-popup", "style"),
         Output("file-browser-target-field", "data"),
         Output("file-browser-current-dir", "data"),
         Output("file-selector-selected-input", "value"),
         Output("file-selector-title", "children")],
        [Input("btn-browse-mean", "n_clicks"),
         Input("btn-browse-sum", "n_clicks"),
         Input("btn-browse-gtf", "n_clicks"),
         Input("btn-browse-geom", "n_clicks"),
         Input("btn-browse-coexp", "n_clicks"),
         Input("btn-browse-fasta", "n_clicks"),
         Input("btn-browse-interpro", "n_clicks"),
         Input("btn-browse-ci", "n_clicks"),
         Input("btn-browse-tool-dist-matrix", "n_clicks"),
         Input("btn-browse-tool-dist-gtf", "n_clicks"),
         Input("btn-browse-tool-dist-meta", "n_clicks"),
         Input("btn-browse-tool-dist-outdir", "n_clicks"),
         Input("btn-browse-tool-coexp-expr", "n_clicks"),
         Input("btn-browse-tool-coexp-gtf", "n_clicks"),
         Input("btn-browse-tool-coexp-outdir", "n_clicks"),
         Input("btn-browse-tool-interpro-fasta", "n_clicks"),
         Input("btn-browse-tool-interpro-gtf", "n_clicks"),
         Input("btn-browse-tool-interpro-outdir", "n_clicks"),
         Input("btn-close-file-selector", "n_clicks"),
         Input("btn-cancel-file-selection", "n_clicks"),
         Input("btn-confirm-file-selection", "n_clicks")],
        [State("input-path-mean", "value"),
         State("input-path-sum", "value"),
         State("input-path-gtf", "value"),
         State("input-path-geom", "value"),
         State("input-path-coexp", "value"),
         State("input-path-fasta", "value"),
         State("input-path-interpro", "value"),
         State("input-path-ci", "value"),
         State("tool-dist-matrix", "value"),
         State("tool-dist-gtf", "value"),
         State("tool-dist-meta", "value"),
         State("tool-dist-outdir", "value"),
         State("tool-coexp-expression", "value"),
         State("tool-coexp-gtf", "value"),
         State("tool-coexp-outdir", "value"),
         State("tool-interpro-fasta", "value"),
         State("tool-interpro-gtf", "value"),
         State("tool-interpro-outdir", "value"),
         State("file-browser-current-dir", "data"),
         State("file-selector-selected-input", "value"),
         State("file-browser-target-field", "data"),
         State("file-selector-popup", "style")]
    )
    def handle_file_selector(browse_mean, browse_sum, browse_gtf, browse_geom, browse_coexp, browse_fasta, browse_interpro, browse_ci,
                             browse_t_dist_matrix, browse_t_dist_gtf, browse_t_dist_meta, browse_t_dist_outdir,
                             browse_t_coexp_expr, browse_t_coexp_gtf, browse_t_coexp_outdir,
                             browse_t_interpro_fasta, browse_t_interpro_gtf, browse_t_interpro_outdir,
                             close_clicks, cancel_clicks, confirm_clicks,
                             path_mean, path_sum, path_gtf, path_geom, path_coexp, path_fasta, path_interpro, path_ci,
                             path_t_dist_matrix, path_t_dist_gtf, path_t_dist_meta, path_t_dist_outdir,
                             path_t_coexp_expr, path_t_coexp_gtf, path_t_coexp_outdir,
                             path_t_interpro_fasta, path_t_interpro_gtf, path_t_interpro_outdir,
                             current_dir, selected_path, target_field, current_style):
        trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]
        if not current_style:
            current_style = {"display": "none"}
        else:
            current_style = current_style.copy()
            
        if not current_dir:
            current_dir = os.getcwd()
            
        if trigger_id.startswith("btn-browse-"):
            field = trigger_id.replace("btn-browse-", "")
            current_style["display"] = "block"
            
            field_paths = {
                "mean": path_mean, "sum": path_sum, "gtf": path_gtf,
                "geom": path_geom, "coexp": path_coexp, "fasta": path_fasta, "interpro": path_interpro, "ci": path_ci,
                "tool-dist-matrix": path_t_dist_matrix,
                "tool-dist-gtf": path_t_dist_gtf,
                "tool-dist-meta": path_t_dist_meta,
                "tool-dist-outdir": path_t_dist_outdir,
                "tool-coexp-expr": path_t_coexp_expr,
                "tool-coexp-gtf": path_t_coexp_gtf,
                "tool-coexp-outdir": path_t_coexp_outdir,
                "tool-interpro-fasta": path_t_interpro_fasta,
                "tool-interpro-gtf": path_t_interpro_gtf,
                "tool-interpro-outdir": path_t_interpro_outdir
            }
            init_path = field_paths.get(field) or ""
            init_path = init_path.strip()
            
            if init_path and os.path.exists(init_path):
                if os.path.isdir(init_path):
                    current_dir = os.path.abspath(init_path)
                else:
                    current_dir = os.path.dirname(os.path.abspath(init_path))
                    
            title_mappings = {
                "mean": "Select Mean Expression TSV File",
                "sum": "Select Sum Expression TSV File",
                "gtf": "Select Exons GTF File",
                "geom": "Select AlphaFold Geometry Directory",
                "coexp": "Select Precomputed Co-expression Directory",
                "fasta": "Select Amino Acid FASTA File",
                "interpro": "Select Precomputed InterPro TSV/XML/JSON Directory",
                "ci": "Select Bootstrap Confidence Intervals File",
                "tool-dist-matrix": "Preprocess: Select Raw Expression Matrix (TSV/CSV) File",
                "tool-dist-gtf": "Preprocess: Select Exons GTF File",
                "tool-dist-meta": "Preprocess: Select Sample Metadata File",
                "tool-dist-outdir": "Preprocess: Select Distribution Output Directory",
                "tool-coexp-expr": "Preprocess: Select Raw Expression Matrix (TSV/CSV) File",
                "tool-coexp-gtf": "Preprocess: Select Exons GTF File",
                "tool-coexp-outdir": "Preprocess: Select Co-expression Output Directory",
                "tool-interpro-fasta": "Preprocess: Select FASTA Sequences File",
                "tool-interpro-gtf": "Preprocess: Select Exons GTF File (Optional filter)",
                "tool-interpro-outdir": "Preprocess: Select InterPro Scan Output Directory"
            }
            title = title_mappings.get(field, f"Select {field.upper()} File")
            return current_style, field, current_dir, init_path, title
            
        elif trigger_id in ["btn-close-file-selector", "btn-cancel-file-selection", "btn-confirm-file-selection"]:
            current_style["display"] = "none"
            
        return current_style, no_update, no_update, no_update, no_update

    @app.callback(
        [Output("input-path-mean", "value"),
         Output("input-path-sum", "value"),
         Output("input-path-gtf", "value"),
         Output("input-path-geom", "value"),
         Output("input-path-coexp", "value"),
         Output("input-path-fasta", "value"),
         Output("input-path-interpro", "value"),
         Output("input-path-ci", "value"),
         Output("tool-dist-matrix", "value"),
         Output("tool-dist-gtf", "value"),
         Output("tool-dist-meta", "value"),
         Output("tool-dist-outdir", "value"),
         Output("tool-coexp-expression", "value"),
         Output("tool-coexp-gtf", "value"),
         Output("tool-coexp-outdir", "value"),
         Output("tool-interpro-fasta", "value"),
         Output("tool-interpro-gtf", "value"),
         Output("tool-interpro-outdir", "value")],
        [Input("btn-confirm-file-selection", "n_clicks")],
        [State("file-selector-selected-input", "value"),
         State("file-browser-target-field", "data"),
         State("input-path-mean", "value"),
         State("input-path-sum", "value"),
         State("input-path-gtf", "value"),
         State("input-path-geom", "value"),
         State("input-path-coexp", "value"),
         State("input-path-fasta", "value"),
         State("input-path-interpro", "value"),
         State("input-path-ci", "value"),
         State("tool-dist-matrix", "value"),
         State("tool-dist-gtf", "value"),
         State("tool-dist-meta", "value"),
         State("tool-dist-outdir", "value"),
         State("tool-coexp-expression", "value"),
         State("tool-coexp-gtf", "value"),
         State("tool-coexp-outdir", "value"),
         State("tool-interpro-fasta", "value"),
         State("tool-interpro-gtf", "value"),
         State("tool-interpro-outdir", "value")]
    )
    def update_settings_input(confirm_clicks, selected_path, target_field,
                              val_mean, val_sum, val_gtf, val_geom, val_coexp, val_fasta, val_interpro, val_ci,
                              val_t_dist_matrix, val_t_dist_gtf, val_t_dist_meta, val_t_dist_outdir,
                              val_t_coexp_expr, val_t_coexp_gtf, val_t_coexp_outdir,
                              val_t_interpro_fasta, val_t_interpro_gtf, val_t_interpro_outdir):
        if not confirm_clicks or not target_field:
            raise PreventUpdate
            
        outputs = [val_mean, val_sum, val_gtf, val_geom, val_coexp, val_fasta, val_interpro, val_ci,
                   val_t_dist_matrix, val_t_dist_gtf, val_t_dist_meta, val_t_dist_outdir,
                   val_t_coexp_expr, val_t_coexp_gtf, val_t_coexp_outdir,
                   val_t_interpro_fasta, val_t_interpro_gtf, val_t_interpro_outdir]
        fields = ["mean", "sum", "gtf", "geom", "coexp", "fasta", "interpro", "ci",
                  "tool-dist-matrix", "tool-dist-gtf", "tool-dist-meta", "tool-dist-outdir",
                  "tool-coexp-expr", "tool-coexp-gtf", "tool-coexp-outdir",
                  "tool-interpro-fasta", "tool-interpro-gtf", "tool-interpro-outdir"]
        
        if target_field in fields:
            idx = fields.index(target_field)
            outputs[idx] = selected_path
            
        return outputs

    @app.callback(
        [Output("tool-console-output", "children"),
         Output("tool-console-container", "style"),
         Output("btn-run-tool", "style"),
         Output("btn-kill-tool", "style"),
         Output("tool-interval", "disabled"),
         Output("tool-progress-container", "style"),
         Output("tool-progress-status-msg", "children"),
         Output("tool-progress-bar-fill", "style")],
        [Input("btn-run-tool", "n_clicks"),
         Input("btn-kill-tool", "n_clicks"),
         Input("tool-interval", "n_intervals")],
        [State("dropdown-select-tool", "value"),
         State("tool-dist-matrix", "value"),
         State("tool-dist-gtf", "value"),
         State("tool-dist-meta", "value"),
         State("tool-dist-outdir", "value"),
         State("tool-dist-meta-sample", "value"),
         State("tool-dist-meta-group", "value"),
         State("tool-dist-cutoff", "value"),
         State("tool-dist-run-bootstrap", "value"),
         State("tool-dist-bootstrap-iter", "value"),
         State("tool-dist-bootstrap-key", "value"),
         State("tool-dist-bootstrap-seed", "value"),
         State("tool-coexp-expression", "value"),
         State("tool-coexp-gtf", "value"),
         State("tool-coexp-outdir", "value"),
         State("tool-coexp-threshold", "value"),
         State("tool-coexp-chunk", "value"),
         State("tool-coexp-method", "value"),
         State("tool-interpro-fasta", "value"),
         State("tool-interpro-gtf", "value"),
         State("tool-interpro-outdir", "value"),
         State("tool-interpro-email", "value"),
         State("tool-interpro-skip", "value"),
         State("tool-map-share", "value"),
         State("tool-map-letter", "value"),
         State("tool-map-user", "value"),
         State("tool-map-pass", "value")],
        prevent_initial_call=True
    )
    def manage_tool_execution(run_clicks, kill_clicks, n_intervals,
                              selected_tool, dist_matrix, dist_gtf, dist_meta, dist_outdir,
                              dist_sample, dist_group, dist_cutoff, dist_run_bootstrap,
                              bootstrap_iter, bootstrap_key, bootstrap_seed,
                              coexp_expr, coexp_gtf, coexp_outdir, coexp_threshold, coexp_chunk, coexp_method,
                              interpro_fasta, interpro_gtf, interpro_outdir, interpro_email, interpro_skip,
                              map_share, map_letter, map_user, map_pass):
        import sys
        import subprocess
        trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]
        log_path = os.path.join(os.getcwd(), "isoform_dashboard", "assets", "tool_run.log")
        os.makedirs(os.path.dirname(log_path), exist_ok=True)

        style_run_enabled = {
            "backgroundColor": "#2ecc71", "color": "white", "border": "none",
            "padding": "10px 20px", "borderRadius": "4px", "cursor": "pointer",
            "fontWeight": "bold", "marginRight": "10px", "fontSize": "14px"
        }
        style_run_disabled = {
            "backgroundColor": "#95a5a6", "color": "#e0e0e0", "border": "none",
            "padding": "10px 20px", "borderRadius": "4px", "cursor": "not-allowed",
            "fontWeight": "bold", "marginRight": "10px", "fontSize": "14px"
        }
        style_kill_visible = {
            "backgroundColor": "#e74c3c", "color": "white", "border": "none",
            "padding": "10px 20px", "borderRadius": "4px", "cursor": "pointer",
            "fontWeight": "bold", "fontSize": "14px", "display": "inline-block"
        }
        style_kill_hidden = {"display": "none"}

        # Helper to read log
        def _read_log_file():
            if os.path.exists(log_path):
                try:
                    with open(log_path, "r", encoding="utf-8", errors="replace") as f:
                        return f.read()
                except Exception as e:
                    return f"Error reading log file: {str(e)}"
            return ""

        # Parse progress
        def _parse_progress_percent(log_text):
            if not log_text:
                return 0
            import re
            matches = re.findall(r'\[PROGRESS_PERCENT\]\s*(\d+)', log_text)
            if not matches:
                return 0
            val = int(matches[-1])
            if selected_tool == 'dist' and dist_run_bootstrap and 'run_ci' in dist_run_bootstrap:
                if "=== Running Step 2/2" in log_text:
                    return int(30 + val * 0.7)
                else:
                    return int(val * 0.3)
            return val

        p = state.get('tool_process')

        if trigger_id == "btn-run-tool" and run_clicks:
            if p is not None and p.poll() is None:
                # Already running
                return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update

            # Construct command
            if selected_tool == 'dist':
                if not dist_matrix or not dist_gtf or not dist_meta or not dist_outdir:
                    return "Error: Raw Expression Matrix, GTF, Metadata, and Output Directory are all required.", {"display": "block"}, style_run_enabled, style_kill_hidden, True, {"display": "none"}, "", {"width": "0%"}
                
                cmd1 = [
                    sys.executable, "-u", "-m", "isoform_distribution.distributions",
                    "--matrix", dist_matrix.strip(),
                    "--gtf", dist_gtf.strip(),
                    "--meta-file", dist_meta.strip(),
                    "--output-dir", dist_outdir.strip(),
                    "--meta-sample-col", (dist_sample or "sample_id").strip(),
                    "--meta-group-col", (dist_group or "region").strip(),
                    "--cutoff-pct", str(dist_cutoff if dist_cutoff is not None else 1.5),
                    "--table-type", "both"
                ]
                
                if dist_run_bootstrap and 'run_ci' in dist_run_bootstrap:
                    cmd2 = [
                        sys.executable, "-u", "-m", "isoform_distribution.bootstrap_isoform_means",
                        "--input", dist_matrix.strip(),
                        "--output-dir", dist_outdir.strip(),
                        "--iterations", str(bootstrap_iter if bootstrap_iter is not None else 1000),
                        "--seed", str(bootstrap_seed if bootstrap_seed is not None else 42),
                        "--include-key", (bootstrap_key or "adult").strip()
                    ]
                    
                    # Write the pipeline runner script to chain them sequentially
                    runner_path = os.path.join(os.getcwd(), "isoform_dashboard", "assets", "pipeline_runner.py")
                    with open(runner_path, "w", encoding="utf-8") as f:
                        f.write(f"""import subprocess
import sys

cmd1 = {repr(cmd1)}
cmd2 = {repr(cmd2)}

print("=== Running Step 1/2: Calculate Distribution Tables ===")
sys.stdout.flush()
res1 = subprocess.run(cmd1)
if res1.returncode != 0:
    print(f"\\nStep 1 failed with exit code {{res1.returncode}}")
    sys.exit(res1.returncode)

print("\\n=== Running Step 2/2: Calculate Bootstrap CIs ===")
sys.stdout.flush()
res2 = subprocess.run(cmd2)
if res2.returncode != 0:
    print(f"\\nStep 2 failed with exit code {{res2.returncode}}")
    sys.exit(res2.returncode)

print("\\n=== Preprocessing Completed Successfully! ===")
""")
                    cmd = [sys.executable, "-u", runner_path]
                else:
                    cmd = cmd1
                    
            elif selected_tool == 'coexp':
                if not coexp_expr or not coexp_gtf or not coexp_outdir:
                    return "Error: Expression TSV, GTF, and Output Directory are all required.", {"display": "block"}, style_run_enabled, style_kill_hidden, True, {"display": "none"}, "", {"width": "0%"}
                cmd = [
                    sys.executable, "-u", "-m", "utilities.precompute_coexpression",
                    "--expression", coexp_expr.strip(),
                    "--gtf", coexp_gtf.strip(),
                    "--outdir", coexp_outdir.strip(),
                    "--threshold", str(coexp_threshold if coexp_threshold is not None else 0.001),
                    "--chunk-size", str(coexp_chunk if coexp_chunk is not None else 250),
                    "--method", (coexp_method or "spearman").strip()
                ]
            elif selected_tool == 'interpro':
                if not interpro_fasta or not interpro_outdir:
                    return "Error: FASTA Sequences File and Output Directory are required.", {"display": "block"}, style_run_enabled, style_kill_hidden, True, {"display": "none"}, "", {"width": "0%"}
                cmd = [
                    sys.executable, "-u", "-m", "interpro.interpro_scan",
                    "--batch",
                    "--fasta", interpro_fasta.strip(),
                    "--output-dir", interpro_outdir.strip(),
                    "--email", (interpro_email or "bm708@cam.ac.uk").strip()
                ]
                if interpro_gtf and interpro_gtf.strip():
                    cmd.extend(["--gtf", interpro_gtf.strip()])
                if interpro_skip and 'skip' in interpro_skip:
                    cmd.append("--skip-existing")
            else: # map-drive
                if sys.platform != 'win32':
                    return "Error: Mapping network drives via the UI is only supported on Windows host servers. For Linux/macOS servers, please use 'mount' or 'sshfs' as described in hosting_instructions.md.", {"display": "block"}, style_run_enabled, style_kill_hidden, True, {"display": "none"}, "", {"width": "0%"}
                if not map_share or not map_letter:
                    return "Error: Remote Share Path and Drive Letter are required.", {"display": "block"}, style_run_enabled, style_kill_hidden, True, {"display": "none"}, "", {"width": "0%"}
                
                share_path = map_share.strip()
                letter = map_letter.strip()
                if not letter.endswith(":"):
                    letter += ":"
                
                cmd_delete = ["net", "use", letter, "/delete", "/y"]
                cmd_map = ["net", "use", letter, share_path]
                if map_pass and map_pass.strip() and map_user and map_user.strip():
                    cmd_map.extend([map_pass.strip(), f"/user:{map_user.strip()}"])
                elif map_user and map_user.strip():
                    cmd_map.extend([f"/user:{map_user.strip()}"])
                cmd_map.append("/persistent:yes")
                
                # Write the runner script
                runner_path = os.path.join(os.getcwd(), "isoform_dashboard", "assets", "drive_mapper_runner.py")
                with open(runner_path, "w", encoding="utf-8") as f:
                    f.write(f"""import subprocess
import sys

print("=== Disconnecting drive {letter} (if already mapped) ===")
sys.stdout.flush()
subprocess.run({repr(cmd_delete)}, shell=True)

print("\\n=== Mapping drive {letter} to {share_path} ===")
sys.stdout.flush()
print("[PROGRESS_PERCENT] 50")
sys.stdout.flush()

res = subprocess.run({repr(cmd_map)}, shell=True)
if res.returncode != 0:
    print(f"\\nMapping failed with exit code {{res.returncode}}")
    sys.exit(res.returncode)

print("[PROGRESS_PERCENT] 100")
print("\\n=== Drive mapped successfully! ===")
""")
                cmd = [sys.executable, "-u", runner_path]

            # Start process
            try:
                log_file = open(log_path, "w", encoding="utf-8")
                log_file.write(f"Executing command: {' '.join(cmd)}\n\n")
                log_file.flush()
                
                # Run the process
                p = subprocess.Popen(
                    cmd,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
                state['tool_process'] = p
                state['tool_log_file'] = log_file
            except Exception as e:
                return f"Failed to start process: {str(e)}", {"display": "block"}, style_run_enabled, style_kill_hidden, True, {"display": "none"}, "", {"width": "0%"}

            return f"Executing command: {' '.join(cmd)}\n...", {"display": "block"}, style_run_disabled, style_kill_visible, False, {"display": "block"}, "Initializing precomputation tool...", {"width": "0%", "height": "100%", "backgroundColor": "#2ecc71", "transition": "width 0.2s"}

        elif trigger_id == "btn-kill-tool" and kill_clicks:
            pct = 0
            if p is not None:
                p.terminate()
                try:
                    p.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    p.kill()
                
                state['tool_process'] = None
                if state.get('tool_log_file'):
                    try:
                        state['tool_log_file'].write("\n[Process terminated by user]\n")
                        state['tool_log_file'].close()
                    except:
                        pass
                    state['tool_log_file'] = None

            log_text = _read_log_file()
            pct = _parse_progress_percent(log_text)
            return log_text, {"display": "block"}, style_run_enabled, style_kill_hidden, True, {"display": "block"}, "Process terminated by user.", {"width": f"{pct}%", "height": "100%", "backgroundColor": "#e74c3c", "transition": "width 0.2s"}

        elif trigger_id == "tool-interval":
            log_text = _read_log_file()
            pct = _parse_progress_percent(log_text)
            if p is not None:
                exit_code = p.poll()
                if exit_code is not None:
                    # Completed
                    state['tool_process'] = None
                    if state.get('tool_log_file'):
                        try:
                            state['tool_log_file'].close()
                        except:
                            pass
                        state['tool_log_file'] = None
                    
                    log_text += f"\n[Process finished with exit code {exit_code}]\n"
                    # Also write it to the log file so it's persisted
                    try:
                        with open(log_path, "a", encoding="utf-8") as lf:
                            lf.write(f"\n[Process finished with exit code {exit_code}]\n")
                    except:
                        pass
                        
                    completed_msg = "Precomputation completed successfully!" if exit_code == 0 else f"Precomputation failed with exit code {exit_code}"
                    bar_color = "#2ecc71" if exit_code == 0 else "#e74c3c"
                    bar_width = "100%" if exit_code == 0 else f"{pct}%"
                    return log_text, {"display": "block"}, style_run_enabled, style_kill_hidden, True, {"display": "block"}, completed_msg, {"width": bar_width, "height": "100%", "backgroundColor": bar_color, "transition": "width 0.2s"}
                else:
                    # Still running
                    return log_text, {"display": "block"}, style_run_disabled, style_kill_visible, False, {"display": "block"}, f"Running... {pct}% completed", {"width": f"{pct}%", "height": "100%", "backgroundColor": "#2ecc71", "transition": "width 0.2s"}
            else:
                return log_text, {"display": "block"}, style_run_enabled, style_kill_hidden, True, {"display": "block"}, f"Idle. Last run: {pct}% completed", {"width": f"{pct}%", "height": "100%", "backgroundColor": "#2ecc71", "transition": "width 0.2s"}

        return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update

    @app.callback(
        Output("file-selector-selected-input", "value", allow_duplicate=True),
        Input("btn-select-current-folder", "n_clicks"),
        State("file-browser-current-dir", "data"),
        prevent_initial_call=True
    )
    def select_current_folder(n_clicks, current_dir):
        if not n_clicks:
            raise PreventUpdate
        return current_dir

    @app.callback(
        [Output("file-browser-list-container", "children"),
         Output("file-browser-breadcrumbs", "children"),
         Output("system-roots-container", "children")],
        [Input("file-browser-current-dir", "data"),
         Input("file-browser-target-field", "data")]
    )
    def render_directory_contents(current_dir, target_field):
        if not current_dir:
            current_dir = os.getcwd()
        
        target_type = 'all'
        if target_field in ["geom", "coexp", "interpro"]:
            target_type = 'dir'
            
        items, error = get_directory_contents(current_dir, target_type)
        
        children = []
        if error:
            children.append(html.Div(f"Error: {error}", style={"color": "red", "padding": "10px"}))
        else:
            for item in items:
                icon = "📁 " if item['is_dir'] else "📄 "
                name = item['name']
                path = item['path']
                
                size_str = ""
                if not item['is_dir'] and item['size'] is not None:
                    size = item['size']
                    if size > 1024 * 1024:
                        size_str = f" ({size / (1024*1024):.1f} MB)"
                    elif size > 1024:
                        size_str = f" ({size / 1024:.1f} KB)"
                    else:
                        size_str = f" ({size} bytes)"
                        
                children.append(
                    html.Div(
                        [
                           html.Span(icon, style={"marginRight": "5px"}),
                           html.Span(name + size_str)
                        ],
                        id={"type": "file-browser-item", "index": path, "is_dir": str(item['is_dir'])},
                        className="file-item",
                        style={"padding": "8px 12px", "borderBottom": "1px solid #eee"}
                    )
                )
                
        breadcrumb_elements = []
        parts = []
        drive, path_tail = os.path.splitdrive(current_dir)
        
        if drive:
            parts.append(drive + "\\")
        elif path_tail.startswith("/"):
            parts.append("/")
            path_tail = path_tail[1:]
            
        if path_tail:
            subparts = [p for p in path_tail.replace("\\", "/").split("/") if p]
            parts.extend(subparts)
            
        accum_path = ""
        for idx, part in enumerate(parts):
            if idx == 0 and drive:
                accum_path = drive + "\\"
            elif idx == 0 and part == "/":
                accum_path = "/"
            else:
                accum_path = os.path.join(accum_path, part)
                
            breadcrumb_elements.append(
                html.Span(
                    part,
                    id={"type": "breadcrumb-segment", "index": accum_path},
                    className="breadcrumb-segment"
                )
            )
            if idx < len(parts) - 1:
                breadcrumb_elements.append(html.Span(" > ", style={"margin": "0 4px", "color": "#7f8c8d"}))
                
        roots = get_system_roots()
        root_buttons = []
        for r in roots:
            root_buttons.append(
                html.Button(
                    r,
                    id={"type": "drive-jump-btn", "index": r},
                    className="drive-btn",
                    n_clicks=0
                )
            )
            
        return children, breadcrumb_elements, root_buttons

    @app.callback(
        [Output("file-browser-current-dir", "data", allow_duplicate=True),
         Output("file-selector-selected-input", "value", allow_duplicate=True)],
        [Input({"type": "file-browser-item", "index": ALL, "is_dir": ALL}, "n_clicks"),
         Input({"type": "breadcrumb-segment", "index": ALL}, "n_clicks"),
         Input({"type": "drive-jump-btn", "index": ALL}, "n_clicks")],
        [State("file-browser-current-dir", "data"),
         State("file-selector-selected-input", "value"),
         State("file-browser-target-field", "data")],
        prevent_initial_call=True
    )
    def navigate_filesystem(item_clicks, breadcrumb_clicks, drive_clicks, current_dir, current_selected, target_field):
        triggered_id = callback_context.triggered_id
        if not triggered_id or not isinstance(triggered_id, dict):
            raise PreventUpdate
            
        target_type = triggered_id.get("type")
        clicked_path = triggered_id.get("index")
        
        if target_type == "file-browser-item":
            is_dir = str(triggered_id.get("is_dir", "False")).lower() == "true"
            if is_dir:
                is_dir_target = target_field in ["geom", "coexp", "interpro"]
                selected_val = clicked_path if is_dir_target else current_selected
                return clicked_path, selected_val
            else:
                return current_dir, clicked_path
                
        elif target_type in ["breadcrumb-segment", "drive-jump-btn"]:
            is_dir_target = target_field in ["geom", "coexp", "interpro"]
            selected_val = clicked_path if is_dir_target else current_selected
            return clicked_path, selected_val
            
        raise PreventUpdate

    @app.callback(
        [Output("progress-interval", "disabled", allow_duplicate=True),
         Output("progress-container", "style", allow_duplicate=True),
         Output("settings-feedback", "style", allow_duplicate=True),
         Output("btn-apply-settings", "disabled"),
         Output("btn-cancel-settings", "disabled"),
         Output("btn-browse-mean", "disabled"),
         Output("btn-browse-sum", "disabled"),
         Output("btn-browse-gtf", "disabled"),
         Output("btn-browse-geom", "disabled"),
         Output("btn-browse-coexp", "disabled"),
         Output("btn-browse-fasta", "disabled"),
         Output("btn-browse-interpro", "disabled"),
         Output("btn-browse-ci", "disabled"),
         Output("btn-open-settings", "disabled"),
         Output("btn-open-preprocess", "disabled")],
        [Input("btn-apply-settings", "n_clicks")],
        [State("input-path-mean", "value"),
         State("input-path-sum", "value"),
         State("input-path-gtf", "value"),
         State("input-path-geom", "value"),
         State("input-path-coexp", "value"),
         State("input-path-fasta", "value"),
         State("input-path-interpro", "value"),
         State("input-path-ci", "value")],
        prevent_initial_call=True
    )
    def apply_data_sources(n_clicks, path_mean, path_sum, path_gtf, path_geom, path_coexp, path_fasta, path_interpro, path_ci):
        if not n_clicks:
            raise PreventUpdate
            
        path_mean = (path_mean or "").strip()
        path_sum = (path_sum or "").strip()
        path_gtf = (path_gtf or "").strip()
        path_geom = (path_geom or "").strip()
        path_coexp = (path_coexp or "").strip()
        path_fasta = (path_fasta or "").strip()
        path_interpro = (path_interpro or "").strip()
        path_ci = (path_ci or "").strip()
        
        state['loading_progress'] = {
            'step': 0,
            'total': 8,
            'msg': 'Initializing loading sequence...',
            'done': False,
            'error': None
        }
        
        import threading
        thread = threading.Thread(
            target=bg_load_data,
            args=(path_mean, path_sum, path_gtf, path_geom, path_coexp, path_fasta, path_interpro, path_ci)
        )
        thread.daemon = True
        thread.start()
        
        return False, {"display": "block"}, {"display": "none"}, True, True, True, True, True, True, True, True, True, True, True, True

    @app.callback(
        [Output("progress-interval", "disabled"),
         Output("progress-container", "style"),
         Output("progress-status-msg", "children"),
         Output("progress-bar-fill", "style"),
         Output("settings-feedback", "children"),
         Output("settings-feedback", "style"),
         Output("welcome-banner", "style"),
         Output("main-dashboard-content", "style"),
         Output("coexpression-network-container", "style"),
         Output("dataset-toggle", "options"),
         Output("data-sources-updated", "data"),
         Output("btn-apply-settings", "disabled", allow_duplicate=True),
         Output("btn-cancel-settings", "disabled", allow_duplicate=True),
         Output("btn-browse-mean", "disabled", allow_duplicate=True),
         Output("btn-browse-sum", "disabled", allow_duplicate=True),
         Output("btn-browse-gtf", "disabled", allow_duplicate=True),
         Output("btn-browse-geom", "disabled", allow_duplicate=True),
         Output("btn-browse-coexp", "disabled", allow_duplicate=True),
         Output("btn-browse-fasta", "disabled", allow_duplicate=True),
         Output("btn-browse-interpro", "disabled", allow_duplicate=True),
         Output("btn-browse-ci", "disabled", allow_duplicate=True),
         Output("btn-open-settings", "disabled", allow_duplicate=True),
         Output("btn-open-preprocess", "disabled", allow_duplicate=True)],
        [Input("progress-interval", "n_intervals")],
        [State("data-sources-updated", "data")],
        prevent_initial_call=True
    )
    def update_progress(n_intervals, current_updated):
        progress = state.get('loading_progress')
        if not progress:
            raise PreventUpdate
            
        if progress.get('error'):
            error_msg = progress['error']
            feedback_style = {"display": "block", "backgroundColor": "#fadbd8", "color": "#78281f"}
            return True, {"display": "none"}, "", {"width": "0%"}, error_msg, feedback_style, no_update, no_update, no_update, no_update, no_update, False, False, False, False, False, False, False, False, False, False, False, False
            
        if progress.get('done'):
            new_updated_time = (current_updated or 0) + 1
            has_sum = state['has_sum']
            gene_coexpression = state['gene_coexpression']
            coexp_style = {"display": "block" if gene_coexpression is not None else "none"}
            toggle_options = [
                {'label': ' Mean', 'value': 'mean'},
                {'label': ' Sum', 'value': 'sum', 'disabled': not has_sum}
            ]
            feedback_style = {"display": "block", "backgroundColor": "#d4efdf", "color": "#196f3d"}
            feedback_msg = "Data sources successfully loaded and applied!"
            
            return True, {"display": "none"}, "", {"width": "100%"}, feedback_msg, feedback_style, {"display": "none"}, {"display": "flex", "width": "100%", "marginBottom": "20px"}, coexp_style, toggle_options, new_updated_time, False, False, False, False, False, False, False, False, False, False, False, False
            
        step = progress.get('step', 0)
        total = progress.get('total', 8)
        msg = progress.get('msg', 'Loading...')
        percentage = int((step / total) * 100) if total > 0 else 0
        
        return False, {"display": "block"}, msg, {"width": f"{percentage}%", "height": "100%", "backgroundColor": "#2ecc71", "transition": "width 0.2s"}, no_update, no_update, no_update, no_update, no_update, no_update, no_update, True, True, True, True, True, True, True, True, True, True, True, True
    
    @app.callback(
        [Output("gene-table", "columns"),
         Output("gene-table", "data")],
        [Input("dataset-toggle", "value"),
         Input("cds-filter-toggle", "value"),
         Input("negative-spearman-toggle", "value"),
         Input("3d-filter-toggle", "value"),
         Input("domain-filter-toggle", "value"),
         Input("data-sources-updated", "data")],
    )
    def update_table(dataset, cds_filter, negative_spearman_filter, three_d_filter, domain_filter, updated_time):
        """Update table based on selected dataset, CDS filter, 3D filter, and domain filter."""
        table_data = state['table_mean_records'] if dataset == 'mean' else state['table_sum_records']

        if not table_data:
            return [], []

        # Get columns from first row
        columns = [{"name": col, "id": col, "editable": False} for col in table_data[0].keys()]
        
        # Set presentation for Gene Name column to render markdown links
        for col in columns:
            if col["id"] == "Gene Name":
                col["presentation"] = "markdown"

        # Apply CDS filter if enabled
        if cds_filter and 'filter_cds' in cds_filter:
            # Filter to only show genes with Has CDS = True
            table_data = [row for row in table_data if row.get("Has CDS") == True]

        # Apply negative Spearman filter if enabled
        if negative_spearman_filter and 'negative_only' in negative_spearman_filter:
            table_data = [row for row in table_data
                          if row.get("Min Spearman") is not None and row.get("Min Spearman") <= 0]

        # Apply 3D structure filter if enabled
        if three_d_filter and 'filter_3d' in three_d_filter:
            has_3d_col = any("Has 3D" in row for row in table_data[:1])
            if not has_3d_col:
                log.warning("update_table: '3D structure' filter triggered but 'Has 3D' column is absent "
                            "– was --geometry-dir passed at startup?")
            before = len(table_data)
            table_data = [row for row in table_data if row.get("Has 3D") is True]
            log.info("update_table: 3D filter reduced rows %d -> %d", before, len(table_data))

        # Apply domain filter if enabled
        if domain_filter and 'filter_domains' in domain_filter:
            has_domains_col = any("Has Domains" in row for row in table_data[:1])
            if not has_domains_col:
                log.warning("update_table: 'domains' filter triggered but 'Has Domains' column is absent "
                            "– check that InterPro results were loaded at startup")
            before = len(table_data)
            table_data = [row for row in table_data if row.get("Has Domains") is True]
            log.info("update_table: domain filter reduced rows %d -> %d", before, len(table_data))

        return columns, table_data

    @app.callback(
        [Output("cds-filter-toggle", "value"),
         Output("cds-filter-toggle", "options")],
        [Input("3d-filter-toggle", "value"),
         Input("domain-filter-toggle", "value")],
        [State("cds-filter-toggle", "value")],
    )
    def sync_cds_filter_with_structures(three_d_filter, domain_filter, current_cds):
        """Force CDS filter on (and disabled) when 3D or domain filters are active."""
        forced_on = bool((three_d_filter and 'filter_3d' in three_d_filter) or
                         (domain_filter and 'filter_domains' in domain_filter))
        options = [{'label': ' Genes with CDS', 'value': 'filter_cds', 'disabled': forced_on}]
        if forced_on:
            return ['filter_cds'], options
        return current_cds or [], options

    @app.callback(
        Output("selected-transcript", "data"),
        [Input("isoform-distribution", "clickData"),
         Input("exon-visualization", "clickData"),
         Input("coexpression-network", "tapNodeData"),
         Input("selected-gene", "data")],
        [State("selected-transcript", "data")]
    )
    def update_selected_transcript(dist_click, exon_click, tap_node, selected_gene, current_transcript):
        """Update selected transcript from distribution, exon visualization, or network clicks.

        Resets to None whenever the selected gene changes so the 3D viewer is
        cleared rather than showing a transcript from the previous gene.
        """
        if not callback_context.triggered:
            return current_transcript

        trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]

        # Gene changed → clear transcript selection unconditionally
        if trigger_id == "selected-gene":
            return None

        # Handle network click
        if trigger_id == "coexpression-network" and tap_node:
            if tap_node.get("type") == "isoform":
                transcript_id = tap_node["id"]
                if transcript_id == current_transcript:
                    return None
                return transcript_id

        # Handle distribution bar click
        if trigger_id == "isoform-distribution" and dist_click:
            if "points" in dist_click and dist_click["points"]:
                point = dist_click["points"][0]
                if "customdata" in point and point["customdata"]:
                    transcript_id = point["customdata"]
                    if transcript_id == current_transcript:
                        return None
                    return transcript_id
                elif "x" in point:
                    transcript_id = point["x"]
                    if transcript_id == current_transcript:
                        return None
                    return transcript_id
        
        # Handle exon visualization click — select/deselect isoform only
        if trigger_id == "exon-visualization" and exon_click:
            if "points" in exon_click and exon_click["points"]:
                point = exon_click["points"][0]
                if "customdata" in point and point["customdata"]:
                    cd = point["customdata"]
                    if isinstance(cd, dict) and cd.get("kind") == "domain":
                        return current_transcript
                    if isinstance(cd, (list, tuple)):
                        if cd and isinstance(cd[0], dict) and cd[0].get("kind") == "domain":
                            return current_transcript
                        if len(cd) > 1 and cd[1] == "domain":
                            return current_transcript
                    if isinstance(cd, str) and cd.startswith("IPR"):
                        return current_transcript
                    transcript_id = cd[0] if isinstance(cd, (list, tuple)) else cd
                    if transcript_id == current_transcript:
                        return None
                    return transcript_id
        
        return current_transcript


    @app.callback(
        Output("selected-gene", "data"),
        [Input("plot-entropy-scatter", "clickData"),
         Input("gene-table", "derived_virtual_selected_rows"),
         Input("coexpression-network", "tapNodeData")],
        [State("selected-gene", "data"),
         State("gene-table", "derived_virtual_data")]
    )
    def update_selected_gene(click_data, selected_rows, network_click, current_gene, table_data):
        """Update selected gene from scatter plot click, table selection, or network node click."""
        # Determine which input triggered the callback
        if not callback_context.triggered:
            return current_gene
        
        trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]
        
        # Handle network node click
        if trigger_id == "coexpression-network" and network_click:
            if network_click.get("type") == "gene":
                gene_id = network_click.get("id")
                if gene_id:
                    return gene_id
            else:
                raise PreventUpdate

        # Handle table row selection
        if trigger_id == "gene-table" and selected_rows and len(selected_rows) > 0:
            row_idx = selected_rows[0]
            if table_data and row_idx < len(table_data):
                gene_id = table_data[row_idx]["Gene ID"]
                return gene_id
        
        # Handle scatter plot click
        if trigger_id == "plot-entropy-scatter" and click_data:
            if "points" not in click_data or not click_data["points"]:
                return current_gene
            
            point = click_data["points"][0]
            gene_id = None
            cd_val = point.get("customdata")
            if cd_val is not None:
                if isinstance(cd_val, (list, tuple)) and cd_val:
                    gene_id = cd_val[0]
                elif isinstance(cd_val, str):
                    gene_id = cd_val
            if not gene_id:
                gene_id = point.get("hovertext") or point.get("text")
            if gene_id:
                return gene_id
        
        return current_gene

    # Cache the base scatter figure (no highlight) per dataset.
    # Rebuilding px.scatter over all genes is expensive; only the gold ring
    # overlay changes when a gene is selected.
    def _filtered_results(results_df, cds_filter, negative_spearman_filter, three_d_filter, domain_filter):
        filtered = results_df
        if filtered is None or filtered.empty:
            return pd.DataFrame()
        if cds_filter and 'filter_cds' in cds_filter and state['genes_with_cds']:
            filtered = filtered[filtered["gene_id"].isin(state['genes_with_cds'])]
        if negative_spearman_filter and 'negative_only' in negative_spearman_filter:
            filtered = filtered[filtered["min_spearman"] <= 0]
        if three_d_filter and 'filter_3d' in three_d_filter:
            if state['genes_with_3d']:
                filtered = filtered[filtered["gene_id"].isin(state['genes_with_3d'])]
            else:
                filtered = filtered.iloc[0:0]
        if domain_filter and 'filter_domains' in domain_filter:
            if state['genes_with_domains']:
                filtered = filtered[filtered["gene_id"].isin(state['genes_with_domains'])]
            else:
                filtered = filtered.iloc[0:0]
        return filtered

    @functools.lru_cache(maxsize=8)
    def _base_scatter(dataset: str, cds_only: bool, negative_only: bool,
                      three_d_only: bool, domains_only: bool):
        results_df = state['results_df_mean'] if dataset == 'mean' else state['results_df_sum']
        if results_df is None or results_df.empty:
            return go.Figure()
        cds_filter = ['filter_cds'] if cds_only else []
        negative_filter = ['negative_only'] if negative_only else []
        three_d_filter = ['filter_3d'] if three_d_only else []
        domain_filter = ['filter_domains'] if domains_only else []
        filtered = _filtered_results(results_df, cds_filter, negative_filter, three_d_filter, domain_filter)
        x_range, y_range = state['scatter_axis_ranges'].get(dataset, (None, None))
        return fig_summed_vs_top_entropy_colored_by_min_spearman(
            filtered,
            selected_gene=None,
            gene_names=state['gene_names'],
            x_range=x_range,
            y_range=y_range,
        )

    @app.callback(
        Output("plot-entropy-scatter", "figure"),
        [Input("selected-gene", "data"),
         Input("dataset-toggle", "value"),
         Input("cds-filter-toggle", "value"),
         Input("negative-spearman-toggle", "value"),
         Input("3d-filter-toggle", "value"),
         Input("domain-filter-toggle", "value"),
         Input("data-sources-updated", "data")],
    )
    def update_scatter_highlight(selected_gene, dataset, cds_filter,
                                 negative_spearman_filter, three_d_filter, domain_filter, updated_time):
        """Update scatter plot to highlight selected gene.

        The expensive base figure is cached per dataset; only the gold selection
        ring is added on top, so clicking a gene is near-instant.
        """
        import copy
        cds_only = bool(cds_filter and 'filter_cds' in cds_filter)
        negative_only = bool(negative_spearman_filter and 'negative_only' in negative_spearman_filter)
        three_d_only = bool(three_d_filter and 'filter_3d' in three_d_filter)
        domains_only = bool(domain_filter and 'filter_domains' in domain_filter)
        
        base_fig = _base_scatter(dataset, cds_only, negative_only, three_d_only, domains_only)
        if base_fig is None:
            return go.Figure()
            
        fig = copy.deepcopy(base_fig)
        x_range, y_range = state['scatter_axis_ranges'].get(dataset, (None, None))
        if x_range is not None:
            fig.update_xaxes(range=x_range, autorange=False)
        if y_range is not None:
            fig.update_yaxes(range=y_range, autorange=False)
        if not selected_gene:
            return fig

        results_df = state['results_df_mean'] if dataset == 'mean' else state['results_df_sum']
        if results_df is None or results_df.empty:
            return fig
            
        results_df = _filtered_results(results_df, cds_filter, negative_spearman_filter,
                                       three_d_filter, domain_filter)
        if not results_df.empty:
            sel = results_df[results_df["gene_id"] == selected_gene]
            if not sel.empty:
                fig.add_trace(go.Scatter(
                    x=sel["summed_isoform_entropy"],
                    y=sel["top_isoform_entropy"],
                    mode='markers',
                    marker=dict(size=20, color='rgba(255,215,0,0)',
                                line=dict(color='gold', width=3)),
                    showlegend=False,
                    hoverinfo='skip',
                    name='Selected',
                ))
        return fig

    @app.callback(
        Output("gene-table", "selected_rows"),
        Input("selected-gene", "data"),
        State("gene-table", "derived_virtual_data")
    )
    def update_table_selection(selected_gene, table_data):
        """Highlight the selected gene in the table."""
        if not selected_gene or not table_data:
            return []
        
        for idx, row in enumerate(table_data):
            if row["Gene ID"] == selected_gene:
                return [idx]
        return []

    # Cache the exon structure figure per (gene_id, dataset).
    # The exon bars never change when only the selected transcript changes;
    # only the highlight overlay does.  We cache the un-highlighted version and
    # reapply the highlight cheaply.
    @functools.lru_cache(maxsize=256)
    def _base_exon_fig(gene_id: str, dataset: str):
        df = state['df_mean'] if dataset == 'mean' else state['df_sum']
        if df is None or df.empty:
            return go.Figure()
        global_col = state['global_col_mean'] if dataset == 'mean' else state['global_col_sum']
        sub = df[df["gene_id"] == gene_id]
        if not sub.empty:
            max_display = 120
            if len(sub) > max_display:
                sub = sub.nlargest(max_display, global_col)
            sub = sub.sort_values(global_col, ascending=False)
            transcript_order = sub["transcript_id"].tolist()
        else:
            transcript_order = None
        return create_exon_visualization(
            gene_id,
            state['isoforms_by_gene'],
            transcript_order=transcript_order,
            selected_transcript=None,      # no highlight in cached version
            af_geometry_mapping=state['af_geometry_mapping'],
            domain_data=state['domain_mapping'],
        )

    @app.callback(
        [Output("exon-visualization", "figure"),
         Output("exon-structure-label", "children")],
        [Input("selected-gene", "data"),
         Input("selected-transcript", "data"),
        Input("selected-domain", "data"),
         Input("dataset-toggle", "value")],
    )
    def render_exon_visualization(gene_id, selected_transcript, selected_domain, dataset):
        """Render exon structure diagram.

        The base figure (no highlight) is cached per gene/dataset so that only
        switching genes triggers a full redraw; selecting a transcript re-renders
        with the highlight but the expensive layout work is already done.
        """
        if not gene_id:
            return go.Figure().update_layout(title=None), "No Exon Structures to show (select a gene from the table)"
        import copy

        def _exon_label_for_gene(target_gene_id: str) -> str:
            isoforms = state['isoforms_by_gene'].get(target_gene_id, {})
            gene_name = state['gene_names'].get(target_gene_id, "")
            gene_display = gene_name if gene_name else target_gene_id
            if not isoforms:
                return f"Exon structure of gene: {gene_display}"
            gene_start = min(exon["exon_start"] for tx in isoforms.values() for exon in tx)
            gene_end = max(exon["exon_end"] for tx in isoforms.values() for exon in tx)
            gene_length = gene_end - gene_start + 1
            has_domains = any(state['domain_mapping'].get(tid) for tid in isoforms.keys())
            domain_suffix = ", domains shown" if has_domains else ""
            return (
                f"Exon structure of gene: {gene_display} "
                f"({gene_length:,} bp, compressed introns{domain_suffix})"
            )

        label = _exon_label_for_gene(gene_id)
        fig = copy.deepcopy(_base_exon_fig(gene_id, dataset))

        if not selected_transcript:
            if selected_domain:
                _highlight_selected_domain(fig, selected_domain)
            return fig, label

        meta = fig.layout.meta or {}
        transcript_order = meta.get("transcript_order") or []
        if selected_transcript not in transcript_order:
            return fig, label

        row_height = meta.get("row_height", Dimensions.ROW_HEIGHT)
        total_width = meta.get("total_width")
        if total_width is None:
            # Fallback: compute max x across traces
            total_width = 0.0
            for trace in fig.data:
                xs = getattr(trace, "x", None)
                if xs is not None:
                    total_width = max(total_width, max(xs))

        idx = transcript_order.index(selected_transcript)
        n_isoforms = len(transcript_order)
        y_pos = (n_isoforms - idx - 1) * row_height

        # Background highlight band
        fig.add_trace(go.Scatter(
            x=[0, total_width, total_width, 0, 0],
            y=[y_pos, y_pos, y_pos + row_height, y_pos + row_height, y_pos],
            fill='toself',
            fillcolor='rgba(0, 191, 255, 0.15)',
            line=dict(color='rgba(0, 191, 255, 0)', width=0),
            mode='lines',
            showlegend=False,
            hoverinfo='skip'
        ))
        

        if selected_domain:
            _highlight_selected_domain(fig, selected_domain)

        return fig, label

    def _highlight_selected_domain(fig, selected_domain):
        """Apply a subtle outline to the selected domain rectangle."""
        if not selected_domain:
            return

        for trace in fig.data:
            customdata = getattr(trace, "customdata", None)
            if not customdata:
                continue
            cd0 = customdata[0] if isinstance(customdata, (list, tuple)) else customdata
            if not isinstance(cd0, dict) or cd0.get("kind") != "domain":
                continue
            if (
                cd0.get("transcript_id") == selected_domain.get("transcript_id")
                and cd0.get("genomic_start") == selected_domain.get("genomic_start")
                and cd0.get("genomic_end") == selected_domain.get("genomic_end")
                and cd0.get("name") == selected_domain.get("name")
            ):
                trace.line.color = "#39FF14"
                trace.line.width = 2.5
                trace.opacity = 0.98

    @app.callback(
        [Output("domain-details-text", "children"),
         Output("domain-details-link", "children"),
        Output("domain-details-wrapper", "style"),
         Output("selected-domain", "data")],
        [Input("exon-visualization", "clickData"),
         Input("selected-gene", "data")],
    )
    def update_domain_details(exon_click, selected_gene):
        """Show copyable domain details when a domain rectangle is clicked."""
        hidden_style = {"display": "none"}
        if callback_context.triggered:
            trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]
            if trigger_id == "selected-gene":
                return "", "", hidden_style, None

        if not exon_click or "points" not in exon_click or not exon_click["points"]:
            return "", "", hidden_style, None

        point = exon_click["points"][0]
        cd = point.get("customdata")
        domain_payload = None

        if isinstance(cd, dict) and cd.get("kind") == "domain":
            domain_payload = cd
        elif isinstance(cd, (list, tuple)) and cd:
            if isinstance(cd[0], dict) and cd[0].get("kind") == "domain":
                domain_payload = cd[0]

        if not domain_payload:
            return "", "", hidden_style, None

        details_lines = [
            f"Name: {domain_payload.get('name', 'N/A')}",
            f"Accession: {domain_payload.get('accession', 'N/A')}",
            f"Library: {domain_payload.get('library', 'N/A')}",
            f"AA: {domain_payload.get('aa_start', 'N/A')} - {domain_payload.get('aa_end', 'N/A')}",
            f"Genomic: {domain_payload.get('genomic_start', 'N/A')} - {domain_payload.get('genomic_end', 'N/A')}",
        ]
        description = domain_payload.get("description")
        if description:
            details_lines.insert(1, f"Description: {description}")

        details_text = "\n".join(details_lines)

        interpro_id = domain_payload.get("interpro_id")
        if interpro_id:
            link_url = f"https://www.ebi.ac.uk/interpro/entry/InterPro/{interpro_id}/"
            link_child = html.A(
                interpro_id,
                href=link_url,
                target="_blank",
                rel="noopener noreferrer",
            )
            details_lines.insert(0, f"InterPro: {interpro_id}")
        else:
            link_child = html.Span("InterPro entry not available")

        visible_style = {
            "flex": "1",
            "padding": "8px",
            "border": "1px solid #ddd",
            "borderRadius": "4px",
            "backgroundColor": "#fafafa",
            "display": "block",
        }
        return details_text, link_child, visible_style, domain_payload

    @app.callback(
        [Output("transcript-details-text", "children"),
         Output("transcript-details-wrapper", "style")],
        [Input("selected-transcript", "data"),
         Input("selected-gene", "data"),
         Input("dataset-toggle", "value")],
    )
    def update_transcript_details(selected_transcript, selected_gene, dataset):
        """Show transcript details when a transcript is selected."""
        hidden_style = {"display": "none"}
        if not selected_transcript or not selected_gene:
            return "", hidden_style

        df = state['df_mean'] if dataset == 'mean' else state['df_sum']
        if df is None or df.empty:
            return "", hidden_style
        global_col = state['global_col_mean'] if dataset == 'mean' else state['global_col_sum']
        sub = df[df["gene_id"] == selected_gene]
        transcript_order = []
        if not sub.empty:
            max_display = 120
            if len(sub) > max_display:
                sub = sub.nlargest(max_display, global_col)
            sub = sub.sort_values(global_col, ascending=False)
            transcript_order = sub["transcript_id"].tolist()

        index_text = "N/A"
        rank_text = "N/A"
        tpm_text = "N/A"
        if selected_transcript in transcript_order:
            index_text = f"{transcript_order.index(selected_transcript) + 1} / {len(transcript_order)}"
            tpm_val = sub.loc[sub["transcript_id"] == selected_transcript, global_col]
            if not tpm_val.empty:
                tpm_text = f"{float(tpm_val.iloc[0]):.2f}"

        gene_name = state['gene_names'].get(selected_gene, "")
        gene_display = f"{selected_gene} ({gene_name})" if gene_name else selected_gene
        details_lines = [
            f"Id: {selected_transcript}",
            f"Gene: {gene_display}",
            f"Index: {index_text}",
            f"TPM: {tpm_text}",
        ]
        details_text = "\n".join(details_lines)

        visible_style = {
            "flex": "1",
            "padding": "8px",
            "border": "1px solid #ddd",
            "borderRadius": "4px",
            "backgroundColor": "#fafafa",
            "display": "block",
        }
        return details_text, visible_style

    @app.callback(
        Output("selection-details-wrapper", "style"),
        [Input("selected-transcript", "data"),
         Input("selected-domain", "data")],
    )
    def toggle_selection_details_row(selected_transcript, selected_domain):
        """Show the details row when either transcript or domain is selected."""
        if not selected_transcript and not selected_domain:
            return {"display": "none"}
        return {
            "marginTop": "8px",
            "display": "flex",
            "gap": "10px",
        }

    @app.callback(
       [Output("isoform-distribution", "figure"),
        Output("selected-gene-label", "children"),
        Output("ci-toggle", "style"),
        Output("ci-toggle", "value")],
       [Input("selected-gene", "data"),
        Input("selected-transcript", "data"),
        Input("ci-toggle", "value"),
        Input("dataset-toggle", "value")],
    )
    def render_isoform_distribution(gene_id, selected_transcript,
                             ci_toggle_value, dataset):
        """Render isoform distribution based on selected gene.
        
        In "mean" mode, the CI checkbox controls whether CIs are displayed.
        In "sum" mode, CIs are never shown and the checkbox is unchecked
        and visually disabled.
        """
        if not gene_id:
            empty_fig = go.Figure().update_layout(title="Click a point in the scatter plot or select from table to see isoform distributions")
            # Hide CI toggle entirely when no gene is selected
            return empty_fig, "No Gene distributions to show (select a gene from the table)", {'display': 'none'}, []

        # Use server-side data instead of browser storage
        df = state['df_mean'] if dataset == 'mean' else state['df_sum']
        if df is None or df.empty:
            return go.Figure().update_layout(title="No dataset loaded"), "No dataset loaded", {'display': 'none'}, []
            
        has_ci = bool(state['ci_dict'] and state['ci_columns'])  # CIs available if we have the data
        global_col = state['global_col_mean'] if dataset == 'mean' else state['global_col_sum']

        # Determine if CIs should be shown based on toggle.
        # In sum mode, CIs are always disabled.
        if dataset == 'sum':
            show_ci = False
        else:
            show_ci = has_ci and ('show_ci' in ci_toggle_value)

        # Determine CI mode from available columns and sample columns
        ci_mode = None
        if show_ci and state['ci_columns']:
            # CI columns are named: ci_lower_{group_name}, ci_upper_{group_name}
            # Extract all group names from CI columns (e.g., "Bipolar" from "ci_lower_Bipolar")
            ci_group_names = set()
            for col in state['ci_columns']:
                # Format: ci_lower_{name} or ci_upper_{name}
                parts = col.split('_', 2)  # Split only on first 2 underscores
                if len(parts) == 3:
                    ci_group_names.add(parts[2])
            
            # Match CI group names against sample columns
            # If sample columns are condition names, they'll match the condition CI groups
            # If sample columns are region names, they'll match the region CI groups
            matching_groups = ci_group_names & set(state['sample_cols'])
            
            if matching_groups:
                # Sample columns directly match a grouping in the CIs
                # Default to condition mode (more common)
                ci_mode = 'condition'
            else:
                # Fall back based on whether we have regional or condition-specific CIs
                # This is a fallback when sample columns don't directly match
                ci_mode = 'global'

        # Create distribution figure
        fig_dist = fig_isoform_sample_panels(
            gene_id, df, state['sample_cols'], show_ci, global_col,
            selected_transcript, state['ci_dict'], ci_mode, state['sample_cols']
        )

        gene_name = state['gene_names'].get(gene_id, "")
        gene_display = gene_name if gene_name else gene_id
        label = f"Isoform Expression  ·  {gene_display}"
        if selected_transcript:
            label += f"  (selected: {selected_transcript})"

        # Configure CI toggle style and value.
        # - Hidden entirely if no CI data.
        # - Disabled and unchecked in sum mode.
        # - Enabled in mean mode.
        if not has_ci:
            ci_toggle_style = {'display': 'none'}
            ci_toggle_value_out = []
        elif dataset == 'sum':
            ci_toggle_style = {
                'display': 'inline-block',
                'pointerEvents': 'none',
                'opacity': 0.4,
            }
            ci_toggle_value_out = []
        else:
            ci_toggle_style = {'display': 'inline-block'}
            ci_toggle_value_out = ci_toggle_value

        return fig_dist, label, ci_toggle_style, ci_toggle_value_out

    # ── Exon selector bar ─────────────────────────────────────────────────

    @app.callback(
        [Output("exon-selector-bar", "children"),
         Output("exon-selector-bar", "style")],
        [Input("selected-transcript", "data"),
         Input("selected-exon", "data")],
    )
    def populate_exon_selector(selected_transcript, selected_exon):
        """Build the row of exon-highlight buttons for the selected transcript.

        Buttons are only shown when the transcript has pre-generated exon viewers.
        The currently-active exon (or 'Full') is styled differently.
        """
        hidden = {"display": "none"}
        if not selected_transcript:
            return [], hidden

        geo = resolve_alphafold_geometry(selected_transcript, state['af_geometry_mapping'])
        if not geo:
            return [], hidden

        exon_htmls = geo.get("exon_htmls", [])
        if not exon_htmls or not any(exon_htmls):
            return [], hidden

        # Which exon is currently highlighted? None → "Full structure" is active
        active_exon_idx = None
        if (selected_exon is not None and
                selected_exon.get("transcript_id") == selected_transcript):
            active_exon_idx = selected_exon.get("exon_idx")

        _btn_base = {
            "padding": "4px 10px",
            "border": "1px solid #aaa",
            "borderRadius": "4px",
            "cursor": "pointer",
            "fontSize": "12px",
            "marginRight": "4px",
            "marginBottom": "4px",
        }
        _active_extra   = {"backgroundColor": "#00BFFF", "color": "white",
                            "fontWeight": "bold", "borderColor": "#0090c0"}
        _inactive_extra = {"backgroundColor": "#f5f5f5", "color": "#333"}

        def _btn_style(is_active):
            return {**_btn_base, **(_active_extra if is_active else _inactive_extra)}

        buttons = [
            html.Button(
                "Full structure",
                id={"type": "exon-btn", "index": -1},
                n_clicks=0,
                style=_btn_style(active_exon_idx is None),
            )
        ]
        # Only create buttons for CDS exons (non-None entries).
        # Label by overall (genomic) exon number — same index used everywhere.
        for genomic_idx, html_path in enumerate(exon_htmls):
            if html_path is None:
                continue
            buttons.append(
                html.Button(
                    f"Exon {genomic_idx + 1}",
                    id={"type": "exon-btn", "index": genomic_idx},
                    n_clicks=0,
                    style=_btn_style(active_exon_idx == genomic_idx),
                )
            )

        bar_style = {
            "display": "flex",
            "flexWrap": "wrap",
            "alignItems": "center",
            "gap": "0",
            "padding": "6px 0 4px 0",
            "borderBottom": "1px solid #ddd",
            "marginBottom": "8px",
        }
        label = html.Span(
            "Highlight exon: ",
            style={"fontSize": "12px", "fontWeight": "bold",
                   "marginRight": "6px", "whiteSpace": "nowrap"},
        )
        return [label, *buttons], bar_style

    @app.callback(
        Output("selected-exon", "data"),
        [Input({"type": "exon-btn", "index": ALL}, "n_clicks"),
         Input("selected-transcript", "data")],
        [State({"type": "exon-btn", "index": ALL}, "id"),
         State("selected-transcript", "data"),
         State("selected-exon", "data")],
        prevent_initial_call=True,
    )
    def exon_button_click(all_n_clicks, selected_transcript_input,
                          all_ids, selected_transcript, current_exon):
        """Update selected-exon when an exon button is clicked or transcript changes."""
        if not callback_context.triggered:
            raise PreventUpdate

        trigger_id = callback_context.triggered[0]["prop_id"]

        # Transcript changed → clear any exon selection
        if "selected-transcript" in trigger_id:
            return None

        # Otherwise a button was clicked
        if not selected_transcript:
            raise PreventUpdate

        import json as _json
        try:
            btn_id = _json.loads(trigger_id.split(".")[0])
        except Exception:
            raise PreventUpdate

        exon_idx = btn_id.get("index")
        if exon_idx is None:
            raise PreventUpdate

        # "Full structure" button (index -1) → clear exon selection
        if exon_idx == -1:
            return None

        # Same button clicked again → deselect (back to full)
        if (current_exon is not None and
                current_exon.get("transcript_id") == selected_transcript and
                current_exon.get("exon_idx") == exon_idx):
            return None

        return {"transcript_id": selected_transcript, "exon_idx": exon_idx}

    # ── In-memory cache for patched AlphaFold HTML files ─────────────────
    # Reading + regex-patching the HTML on every exon button click is slow.
    # Cache the result keyed by absolute file path.
    _html_content_cache: dict = {}

    import re as _re

    _STYLE_PATCH = (
        "<style>"
        "html,body{margin:0;padding:0;background:#ffffff!important;"
        "overflow:hidden;width:100%;height:100%;}"
        "body{display:flex;flex-direction:column;"
        "align-items:center;justify-content:center;}"
        "div[id^='3dmolviewer']{"
        "width:100%!important;max-width:100%!important;"
        "margin:0 auto!important;}"
        "canvas{display:block;margin:auto;}"
        "</style>"
    )
    _JS_PATCH = (
        "<script>"
        "window.addEventListener('load',function(){"
        "var els=document.querySelectorAll(\"div[id^='3dmolviewer']\");"
        "els.forEach(function(el){"
        "el.style.width='100%';"
        "el.style.maxWidth='100%';"
        "el.style.margin='0 auto';"
        "});"
        "});"
        "</script>"
    )

    def _load_patched_html(path: str) -> str:
        """Return patched HTML content, reading from disk only on first access."""
        if path not in _html_content_cache:
            content = open(path, "r").read()
            # Remove dark legend baked in by extract_3d_geometry.py
            content = _re.sub(
                r'<div[^>]*background:#1a1a2e[^>]*>.*?pLDDT colour key.*?</div>',
                '',
                content,
                flags=_re.DOTALL,
            )
            content = content.replace("</head>", _STYLE_PATCH + _JS_PATCH + "</head>", 1)
            _html_content_cache[path] = content
        return _html_content_cache[path]

    @app.callback(
        Output("glb-viewer-container", "children"),
        [Input("selected-transcript", "data"),
         Input("selected-gene", "data"),
         Input("selected-exon", "data")],
    )
    def update_glb_viewer(selected_transcript, selected_gene, selected_exon):
        """Update 3D structure viewer with AlphaFold HTML viewer.

        When an exon is selected (selected_exon store is set and matches the
        current transcript), load the pre-generated exon-highlighted HTML
        instead of the full-protein viewer.
        """
        if not selected_transcript:
            return html.P("Select a transcript with a 3D model to view it here",
                          style={"textAlign": "center", "color": "#999", "padding": "20px"})

        if not state['af_geometry_mapping'] or not selected_gene:
            return html.P(
                f"No 3D model available for {selected_transcript}",
                style={"textAlign": "center", "color": "#999", "padding": "20px"},
            )

        geo = resolve_alphafold_geometry(selected_transcript, state['af_geometry_mapping'])

        if not geo:
            return html.P(
                f"No 3D model available for {selected_transcript}",
                style={"textAlign": "center", "color": "#999", "padding": "20px"},
            )

        # ── Determine which HTML file to serve ───────────────────────────
        # Default: full-protein viewer
        html_path = geo.get("html")
        # If this entry is a sequence-identity alias, note the original transcript
        if geo.get("is_sequence_alias"):
            alias_of = geo.get("alias_of", "")
            viewer_label = (
                f"AlphaFold 3D – {selected_transcript} "
                f"(identical sequence to {alias_of})"
            )
        else:
            viewer_label = f"AlphaFold 3D – {selected_transcript}"

        # Override with exon-specific HTML if an exon is selected for this transcript
        if (selected_exon is not None and
                selected_exon.get("transcript_id") == selected_transcript):
            exon_idx = selected_exon.get("exon_idx")
            exon_htmls = geo.get("exon_htmls", [])
            if (exon_idx is not None and
                    exon_idx < len(exon_htmls) and
                    exon_htmls[exon_idx] and
                    os.path.isfile(exon_htmls[exon_idx])):
                html_path = exon_htmls[exon_idx]
                viewer_label = f"{viewer_label} · Exon {exon_idx + 1}"

        # ── AlphaFold HTML viewer ─────────────────────────────────────────
        if html_path and os.path.isfile(html_path):
            html_content = _load_patched_html(html_path)

            children = [
                html.P(
                    viewer_label,
                    style={"fontWeight": "bold", "marginBottom": "4px", "fontSize": "13px"},
                ),
            ]
            children += [
                html.Iframe(
                    srcDoc=html_content,
                    style={"width": "100%", "height": "460px",
                           "border": "none", "display": "block"},
                ),
                # pLDDT colour legend
                html.Div([
                    html.Span("pLDDT: ", style={"fontWeight": "bold", "marginRight": "6px"}),
                    html.Span("■", style={"color": "royalblue",       "marginRight": "3px"}),
                    html.Span("≥90 Very high",   style={"marginRight": "10px"}),
                    html.Span("■", style={"color": "cornflowerblue", "marginRight": "3px"}),
                    html.Span("70–89 Confident", style={"marginRight": "10px"}),
                    html.Span("■", style={"color": "#bbbb00",        "marginRight": "3px"}),
                    html.Span("50–69 Low",        style={"marginRight": "10px"}),
                    html.Span("■", style={"color": "orange",         "marginRight": "3px"}),
                    html.Span("<50 Very low"),
                ], style={"fontSize": "11px", "marginTop": "4px",
                          "fontFamily": "Arial, sans-serif"}),
            ]
            return html.Div(children)

        # HTML not present (--no-viewer was used) — show CSV summary instead
        if geo.get("csv") and os.path.isfile(geo["csv"]):
            try:
                df_geo = pd.read_csv(geo["csv"])
                n_res = len(df_geo)
                mean_plddt = (df_geo["mean_plddt"].mean()
                              if "mean_plddt" in df_geo.columns else None)
                pct_conf = (
                    (df_geo["mean_plddt"] >= 70).mean() * 100
                    if "mean_plddt" in df_geo.columns else None
                )
                return html.Div([
                    html.P(f"AlphaFold geometry: {selected_transcript}",
                           style={"fontWeight": "bold", "textAlign": "center"}),
                    html.P(
                        (f"Residues: {n_res}  |  Mean pLDDT: {mean_plddt:.1f}  |  "
                         f"Confident (≥70): {pct_conf:.0f}%")
                        if mean_plddt else f"Residues: {n_res}",
                        style={"textAlign": "center", "color": "#555"},
                    ),
                    html.P(
                        "Interactive viewer not available. "
                        "Re-run extract_3d_geometry.py without --no-viewer.",
                        style={"textAlign": "center", "color": "#999", "fontSize": "12px"},
                    ),
                ])
            except Exception:
                pass

        return html.P(
            f"No 3D model available for {selected_transcript}",
            style={"textAlign": "center", "color": "#999", "padding": "20px"},
        )

    @app.callback(
        Output("protein-sequence-container", "children"),
        [Input("selected-transcript", "data"),
         Input("selected-exon", "data")],
    )
    def update_protein_sequence(selected_transcript, selected_exon):
        """Update protein sequence display, highlighting the selected exon's residues."""
        if not selected_transcript or not state['protein_sequences']:
            return [
                html.Div([
                    html.Div([
                        dcc.Clipboard(
                            id="protein-sequence-clipboard",
                            target_id="protein-sequence-text",
                            style={"float": "right"},
                            title="Copy protein sequence",
                        ),
                    ], style={"marginBottom": "6px"}),
                    html.Pre(
                        id="protein-sequence-text",
                        style={
                            "whiteSpace": "pre-wrap",
                            "margin": 0,
                            "fontFamily": "monospace",
                            "fontSize": "12px",
                            "wordBreak": "break-all",
                            "backgroundColor": "#f9f9f9",
                            "border": "none",
                            "padding": 0,
                        },
                        children="Select a transcript to view its protein sequence"
                    )
                ])
            ]

        if selected_transcript not in state['protein_sequences']:
            return [
                html.Div([
                    html.Div([
                        dcc.Clipboard(
                            id="protein-sequence-clipboard",
                            target_id="protein-sequence-text",
                            style={"float": "right"},
                            title="Copy protein sequence",
                        ),
                    ], style={"marginBottom": "6px"}),
                    html.Pre(
                        id="protein-sequence-text",
                        style={
                            "whiteSpace": "pre-wrap",
                            "margin": 0,
                            "fontFamily": "monospace",
                            "fontSize": "12px",
                            "wordBreak": "break-all",
                            "backgroundColor": "#f9f9f9",
                            "border": "none",
                            "padding": 0,
                        },
                        children=f"No protein sequence available for {selected_transcript}"
                    )
                ])
            ]

        sequence = state['protein_sequences'][selected_transcript]

        # --- Determine highlight range (1-based, inclusive) ---
        hi_start = hi_end = None  # None → full structure (no highlight)
        exon_label = None

        if (selected_exon is not None and
                selected_exon.get("transcript_id") == selected_transcript):
            exon_idx = selected_exon.get("exon_idx")
            if exon_idx is not None:
                exons_sorted = state['transcript_to_exons'].get(selected_transcript)
                if exons_sorted is not None and exon_idx < len(exons_sorted):
                    res_range = compute_exon_residue_range(exon_idx, exons_sorted)
                    if res_range is not None:
                        hi_start, hi_end = res_range
                        exon_label = f"Exon {exon_idx + 1}"

        # --- Build colour-coded sequence spans (60 aa per row) ---
        chunk_size = 60
        rows = []
        for chunk_start in range(0, len(sequence), chunk_size):
            chunk = sequence[chunk_start:chunk_start + chunk_size]
            row_spans = []
            for char_i, aa in enumerate(chunk):
                res_1based = chunk_start + char_i + 1  # 1-based residue number
                if hi_start is not None and hi_start <= res_1based <= hi_end:
                    row_spans.append(html.Span(
                        aa,
                        style={
                            "backgroundColor": "#FFD700",
                            "color": "#000",
                            "fontWeight": "bold",
                            "borderRadius": "2px",
                        }
                    ))
                else:
                    color = "#aaa" if hi_start is not None else "#333"
                    row_spans.append(html.Span(aa, style={"color": color}))
            rows.append(row_spans)

        # Build one sequence element for both display and clipboard target.
        highlighted_children = []
        for i, row_spans in enumerate(rows):
            highlighted_children.extend(row_spans)
            if i < len(rows) - 1:
                highlighted_children.append("\n")

        # --- Header line ---
        if hi_start is not None:
            exon_aa = hi_end - hi_start + 1
            subtitle = (f"{exon_label}  ·  {exon_aa} aa"
                        f"  ·  residues {hi_start}–{hi_end}"
                        f"  of  {len(sequence)} aa")
            subtitle_color = "#555"
        else:
            subtitle = f"{len(sequence)} aa"
            subtitle_color = "#333"

        return [
            html.Div([
                html.Div([
                    dcc.Clipboard(
                        id="protein-sequence-clipboard",
                        target_id="protein-sequence-text",
                        style={"float": "right"},
                        title="Copy protein sequence",
                    ),
                ], style={"marginBottom": "6px"}),
                html.P([
                    html.Span(selected_transcript,
                              style={"fontWeight": "bold", "fontFamily": "Arial, sans-serif"}),
                    html.Span(f"  ·  {subtitle}",
                              style={"fontWeight": "bold", "fontFamily": "Arial, sans-serif",
                                     "marginLeft": "6px", "color": subtitle_color}),
                ], style={"marginBottom": "8px"}),
                html.Pre(
                    id="protein-sequence-text",
                    style={
                        "whiteSpace": "pre-wrap",
                        "margin": 0,
                        "fontFamily": "monospace",
                        "fontSize": "12px",
                        "wordBreak": "break-all",
                        "letterSpacing": "0.5px",
                        "lineHeight": "1.6",
                        "backgroundColor": "#f9f9f9",
                        "border": "none",
                        "padding": 0,
                    },
                    children=highlighted_children,
                ),
            ])
        ]

    @app.callback(
        [Output("expanded-network-genes", "data"),
         Output("expanded-gene-position", "data")],
        [Input("coexpression-network", "tapNode"),
         Input("selected-gene", "data")],
        State("expanded-network-genes", "data"),
    )
    def handle_network_click(tap_node, selected_gene, expanded_genes):
        """Manage gene-node expansion in the network.

        Two triggers:
        - selected-gene changes: clear all expansions so the new ego network
          starts clean with no stale compound nodes.
        - tapNode (user clicks a node): enforce single-expansion and capture
          the node's current rendered position so the compound box can be
          offset from that position when re-drawn.
        """
        if not callback_context.triggered:
            raise PreventUpdate

        trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]

        # Gene selection changed
        if trigger_id == "selected-gene":
            # If the newly selected gene is already in the expanded genes list
            # (which happens if the user clicked the node in the network itself),
            # do not wipe the expansions. Otherwise, wipe expansions for a clean start.
            if expanded_genes and selected_gene in expanded_genes:
                return no_update, no_update
            return [], None

        # Node tapped in the network
        if not tap_node:
            raise PreventUpdate

        node_data = tap_node.get('data', {})
        position  = tap_node.get('position', None)  # {'x': ..., 'y': ...}

        if node_data.get('type') != 'gene':
            raise PreventUpdate

        gene_id = node_data['id']
        expanded_genes = expanded_genes or []

        if gene_id in expanded_genes:
            # Clicking the already-open gene collapses it
            return [], None
        else:
            # Replace whatever was open with just this gene; record its position
            return [gene_id], position

    @app.callback(
        [Output("coexpression-network", "elements"),
         Output("coexpression-network", "layout")],
        [Input("selected-gene", "data"),
         Input("expanded-network-genes", "data"),
         Input("dataset-toggle", "value"),
         Input("network-threshold-input", "value"),
         Input("network-neighbors-input", "value")],
        State("expanded-gene-position", "data"),
    )
    def update_network(selected_gene, expanded_genes, dataset, threshold, max_neighbors, expanded_gene_position):
        """Update network elements based on selected gene, expanded nodes, threshold, and max neighbors."""
        cyto_layout = {
            'name': 'cose',
            'idealEdgeLength': 60,
            'nodeOverlap': 10,
            'refresh': 20,
            'fit': True,
            'padding': 20,
            'randomize': False,
            'componentSpacing': 60,
            'nodeRepulsion': 3000,
            'edgeElasticity': 80,
            'nestingFactor': 0.1,
            'gravity': 100,
            'numIter': 1000,
            'initialTemp': 200,
            'coolingFactor': 0.95,
            'minTemp': 1.0,
            'animate': True,
        }

        if state['gene_coexpression'] is None:
            return [], cyto_layout

        # Pass the active expression DataFrame so the isoform list is filtered
        # by expression level — matching the exon-panel's nlargest(120, global_col) cap.
        active_df  = state['df_mean'] if dataset == 'mean' else state['df_sum']
        if active_df is None or active_df.empty:
            return [], cyto_layout
            
        active_col = state['global_col_mean'] if dataset == 'mean' else state['global_col_sum']

        # Handle None threshold and neighbors values from layout gracefully
        threshold_val = float(threshold) if threshold is not None else 0.3
        max_neighbors_val = int(max_neighbors) if max_neighbors is not None else 10

        elements = generate_network_elements(
            gene_coexpression=state['gene_coexpression'],
            gene_coexpression_idx=state['gene_coexpression_idx'],
            isoform_coexpression=state['isoform_coexpression'],
            isoform_coexpression_idx=state['isoform_coexpression_idx'],
            target_gene=selected_gene,
            expanded_genes=expanded_genes,
            top_k_genes=max_neighbors_val,
            top_k_isoforms=10,
            isoforms_by_gene=state['isoforms_by_gene'],
            gene_names=state['gene_names'],
            df_expression=active_df,
            global_col=active_col,
            max_isoforms=120,
            expanded_gene_position=expanded_gene_position,
            threshold=threshold_val,
        )
        return elements, cyto_layout

    @app.callback(
        [Output("coexpression-network", "stylesheet"),
         Output("coexpression-network", "layout")],
        Input("selected-transcript", "data"),
    )
    def update_network_stylesheet(selected_transcript):
        """Dynamically update coexpression network stylesheet to highlight the selected transcript
        without triggering a full physics layout re-run.
        """
        from .coexpression_network import get_cytoscape_stylesheet
        base_style = get_cytoscape_stylesheet()
        
        # When styling/highlighting changes, lock positions to 'preset' to prevent
        # Cytoscape from re-running the physics (cose) layout and shifting elements.
        preset_layout = {'name': 'preset'}

        if not selected_transcript:
            return base_style, no_update

        # Append a rule specifically for the selected transcript ID
        highlight_rule = {
            'selector': f'node[id = "{selected_transcript}"]',
            'style': {
                'background-color': '#6D28D9',
                'background-opacity': 1.0,
                'border-width': 1.5,
                'border-color': '#ffffff',
                'width': '18px',
                'height': '18px',
                'font-size': '8px',
                'color': '#ffffff',
                'text-outline-color': '#4C1D95',
                'text-outline-width': 1,
                'z-index': 10,
            }
        }
        return base_style + [highlight_rule], preset_layout