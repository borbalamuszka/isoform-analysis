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

log = logging.getLogger(__name__)

from .data_processing import (
    prepare_table_data,
    compute_min_spearman_per_gene,
    compute_gene_ranking,
    compute_gene_ranking_by_expression,
    gene_has_cds,
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
from .config import Colors, Dimensions

# Residue-range helper (pure function, no heavy deps)
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utilities.extract_3d_geometry import compute_exon_residue_range


def create_app(df_mean: pd.DataFrame, df_sum: pd.DataFrame, results_df_mean: pd.DataFrame,
               results_df_sum: pd.DataFrame, sample_cols, ci_df: pd.DataFrame, ci_columns: list,
               global_col_mean: str, global_col_sum: str, isoforms_by_gene, gene_names: dict,
               has_sum: bool,
               geometry_dir: str = None,
               protein_sequences: dict = None,
               domain_mapping: dict = None,
               default_ranking: str = "spearman"):
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

    Returns:
        Configured Dash app
    """
    protein_sequences = protein_sequences or {}
    domain_mapping = domain_mapping or {}
    app = Dash(__name__)
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

    # Prepare table data for both datasets
    table_df_mean = prepare_table_data(results_df_mean, isoforms_by_gene, gene_names,
                                       af_geometry_mapping=af_geometry_mapping,
                                       default_ranking=default_ranking,
                                       protein_sequences=protein_sequences,
                                       domain_mapping=domain_mapping)
    table_df_sum = (prepare_table_data(results_df_sum, isoforms_by_gene, gene_names,
                                       af_geometry_mapping=af_geometry_mapping,
                                       default_ranking=default_ranking,
                                       protein_sequences=protein_sequences,
                                       domain_mapping=domain_mapping)
                    if has_sum else None)
    
    # Pre-compute min_spearman and rank for both datasets
    results_df_mean["min_spearman"] = compute_min_spearman_per_gene(results_df_mean)
    results_df_mean["rank"] = compute_gene_ranking(results_df_mean)
    results_df_mean["rank_by_expression"] = compute_gene_ranking_by_expression(results_df_mean)
    
    if has_sum:
        results_df_sum["min_spearman"] = compute_min_spearman_per_gene(results_df_sum)
        results_df_sum["rank"] = compute_gene_ranking(results_df_sum)
        results_df_sum["rank_by_expression"] = compute_gene_ranking_by_expression(results_df_sum)

    app.layout = html.Div([
        dcc.Store(id="selected-domain"),
        # Dataset toggle at the top
        html.Div([
            html.H3("Isoform Entropy Dashboard", style={**global_style, "marginBottom": "10px"}),
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
        
        # Split view container
        html.Div([
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
                        options=[{'label': ' Only genes with 3D structure', 'value': 'filter_3d'}],
                        value=[],
                        style={'display': 'inline-block', 'marginRight': '20px'}
                    ),
                    dcc.Checklist(
                        id='domain-filter-toggle',
                        options=[{'label': ' Only genes with domains', 'value': 'filter_domains'}],
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
        ], style={"display": "flex", "width": "100%"}),
        
        dcc.Store(id="selected-gene"),
        dcc.Store(id="selected-transcript"),
        dcc.Store(id="selected-exon"),
    ], style={"margin": "20px", **global_style})


    # Register all callbacks
    _register_callbacks(app, isoforms_by_gene, df_mean, df_sum,
                       results_df_mean, results_df_sum,
                       table_df_mean, table_df_sum,
                       ci_df, ci_columns, global_col_mean, global_col_sum, sample_cols,
                       af_geometry_mapping=af_geometry_mapping,
                       protein_sequences=protein_sequences,
                       domain_mapping=domain_mapping,
                       gene_names=gene_names)

    return app


def _register_callbacks(app, isoforms_by_gene, df_mean, df_sum,
                       results_df_mean, results_df_sum,
                       table_df_mean, table_df_sum,
                       ci_df, ci_columns, global_col_mean, global_col_sum, sample_cols,
                       af_geometry_mapping=None, protein_sequences=None,
                       domain_mapping=None,
                       gene_names=None):
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
    protein_sequences = protein_sequences or {}
    domain_mapping = domain_mapping or {}
    gene_names = gene_names or {}

    _genes_with_cds = set()
    if isoforms_by_gene:
        _genes_with_cds = {gene_id for gene_id in isoforms_by_gene.keys()
                           if gene_has_cds(gene_id, isoforms_by_gene)}

    _genes_with_3d = set()
    if isoforms_by_gene and af_geometry_mapping:
        for gene_id, transcripts in isoforms_by_gene.items():
            for transcript_id in transcripts.keys():
                key = transcript_id.replace(".", "").lower()
                if key in af_geometry_mapping:
                    _genes_with_3d.add(gene_id)
                    break

    _genes_with_domains = set()
    if isoforms_by_gene and domain_mapping:
        for gene_id, transcripts in isoforms_by_gene.items():
            for transcript_id in transcripts.keys():
                if domain_mapping.get(transcript_id):
                    _genes_with_domains.add(gene_id)
                    break

    def _axis_ranges_for_results(results_df: pd.DataFrame):
        if results_df is None or results_df.empty:
            return None, None
        x_min = results_df["summed_isoform_entropy"].min()
        x_max = results_df["summed_isoform_entropy"].max()
        y_min = results_df["top_isoform_entropy"].min()
        y_max = results_df["top_isoform_entropy"].max()
        return [x_min, x_max], [y_min, y_max]

    _scatter_axis_ranges = {
        "mean": _axis_ranges_for_results(results_df_mean),
        "sum": _axis_ranges_for_results(results_df_sum),
    }

    # Pre-convert table data to list-of-dicts once (server-side, never sent to browser)
    _table_mean_records = table_df_mean.to_dict('records') if table_df_mean is not None else []
    _table_sum_records  = table_df_sum.to_dict('records')  if table_df_sum  is not None else []

    # Pre-convert CI DataFrame to a plain dict keyed by transcript_id so callbacks
    # never need to round-trip it through the browser.
    _ci_dict = ci_df.to_dict('index') if ci_df is not None else {}

    # Flat transcript → sorted-exon-list lookup (avoids nested loop in every callback)
    _transcript_to_exons: dict = {}
    for _gene_exons in isoforms_by_gene.values():
        for _tid, _exons in _gene_exons.items():
            _transcript_to_exons[_tid] = sorted(_exons, key=lambda e: e["exon_start"])
    
    @app.callback(
        [Output("gene-table", "columns"),
         Output("gene-table", "data")],
        [Input("dataset-toggle", "value"),
         Input("cds-filter-toggle", "value"),
         Input("negative-spearman-toggle", "value"),
         Input("3d-filter-toggle", "value"),
         Input("domain-filter-toggle", "value")],
    )
    def update_table(dataset, cds_filter, negative_spearman_filter, three_d_filter, domain_filter):
        """Update table based on selected dataset, CDS filter, 3D filter, and domain filter."""
        table_data = _table_mean_records if dataset == 'mean' else _table_sum_records

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
        options = [{'label': ' Only genes with CDS', 'value': 'filter_cds', 'disabled': forced_on}]
        if forced_on:
            return ['filter_cds'], options
        return current_cds or [], options

    @app.callback(
        Output("selected-transcript", "data"),
        [Input("isoform-distribution", "clickData"),
         Input("exon-visualization", "clickData"),
         Input("selected-gene", "data")],
        [State("selected-transcript", "data")]
    )
    def update_selected_transcript(dist_click, exon_click, selected_gene, current_transcript):
        """Update selected transcript from distribution or exon visualization clicks.

        Resets to None whenever the selected gene changes so the 3D viewer is
        cleared rather than showing a transcript from the previous gene.
        """
        if not callback_context.triggered:
            return current_transcript

        trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]

        # Gene changed → clear transcript selection unconditionally
        if trigger_id == "selected-gene":
            return None

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
         Input("gene-table", "derived_virtual_selected_rows")],
        [State("selected-gene", "data"),
         State("gene-table", "derived_virtual_data")]
    )
    def update_selected_gene(click_data, selected_rows, current_gene, table_data):
        """Update selected gene from scatter plot click or table selection."""
        # Determine which input triggered the callback
        if not callback_context.triggered:
            return current_gene
        
        trigger_id = callback_context.triggered[0]["prop_id"].split(".")[0]
        
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
        if cds_filter and 'filter_cds' in cds_filter and _genes_with_cds:
            filtered = filtered[filtered["gene_id"].isin(_genes_with_cds)]
        if negative_spearman_filter and 'negative_only' in negative_spearman_filter:
            filtered = filtered[filtered["min_spearman"] <= 0]
        if three_d_filter and 'filter_3d' in three_d_filter:
            if _genes_with_3d:
                filtered = filtered[filtered["gene_id"].isin(_genes_with_3d)]
            else:
                filtered = filtered.iloc[0:0]
        if domain_filter and 'filter_domains' in domain_filter:
            if _genes_with_domains:
                filtered = filtered[filtered["gene_id"].isin(_genes_with_domains)]
            else:
                filtered = filtered.iloc[0:0]
        return filtered

    @functools.lru_cache(maxsize=8)
    def _base_scatter(dataset: str, cds_only: bool, negative_only: bool,
                      three_d_only: bool, domains_only: bool):
        results_df = results_df_mean if dataset == 'mean' else results_df_sum
        cds_filter = ['filter_cds'] if cds_only else []
        negative_filter = ['negative_only'] if negative_only else []
        three_d_filter = ['filter_3d'] if three_d_only else []
        domain_filter = ['filter_domains'] if domains_only else []
        filtered = _filtered_results(results_df, cds_filter, negative_filter, three_d_filter, domain_filter)
        x_range, y_range = _scatter_axis_ranges.get(dataset, (None, None))
        return fig_summed_vs_top_entropy_colored_by_min_spearman(
            filtered,
            selected_gene=None,
            gene_names=gene_names,
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
         Input("domain-filter-toggle", "value")],
    )
    def update_scatter_highlight(selected_gene, dataset, cds_filter,
                                 negative_spearman_filter, three_d_filter, domain_filter):
        """Update scatter plot to highlight selected gene.

        The expensive base figure is cached per dataset; only the gold selection
        ring is added on top, so clicking a gene is near-instant.
        """
        import copy
        cds_only = bool(cds_filter and 'filter_cds' in cds_filter)
        negative_only = bool(negative_spearman_filter and 'negative_only' in negative_spearman_filter)
        three_d_only = bool(three_d_filter and 'filter_3d' in three_d_filter)
        domains_only = bool(domain_filter and 'filter_domains' in domain_filter)
        fig = copy.deepcopy(_base_scatter(dataset, cds_only, negative_only, three_d_only, domains_only))
        if not selected_gene:
            return fig

        results_df = results_df_mean if dataset == 'mean' else results_df_sum
        results_df = _filtered_results(results_df, cds_filter, negative_spearman_filter,
                                       three_d_filter, domain_filter)
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
        df = df_mean if dataset == 'mean' else df_sum
        global_col = global_col_mean if dataset == 'mean' else global_col_sum
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
            isoforms_by_gene,
            transcript_order=transcript_order,
            selected_transcript=None,      # no highlight in cached version
            af_geometry_mapping=af_geometry_mapping,
            domain_data=domain_mapping,
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
            return go.Figure().update_layout(title=None), "No gene selected"
        import copy

        def _exon_label_for_gene(target_gene_id: str) -> str:
            isoforms = isoforms_by_gene.get(target_gene_id, {})
            gene_name = gene_names.get(target_gene_id, "")
            gene_display = gene_name if gene_name else target_gene_id
            if not isoforms:
                return f"Exon structure of gene: {gene_display}"
            gene_start = min(exon["exon_start"] for tx in isoforms.values() for exon in tx)
            gene_end = max(exon["exon_end"] for tx in isoforms.values() for exon in tx)
            gene_length = gene_end - gene_start + 1
            has_domains = any(domain_mapping.get(tid) for tid in isoforms.keys())
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

        df = df_mean if dataset == 'mean' else df_sum
        global_col = global_col_mean if dataset == 'mean' else global_col_sum
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

        gene_name = gene_names.get(selected_gene, "")
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
            return empty_fig, "No gene selected", {'display': 'none'}, []

        # Use server-side data instead of browser storage
        df = df_mean if dataset == 'mean' else df_sum
        has_ci = bool(_ci_dict and ci_columns)  # CIs available if we have the data
        global_col = global_col_mean if dataset == 'mean' else global_col_sum

        # Determine if CIs should be shown based on toggle.
        # In sum mode, CIs are always disabled.
        if dataset == 'sum':
            show_ci = False
        else:
            show_ci = has_ci and ('show_ci' in ci_toggle_value)

        # Determine CI mode from available columns and sample columns
        ci_mode = None
        if show_ci and ci_columns:
            # CI columns are named: ci_lower_{group_name}, ci_upper_{group_name}
            # Extract all group names from CI columns (e.g., "Bipolar" from "ci_lower_Bipolar")
            ci_group_names = set()
            for col in ci_columns:
                # Format: ci_lower_{name} or ci_upper_{name}
                parts = col.split('_', 2)  # Split only on first 2 underscores
                if len(parts) == 3:
                    ci_group_names.add(parts[2])
            
            # Match CI group names against sample columns
            # If sample columns are condition names, they'll match the condition CI groups
            # If sample columns are region names, they'll match the region CI groups
            matching_groups = ci_group_names & set(sample_cols)
            
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
            gene_id, df, sample_cols, show_ci, global_col,
            selected_transcript, _ci_dict, ci_mode, sample_cols
        )

        gene_name = gene_names.get(gene_id, "")
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

        geo = resolve_alphafold_geometry(selected_transcript, af_geometry_mapping)
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

        if not af_geometry_mapping or not selected_gene:
            return html.P(
                f"No 3D model available for {selected_transcript}",
                style={"textAlign": "center", "color": "#999", "padding": "20px"},
            )

        geo = resolve_alphafold_geometry(selected_transcript, af_geometry_mapping)

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
        if not selected_transcript or not protein_sequences:
            return [
                html.Div([
                    html.Div([
                        html.Span("Amino Acid Sequence", style={"fontWeight": "bold"}),
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

        if selected_transcript not in protein_sequences:
            return [
                html.Div([
                    html.Div([
                        html.Span("Amino Acid Sequence", style={"fontWeight": "bold"}),
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

        sequence = protein_sequences[selected_transcript]

        # --- Determine highlight range (1-based, inclusive) ---
        hi_start = hi_end = None  # None → full structure (no highlight)
        exon_label = None

        if (selected_exon is not None and
                selected_exon.get("transcript_id") == selected_transcript):
            exon_idx = selected_exon.get("exon_idx")
            if exon_idx is not None:
                exons_sorted = _transcript_to_exons.get(selected_transcript)
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
            spans = []
            for char_i, aa in enumerate(chunk):
                res_1based = chunk_start + char_i + 1  # 1-based residue number
                if hi_start is not None and hi_start <= res_1based <= hi_end:
                    spans.append(html.Span(
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
                    spans.append(html.Span(aa, style={"color": color}))
            rows.append(html.Div(
                spans,
                style={"marginBottom": "4px", "letterSpacing": "0.5px",
                       "lineHeight": "1.6"},
            ))

        # Join all rows for display and for copying
        display_seq = "\n".join(["".join([span.props['children'] if hasattr(span, 'props') else span.children for span in row.children]) if hasattr(row, 'children') else "" for row in rows])
        copy_seq = sequence

        return [
            html.Div([
                html.Div([
                    html.Span(
                        f"Amino Acid Sequence ({len(sequence)} residues)",
                        style={"fontWeight": "bold"},
                    ),
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
                    children=copy_seq
                )
            ])
        ]

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

        return html.Div([
            html.P([
                html.Span(selected_transcript,
                          style={"fontWeight": "bold", "fontFamily": "Arial, sans-serif"}),
                html.Span(f"  ·  {subtitle}",
                          style={"fontWeight": "bold", "fontFamily": "Arial, sans-serif",
                                 "marginLeft": "6px"}),
            ], style={"marginBottom": "8px"}),
            html.Div(rows),
        ])