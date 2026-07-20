"""Isoform distribution panel visualizations."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from isoform_distribution.utils import parse_sample_name
from .config import Colors, Dimensions


def fig_isoform_sample_panels(gene_id: str, df: pd.DataFrame, sample_cols, has_ci: bool,
                              global_col: str, selected_transcript=None, ci_data=None,
                              ci_mode=None, all_sample_cols=None):
    """Create isoform distribution panels with optional transcript highlighting and CIs.

    Optimized for faster rendering.
    """
    if df is None or df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No expression data loaded.<br>Configure a Mean/Sum Expression TSV file to view the isoform expression plot.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=12, color="#7f8c8d")
        )
        return fig

    sub = df[df["gene_id"] == gene_id]
    if sub.empty:
        return go.Figure().update_layout(title=f"No isoforms found for {gene_id}")

    max_display = 120
    truncated = len(sub) > max_display
    if truncated:
        sub = sub.nlargest(max_display, global_col)

    sub = sub.sort_values(global_col, ascending=False)
    transcripts = sub["transcript_id"].tolist()
    # Use transcript IDs directly as bar x-labels (the sort order is preserved).
    global_vals = sub[global_col].astype(float).tolist()

    # Pre-compute colors to avoid repeated conditionals
    bg_colors = [Colors.TRANSCRIPT_HIGHLIGHT_BG if t == selected_transcript else 'lightgrey'
                 for t in transcripts]
    bg_line_colors = [Colors.TRANSCRIPT_HIGHLIGHT if t == selected_transcript else 'lightgrey'
                      for t in transcripts]
    bg_line_widths = [2 if t == selected_transcript else 0 for t in transcripts]
    bar_colors = [Colors.TRANSCRIPT_HIGHLIGHT if t == selected_transcript else '#1f77b4' for t in transcripts]
    bar_line_widths = [2 if t == selected_transcript else 0 for t in transcripts]

    fig = make_subplots(
        rows=1,
        cols=len(sample_cols),
        shared_yaxes=True,
        horizontal_spacing=0.02,
        subplot_titles=[s for s in sample_cols]
    )

    for i, sample in enumerate(sample_cols, start=1):
        sample_vals = sub[sample].astype(float).tolist()

        # Background bars
        fig.add_trace(
            go.Bar(
                x=transcripts,
                y=global_vals,
                marker=dict(
                    color=bg_colors,
                    line=dict(color=bg_line_colors, width=bg_line_widths)
                ),
                name=global_col,
                customdata=transcripts,
                showlegend=False,
                hovertemplate='Transcript: %{customdata}<br>' + global_col + ': %{y:.2f}<extra></extra>'
            ),
            row=1, col=i
        )

        # Main bars
        fig.add_trace(
            go.Bar(
                x=transcripts,
                y=sample_vals,
                name=sample,
                marker=dict(
                    color=bar_colors,
                    line=dict(color=Colors.TRANSCRIPT_HIGHLIGHT, width=bar_line_widths)
                ),
                customdata=transcripts,
                showlegend=False,
                hovertemplate='Transcript: %{customdata}<br>' + sample + ': %{y:.2f}<extra></extra>'
            ),
            row=1, col=i
        )

        # Only add CIs if needed (expensive operation)
        if has_ci and ci_data and ci_mode:
            # Determine which CI columns to use based on mode and sample
            ci_lower_vals = []
            ci_upper_vals = []

            for transcript in transcripts:
                transcript_ci = ci_data.get(transcript, {})

                if ci_mode == 'global':
                    lower = transcript_ci.get('ci_lower', None)
                    upper = transcript_ci.get('ci_upper', None)
                else:
                    if ci_mode == 'region':
                        region, _ = parse_sample_name(sample)
                        group = region if region else sample
                    else:  # condition
                        _, condition = parse_sample_name(sample)
                        group = condition if condition else sample

                    lower = transcript_ci.get(f'ci_lower_{group}', None)
                    upper = transcript_ci.get(f'ci_upper_{group}', None)

                ci_lower_vals.append(lower if lower is not None else float('nan'))
                ci_upper_vals.append(upper if upper is not None else float('nan'))

            # Build error bar arrays centered on the sample data values
            ci_center = np.array([
                val if not np.isnan(val) else float('nan')
                for val in sample_vals
            ])
            err_plus = np.array([
                hi - val if not (np.isnan(hi) or np.isnan(val)) else float('nan')
                for hi, val in zip(ci_upper_vals, ci_center)
            ])
            err_minus = np.array([
                val - lo if not (np.isnan(lo) or np.isnan(val)) else float('nan')
                for lo, val in zip(ci_lower_vals, ci_center)
            ])

            # Single scatter trace with asymmetric error bars replaces 2*N traces
            fig.add_trace(
                go.Scatter(
                    x=transcripts,
                    y=ci_center.tolist(),
                    mode='markers',
                    marker=dict(symbol='line-ew', size=10, color='black',
                                line=dict(width=2)),
                    error_y=dict(
                        type='data',
                        symmetric=False,
                        array=err_plus.tolist(),
                        arrayminus=err_minus.tolist(),
                        color='black',
                        thickness=1.5,
                        width=0,
                    ),
                    showlegend=False,
                    hovertemplate=(
                        'Transcript: %{customdata[0]}<br>'
                        f'CI [{sample}]: [%{{customdata[1]:.2f}}, %{{customdata[2]:.2f}}]'
                        '<extra></extra>'
                    ),
                    customdata=list(zip(transcripts, ci_lower_vals, ci_upper_vals)),
                ),
                row=1, col=i,
            )

    # Slant tick labels so long transcript IDs are readable without overlap.
    fig.update_xaxes(tickangle=-45, tickfont=dict(size=9))

    fig.update_yaxes(
        title_text="Transcripts Per Million (TPM)",
        title_font=dict(size=16),
        tickfont=dict(size=16),
        row=1,
        col=1,
    )

    fig.update_layout(
        barmode="overlay",
        margin=dict(l=Dimensions.MARGIN_LEFT, r=Dimensions.MARGIN_RIGHT,
                    t=Dimensions.MARGIN_TOP, b=120),
        height=Dimensions.CHART_HEIGHT_ISOFORMS,
        showlegend=False,
        font=dict(family="Arial, sans-serif", size=13, color="#333333"),
    )


    return fig
