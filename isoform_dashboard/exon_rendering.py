"""Exon and domain track rendering."""

import plotly.graph_objects as go

from .config import Colors, Dimensions
from .gtf_parser import build_compressed_mapping
from .domain_rendering import add_domain_traces, build_gene_domain_color_map
from .visualization_utils import has_cds


def create_exon_visualization(gene_id, isoforms_by_gene, transcript_order=None,
                              selected_transcript=None,
                              af_geometry_mapping=None,
                              exon_scale_bp_per_unit=None, intron_width_units=None,
                              domain_data=None):
    """Create interactive exon visualization using Plotly.

    Args:
        gene_id: Gene ID to visualize
        isoforms_by_gene: Dictionary mapping gene_id -> transcript_id -> list of exons
        transcript_order: Optional list of transcript IDs to control display order
        selected_transcript: Transcript ID to highlight (optional)
        af_geometry_mapping: Output of build_alphafold_geometry_mapping() (optional).
        exon_scale_bp_per_unit: Scaling factor for exon size
        intron_width_units: Fixed width for introns
        domain_data: Optional dictionary mapping transcript_id -> list of domain dicts.
                     Each domain dict should have keys: name, type, genomic_start, genomic_end.

    Returns:
        Plotly figure
    """
    # Use defaults from config if not provided
    if exon_scale_bp_per_unit is None:
        exon_scale_bp_per_unit = Dimensions.EXON_SCALE_BP_PER_UNIT
    if intron_width_units is None:
        intron_width_units = Dimensions.INTRON_WIDTH_UNITS

    if gene_id not in isoforms_by_gene:
        return go.Figure().update_layout(title=f"No exon data for {gene_id}")

    isoforms = isoforms_by_gene[gene_id]
    af_geometry_mapping = af_geometry_mapping or {}

    # Build compressed coordinate mapping
    mapping = build_compressed_mapping(
        isoforms,
        exon_scale_bp_per_unit=exon_scale_bp_per_unit,
        intron_width_units=intron_width_units
    )
    map_coord = mapping["map_coord"]
    total_width = mapping["total_width"]

    gene_start = min(exon["exon_start"] for tx in isoforms.values() for exon in tx)
    gene_end = max(exon["exon_end"] for tx in isoforms.values() for exon in tx)
    gene_length = gene_end - gene_start + 1

    strand = next(
        (exon.get("strand") for tx in isoforms.values() for exon in tx if exon.get("strand")),
        "+",
    )
    arrow_text = "←" if strand == "-" else "→"

    if transcript_order:
        sorted_isoforms = [(tid, isoforms[tid]) for tid in transcript_order if tid in isoforms]
    else:
        sorted_isoforms = sorted(isoforms.items())

    max_label_chars = 22
    label_texts = []
    for idx, (tid, _) in enumerate(sorted_isoforms):
        label = f"{idx + 1}. {tid}"
        if len(label) > max_label_chars:
            label = f"{label[:max_label_chars - 1]}…"
        label_texts.append(label)

    fig = go.Figure()
    global_domain_colors = build_gene_domain_color_map(
        domain_data,
        [tid for tid, _ in sorted_isoforms],
    )

    exon_bar_height = Dimensions.EXON_BAR_HEIGHT
    domain_bar_height = min(Dimensions.DOMAIN_BAR_HEIGHT, exon_bar_height * 0.38)
    domain_gap = 0.0
    has_domains = bool(domain_data) and any(
        domain_data.get(tid) for tid, _ in sorted_isoforms
    )
    row_height = Dimensions.ROW_HEIGHT
    if has_domains:
        min_row_height = exon_bar_height + 2 * (domain_bar_height + domain_gap)
        row_height = max(row_height, min_row_height)
    n_isoforms = len(sorted_isoforms)

    for idx, (transcript_id, exons) in enumerate(sorted_isoforms):
        exons = sorted(exons, key=lambda e: e["exon_start"])
        y_pos = (n_isoforms - idx - 1) * row_height
        # Determine if this transcript is selected
        is_selected = (transcript_id == selected_transcript)

        # Add background for selected transcripts
        if is_selected:
            # Cyan/teal background for selected
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

        exon_label_annotations = []
        # Collect exon centres for click-target trace (drawn after polygons)
        click_xs: list = []
        click_ys: list = []
        click_custom: list = []
        click_texts: list = []

        UTR_COLOR = Colors.EXON_UTR
        CDS_COLOR = Colors.EXON_CDS

        for exon_idx, exon in enumerate(exons):
            x0 = map_coord(exon["exon_start"])
            x1 = map_coord(exon["exon_end"])

            has_cds_exon = has_cds(exon)

            # Add border for selected transcript
            if is_selected:
                border_color = Colors.TRANSCRIPT_HIGHLIGHT
                border_width = 2
            else:
                border_color = None
                border_width = 0

            y_top    = y_pos + (row_height + exon_bar_height) / 2
            y_bottom = y_pos + (row_height - exon_bar_height) / 2

            def _rect_trace(rx0, rx1, color, bd_color=None, bd_width=0,
                            _yt=y_top, _yb=y_bottom, _tid=transcript_id):
                """Return a filled rectangle trace."""
                return go.Scatter(
                    x=[rx0, rx1, rx1, rx0, rx0],
                    y=[_yb, _yb, _yt, _yt, _yb],
                    fill='toself',
                    fillcolor=color,
                    line=dict(color=bd_color or color, width=bd_width),
                    mode='lines',
                    hoverinfo='skip',
                    showlegend=False,
                    legendgroup=_tid
                )

            if has_cds_exon:
                cds_x0 = map_coord(exon["cds_start"])
                cds_x1 = map_coord(exon["cds_end"])

                # Draw full exon in UTR blue first (UTR overhangs on either side)
                fig.add_trace(_rect_trace(x0, x1, UTR_COLOR,
                                          border_color, border_width))
                # Draw only the CDS region in orange on top
                fig.add_trace(_rect_trace(cds_x0, cds_x1, CDS_COLOR,
                                          border_color, border_width))
            else:
                # Pure UTR exon
                fig.add_trace(_rect_trace(x0, x1, UTR_COLOR,
                                          border_color, border_width))

            # Build hover text
            if has_cds_exon:
                hover_text = (f"<b>{transcript_id}</b><br>"
                              f"Exon {exon_idx + 1} (CDS): "
                              f"{exon['exon_start']:,} - {exon['exon_end']:,}"
                              f"<br>CDS: {exon['cds_start']:,} - {exon['cds_end']:,}")
            else:
                hover_text = (f"<b>{transcript_id}</b><br>"
                              f"Exon {exon_idx + 1} (UTR): "
                              f"{exon['exon_start']:,} - {exon['exon_end']:,}")

            # Exon number label: collect positions to draw on top later
            exon_width = x1 - x0
            cx = x0 + min(0.4, exon_width * 0.2)
            cy = y_pos + row_height / 2
            exon_label_annotations.append((cx, cy, str(exon_idx + 1)))

            # Centre point for click-target trace
            # Add markers at 25 %, 50 %, and 75 % of the exon width so that
            # clicks anywhere along wide CDS exons are within range of a marker.
            for frac in (0.25, 0.5, 0.75):
                click_xs.append(x0 + frac * (x1 - x0))
                click_ys.append(cy)
                click_custom.append([transcript_id])
                click_texts.append(hover_text)

        # One transparent marker layer per exon — Plotly registers clicks on markers.
        # Multiple markers per exon ensure clicks anywhere along wide CDS bars are captured.
        if click_xs:
            fig.add_trace(go.Scatter(
                x=click_xs,
                y=click_ys,
                mode='markers',
                marker=dict(size=24, opacity=0, color='rgba(0,0,0,0)'),
                text=click_texts,
                hoverinfo='text',
                customdata=click_custom,
                showlegend=False,
                legendgroup=transcript_id
            ))

    # Render protein domains above exons if available
        if domain_data and transcript_id in domain_data:
            add_domain_traces(
                fig=fig,
                domains=domain_data[transcript_id],
                exons=exons,
                map_coord=map_coord,
                transcript_id=transcript_id,
                y_pos=y_pos,
                exon_bar_height=exon_bar_height,
                row_height=row_height,
                domain_color_map=global_domain_colors,
            )

        # Draw exon indices after domains so labels stay visible on top.
        for label_x, label_y, label_text in exon_label_annotations:
            fig.add_annotation(
                x=label_x,
                y=label_y,
                text=label_text,
                showarrow=False,
                xanchor='left',
                yanchor='middle',
                font=dict(size=9, color='white', family='Arial'),
            )

        # Add transcript label
        label_offset = total_width * 0.015
        label_color = Colors.TRANSCRIPT_HIGHLIGHT if is_selected else 'black'
        label_weight = 'bold' if is_selected else 'normal'
        label_text = label_texts[idx] if idx < len(label_texts) else f"{idx + 1}. {transcript_id}"

        fig.add_annotation(
            x=-label_offset,
            y=y_pos + row_height / 2,
            text=label_text,
            showarrow=False,
            xanchor='right',
            yanchor='middle',
            font=dict(size=12, family="Arial",
                     weight=label_weight,
                     color=label_color)
        )

    fig.update_layout(
        title=None,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='white',
        height=max(n_isoforms * 40, 200),
        margin=dict(l=0, r=0, t=20, b=20),
        hovermode='closest',
        meta={
            "total_width": float(total_width),
            "row_height": float(row_height),
            "transcript_order": [tid for tid, _ in sorted_isoforms],
            "label_texts": label_texts,
        }
    )

    top_y = n_isoforms * row_height + 0.6
    arrow_label = (
        "Strand direction "
        f"<span style='font-size:16px; font-weight:700;'>{arrow_text}</span>"
    )
    fig.add_annotation(
        x=total_width * 0.5,
        y=top_y,
        text=arrow_label,
        showarrow=False,
        yanchor='bottom',
        font=dict(size=12, color='#333')
    )

    return fig
