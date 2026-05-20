"""Domain rendering helpers for exon visualizations."""

import plotly.graph_objects as go

from .config import Colors, Dimensions
from .visualization_utils import has_cds

def build_gene_domain_color_map(domain_data, transcript_ids):
    """Assign stable colors to domain names within a single gene."""
    if not domain_data or not transcript_ids:
        return {}

    palette = Colors.DOMAIN_NAME_COLORS
    domain_names = []
    for transcript_id in transcript_ids:
        domains = domain_data.get(transcript_id)
        if not domains:
            continue
        for domain in domains:
            name = domain.get("name")
            if name and name not in domain_names:
                domain_names.append(name)

    color_map = {
        name: palette[idx % len(palette)]
        for idx, name in enumerate(sorted(domain_names))
    }

    return color_map


def add_domain_traces(
    fig,
    domains,
    exons,
    map_coord,
    transcript_id,
    y_pos,
    exon_bar_height,
    row_height,
    domain_color_map,
):
    """Render domain rectangles clipped to CDS parts of overlapping exons."""
    if not domains:
        return

    exon_y_top = y_pos + (row_height + exon_bar_height) / 2
    domain_bar_height = min(Dimensions.DOMAIN_BAR_HEIGHT, exon_bar_height * 0.38)
    domain_gap = 0.0
    domain_y_bottom = exon_y_top + domain_gap
    domain_y_top = domain_y_bottom + domain_bar_height
    domain_name_to_color = domain_color_map

    for domain in domains:
        domain_color = domain_name_to_color[domain["name"]]
        domain_type = domain.get("type", "DOMAIN")
        interpro_id = domain.get("interpro_id") or ""
        dom_min = min(domain["genomic_start"], domain["genomic_end"])
        dom_max = max(domain["genomic_start"], domain["genomic_end"])
        domain_segments = []

        for exon in exons:
            if not has_cds(exon):
                continue

            overlap_start = max(dom_min, exon["cds_start"])
            overlap_end = min(dom_max, exon["cds_end"])
            if overlap_end <= overlap_start:
                continue

            dom_x0 = map_coord(overlap_start)
            dom_x1 = map_coord(overlap_end)
            domain_segments.append((dom_x0, dom_x1))

            description = domain.get('description') or 'N/A'
            hover_text = (
                f"<b>{domain['name']}</b><br>"
                f"{description}<br>"
                f"Type: {domain_type}<br>"
                f"Position: {domain['genomic_start']:,} - {domain['genomic_end']:,} bp<br>"
                f"AA: {domain['aa_start']}-{domain['aa_end']}<br>"
                f"Accession: {domain.get('accession', 'N/A')}<br>"
                f"Library: {domain.get('library', 'N/A')}"
            )

            domain_customdata = [
                {
                    "interpro_id": interpro_id,
                    "kind": "domain",
                    "transcript_id": transcript_id,
                    "name": domain.get("name", ""),
                    "description": domain.get("description", ""),
                    "accession": domain.get("accession", ""),
                    "library": domain.get("library", ""),
                    "aa_start": domain.get("aa_start"),
                    "aa_end": domain.get("aa_end"),
                    "genomic_start": domain.get("genomic_start"),
                    "genomic_end": domain.get("genomic_end"),
                }
                for _ in range(5)
            ]

            fig.add_trace(go.Scatter(
                x=[dom_x0, dom_x1, dom_x1, dom_x0, dom_x0],
                y=[domain_y_bottom, domain_y_bottom, domain_y_top, domain_y_top, domain_y_bottom],
                fill='toself',
                fillcolor=domain_color,
                line=dict(color=domain_color, width=0.8),
                mode='lines',
                hoveron='fills',
                text=hover_text,
                hoverinfo='text',
                customdata=domain_customdata,
                showlegend=False,
                legendgroup=transcript_id
            ))

        if len(domain_segments) > 1:
            domain_segments.sort(key=lambda seg: seg[0])
            y_mid = (domain_y_top + domain_y_bottom) / 2
            for (prev_start, prev_end), (next_start, next_end) in zip(
                domain_segments, domain_segments[1:]
            ):
                fig.add_trace(go.Scatter(
                    x=[prev_end, next_start],
                    y=[y_mid, y_mid],
                    mode='lines',
                    line=dict(color=domain_color, width=1, dash='5px,3px'),
                    hoverinfo='skip',
                    showlegend=False,
                    legendgroup=transcript_id
                ))
