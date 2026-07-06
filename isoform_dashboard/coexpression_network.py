import numpy as np
import pandas as pd
import dash_cytoscape as cyto
from dash import html


def generate_network_elements(gene_coexpression=None, gene_coexpression_idx=None,
                              isoform_coexpression=None, isoform_coexpression_idx=None,
                              target_gene=None, expanded_genes=None,
                              top_k_genes=10, top_k_isoforms=10, isoforms_by_gene=None,
                              gene_names=None, df_expression=None, global_col=None,
                              max_isoforms=120, expanded_gene_position=None,
                              threshold=0.1):
    """Generate Cytoscape elements using precomputed sparse matrices.

    If target_gene is provided, a localized subgraph centred on that gene is
    built.  If expanded_genes is provided, those genes are rendered as Cytoscape
    compound nodes containing their isoform children.

    Args:
        gene_coexpression:      Sparse gene-level coexpression matrix.
        gene_coexpression_idx:  Row/col index labels for the gene matrix.
        isoform_coexpression:   Sparse isoform-level coexpression matrix.
        isoform_coexpression_idx: Row/col index labels for the isoform matrix.
        target_gene:            Gene ID to centre the network on.
        expanded_genes:         Set/list of gene IDs currently expanded as
                                compound nodes showing their isoforms.
        top_k_genes:            Budget for gene-level edges to display.
        top_k_isoforms:         Budget for isoform-level edges within an
                                expanded gene.
        isoforms_by_gene:       Dict gene_id → {transcript_id: [exon list]}.
        gene_names:             Dict gene_id → human-readable gene name.
        df_expression:          DataFrame with gene_id as index column and
                                transcript_id as a column, used to filter
                                isoforms by expression level.  May be None.
        global_col:             Column in df_expression holding mean/sum
                                expression per isoform (used for ranking).
        max_isoforms:           Maximum isoforms to show per expanded gene
                                (matches the exon-panel cap, default 120).
        expanded_gene_position: Dictionary {'x': float, 'y': float} of the gene
                                node's position before it was expanded. Used to
                                offset the compound box so it doesn't obscure
                                incident gene-gene edges.
    """
    elements = []
    if expanded_genes is None:
        expanded_genes = set()
    else:
        expanded_genes = set(expanded_genes)

    if gene_coexpression is None or gene_coexpression_idx is None:
        return []

    # ── Clean index labels (handle bytes from numpy serialisation) ───────────
    gene_idx_clean = [
        g.decode('utf-8') if isinstance(g, bytes) else str(g)
        for g in gene_coexpression_idx
    ]

    # Build bidirectional lookup: both versioned ("ENSG.14") and bare ("ENSG")
    # IDs resolve to the same matrix row so callers don't need to strip versions.
    gene_to_idx: dict[str, int] = {}
    for i, g in enumerate(gene_idx_clean):
        gene_to_idx[g] = i
        if '.' in g:
            gene_to_idx[g.split('.')[0]] = i

    gene_nodes: set[str] = set()
    gene_edges: list[dict] = []
    actual_target_gene: str | None = None

    # ── Build gene-level network ─────────────────────────────────────────────
    if target_gene:
        target_idx = gene_to_idx.get(target_gene)
        if target_idx is None:
            target_idx = gene_to_idx.get(target_gene.split('.')[0])

        if target_idx is not None:
            actual_target_gene = gene_idx_clean[target_idx]

            row = gene_coexpression.getrow(target_idx).toarray()[0]
            abs_row = np.abs(row)
            abs_row[target_idx] = 0  # exclude self-correlation

            # Select top-K neighbours with absolute correlation >= threshold
            candidate_indices = np.where(abs_row >= threshold)[0]
            sorted_candidates = candidate_indices[np.argsort(abs_row[candidate_indices])[::-1]]
            top_indices = list(sorted_candidates[:top_k_genes])

            selected_indices = [target_idx] + top_indices
            selected_genes = [gene_idx_clean[i] for i in selected_indices]

            # Populate gene_nodes from the candidate list
            gene_nodes.update(selected_genes)

            # Build all pairwise edges within the selected subgraph with threshold >= threshold
            gene_edges = []
            for i in range(len(selected_indices)):
                for j in range(i + 1, len(selected_indices)):
                    idx_i, idx_j = selected_indices[i], selected_indices[j]
                    g_i, g_j = gene_idx_clean[idx_i], gene_idx_clean[idx_j]
                    # skip edges where BOTH endpoints are expanded compound nodes
                    if g_i in expanded_genes and g_j in expanded_genes:
                        continue
                    val = gene_coexpression[idx_i, idx_j]
                    abs_val = abs(float(val))
                    if val != 0 and not np.isnan(val) and abs_val >= threshold:
                        gene_edges.append({
                            'source': g_i,
                            'target': g_j,
                            'weight': float(val),
                            'abs_weight': abs_val,
                            'weight_label': f"{float(val):+.3f}",
                        })

        else:
            # Gene not in matrix – render it as an isolated node (no edges)
            gene_nodes.add(target_gene)

    # ── Emit gene nodes (and isoform children for expanded genes) ────────────
    for gene in gene_nodes:
        is_expanded = gene in expanded_genes
        node_elem = {
            'data': {
                'id': gene,
                'label': gene_names.get(gene, gene) if gene_names else gene,
                'type': 'gene',
                'expanded': is_expanded,
            },
            'classes': 'gene-node' + (' expanded' if is_expanded else '')
        }

        # If this is the expanded gene and we captured its position, calculate
        # a base position offset by 150px. We do NOT set this position on the
        # compound parent node itself (which is invalid/ignored in Cytoscape),
        # but rather distribute the isoform child nodes around it.
        base_pos = None
        if is_expanded and expanded_gene_position:
            base_pos = {
                'x': expanded_gene_position['x'] + 150,
                'y': expanded_gene_position['y'] + 150
            }

        elements.append(node_elem)

        if not is_expanded:
            continue
        if isoform_coexpression is None or isoform_coexpression_idx is None:
            continue
        if not isoforms_by_gene:
            continue

        # Build isoform lookup
        iso_idx_clean = [
            iso.decode('utf-8') if isinstance(iso, bytes) else str(iso)
            for iso in isoform_coexpression_idx
        ]
        iso_to_idx: dict[str, int] = {}
        for i, iso in enumerate(iso_idx_clean):
            iso_to_idx[iso] = i
            if '.' in iso:
                iso_to_idx[iso.split('.')[0]] = i

        # Start with all isoforms in the GTF for this gene
        gtf_isoforms = list(isoforms_by_gene.get(gene, {}).keys())

        # FIX 2 – filter & rank by expression to match the exon-panel cap.
        # The exon panel uses nlargest(max_display, global_col) from the active
        # DataFrame; replicate that here so both widgets show the same isoforms.
        if df_expression is not None and global_col is not None:
            gene_sub = df_expression[df_expression['gene_id'] == gene]
            if not gene_sub.empty:
                if len(gene_sub) > max_isoforms:
                    gene_sub = gene_sub.nlargest(max_isoforms, global_col)
                expressed_ids = set(gene_sub['transcript_id'].tolist())
                # Keep only GTF isoforms that also appear in the expression data,
                # preserving expression-ranked order.
                ordered = gene_sub.sort_values(global_col, ascending=False)['transcript_id'].tolist()
                gtf_set  = set(gtf_isoforms)
                gtf_isoforms = [t for t in ordered if t in gtf_set]

        # Determine dynamic size constraints based on the number of isoforms to prevent
        # the compound node box from sprawling too large.
        num_isos = len(gtf_isoforms)
        if num_isos > 20:
            node_size = '7px'
            font_size = '4.5px'
            padding = '3px'
        elif num_isos > 8:
            node_size = '10px'
            font_size = '6px'
            padding = '5px'
        else:
            node_size = '14px'
            font_size = '7px'
            padding = '8px'

        # Apply the compact padding dynamically to the compound parent gene node style
        node_elem['style'] = {'padding': padding}

        # Resolve each isoform ID to a matrix index to check for edges
        valid_iso_data: list[tuple[str, int]] = []
        for iso in gtf_isoforms:
            idx = iso_to_idx.get(iso) or iso_to_idx.get(iso.split('.')[0])
            if idx is not None:
                valid_iso_data.append((iso_idx_clean[idx], idx))

        # Build list of valid edges first to find which isoforms have links
        valid_indices = [idx for _, idx in valid_iso_data]
        edges_to_render = []
        connected_isoforms = set()
        for i in range(len(valid_indices)):
            for j in range(i + 1, len(valid_indices)):
                val = isoform_coexpression[valid_indices[i], valid_indices[j]]
                abs_val = abs(float(val))
                if val != 0 and not np.isnan(val) and abs_val >= threshold:
                    name_i = iso_idx_clean[valid_indices[i]]
                    name_j = iso_idx_clean[valid_indices[j]]
                    connected_isoforms.add(name_i)
                    connected_isoforms.add(name_j)
                    edges_to_render.append({
                        'source': name_i,
                        'target': name_j,
                        'weight': float(val),
                        'abs_weight': abs_val,
                        'weight_label': f"{float(val):+.3f}",
                    })

        # Select a hub isoform node to connect singletons to.
        # This keeps them tightly clustered together without requiring any locked anchor nodes.
        hub_name = None
        if connected_isoforms:
            hub_name = list(connected_isoforms)[0]
        elif gtf_isoforms:
            hub_name = gtf_isoforms[0]

        # Emit isoform child nodes (always render all above-expression-threshold isoforms,
        # even if they are not in the coexpression matrix and are singletons).
        for idx_iso, iso_name in enumerate(gtf_isoforms):
            is_singleton = iso_name not in connected_isoforms

            iso_elem = {
                'data': {
                    'id': iso_name,
                    'label': iso_name,
                    'parent': gene,
                    'type': 'isoform',
                },
                'classes': 'isoform-node',
                'style': {
                    'width': node_size,
                    'height': node_size,
                    'font-size': font_size,
                }
            }

            # Position child nodes in a tight Fermat (golden angle) spiral layout around
            # the offset center. This provides mathematically optimal packing density.
            if base_pos:
                angle = 2.39996 * idx_iso
                scale_factor = 4.0 if num_isos > 20 else (6.0 if num_isos > 8 else 8.0)
                radius = scale_factor * np.sqrt(idx_iso + 1)
                iso_elem['position'] = {
                    'x': base_pos['x'] + radius * np.cos(angle),
                    'y': base_pos['y'] + radius * np.sin(angle)
                }

            elements.append(iso_elem)

            # If this is a singleton (and we have a hub), connect it to the hub isoform via
            # an invisible spring edge. This pulls the singleton to the hub, keeping the
            # compound box compact, while remaining 100% grabbable and unlocked.
            if is_singleton and hub_name and iso_name != hub_name:
                elements.append({
                    'data': {
                        'source': iso_name,
                        'target': hub_name,
                        'edge_type': 'anchor_link'
                    },
                    'classes': 'anchor-edge',
                    'style': {
                        'width': 1,
                        'opacity': 0.0,
                        'events': 'no'
                    }
                })

        # Emit the isoform-level edges
        for edge in edges_to_render:
            elements.append({
                'data': {
                    'source': edge['source'],
                    'target': edge['target'],
                    'weight': edge['weight'],
                    'abs_weight': edge['abs_weight'],
                    'weight_label': edge.get('weight_label', ''),
                    'edge_type': 'isoform_isoform',
                },
                'classes': 'isoform-edge' + (' positive' if edge['weight'] > 0 else ' negative')
            })

    # ── Emit gene-level edges ─────────────────────────────────────────────────
    for edge in gene_edges:
        elements.append({
            'data': {
                'source': edge['source'],
                'target': edge['target'],
                'weight': edge['weight'],
                'abs_weight': edge['abs_weight'],
                'weight_label': edge.get('weight_label', f"{edge['weight']:+.3f}"),
                'edge_type': 'gene_gene',
            },
            'classes': 'gene-edge' + (' positive' if edge['weight'] > 0 else ' negative')
        })

    return elements

def get_cytoscape_stylesheet():
    return [
        # ── Base node style ───────────────────────────────────────────────────
        {
            'selector': 'node',
            'style': {
                'label': 'data(label)',
                'font-size': '10px',
                'text-valign': 'center',
                'text-halign': 'center',
                'color': '#fff',
                'text-outline-width': 1,
                'text-outline-color': '#555',
            }
        },
        # ── Collapsed gene nodes ──────────────────────────────────────────────
        {
            'selector': '.gene-node',
            'style': {
                'background-color': '#0074D9',
                'width': '40px',
                'height': '40px',
                'font-weight': 'bold',
                'font-size': '12px',
                'z-index': 10,
            }
        },
        # ── Expanded gene nodes (compound / parent) ───────────────────────────
        # background-opacity is the Cytoscape-idiomatic way to control fill
        # transparency; rgba alpha in background-color is silently ignored by
        # Cytoscape.js.  Set opacity low so neighbour nodes/edges show through.
        {
            'selector': '.gene-node.expanded',
            'style': {
                'background-color': '#0074D9',
                'background-opacity': 0.05,
                'border-width': 1.5,
                'border-color': '#0074D9',
                'border-opacity': 0.5,
                'border-style': 'dashed',
                'color': '#0074D9',
                'text-valign': 'top',
                'text-halign': 'center',
                'text-outline-width': 0,
                'font-size': '11px',
                'font-weight': 'bold',
                'padding': '8px',
                'compound-sizing-wrt-labels': 'include',
                'z-index': 1,
            }
        },
        # ── Isoform child nodes ───────────────────────────────────────────────
        # Smaller than gene nodes; semi-transparent so the box stays light.
        {
            'selector': '.isoform-node',
            'style': {
                'background-color': '#27AE60',
                'background-opacity': 0.65,
                'width': '14px',
                'height': '14px',
                'font-size': '7px',
                'color': '#fff',
                'text-outline-width': 1,
                'text-outline-color': '#1B5E20',
                'z-index': 5,
            }
        },

        # ── Gene–gene edges ───────────────────────────────────────────────────
        # Opacity scales with |weight| so strong correlations are visually
        # prominent and weak ones recede without disappearing entirely.
        # z-index 999 ensures gene-gene edges always render on top of the
        # expanded compound-node background (which has z-index 1).
        {
            'selector': '.gene-edge',
            'style': {
                'width': 'mapData(abs_weight, 0, 1, 1.0, 5.0)',
                'opacity': 'mapData(abs_weight, 0, 1, 0.35, 0.90)',
                'z-index': 999,
            }
        },
        # Show the numeric weight on a gene-gene edge when it is selected,
        # so users can hover-click to read the exact correlation value.
        {
            'selector': '.gene-edge:selected',
            'style': {
                'label': 'data(weight_label)',
                'font-size': '9px',
                'color': '#333',
                'text-background-color': '#fffff0',
                'text-background-opacity': 0.85,
                'text-background-padding': '2px',
                'text-background-shape': 'roundrectangle',
                'text-border-color': '#aaa',
                'text-border-width': 0.5,
                'text-border-opacity': 1,
            }
        },
        # ── Isoform–isoform edges ─────────────────────────────────────────────
        # Deliberately lighter (lower max opacity) than gene–gene edges so the
        # isoform sub-network reads as secondary / internal detail.
        {
            'selector': '.isoform-edge',
            'style': {
                'width': 'mapData(abs_weight, 0, 1, 0.5, 3.0)',
                'opacity': 'mapData(abs_weight, 0, 1, 0.15, 0.55)',
                'z-index': 3,
            }
        },
        {
            'selector': '.isoform-edge:selected',
            'style': {
                'label': 'data(weight_label)',
                'font-size': '8px',
                'color': '#333',
                'text-background-color': '#fffff0',
                'text-background-opacity': 0.85,
                'text-background-padding': '2px',
                'text-background-shape': 'roundrectangle',
            }
        },
        # ── Positive / negative correlation colour coding ──────────────────────
        {
            'selector': '.positive',
            'style': {'line-color': '#C0392B'}   # deep red
        },
        {
            'selector': '.negative',
            'style': {'line-color': '#2471A3'}   # deep blue
        },
    ]

def create_coexpression_widget():
    return cyto.Cytoscape(
        id='coexpression-network',
        layout={
            'name': 'cose',
            # Pull nodes closer together for a denser, less sprawling layout
            'idealEdgeLength': 60,
            'nodeOverlap': 10,
            'refresh': 20,
            'fit': True,
            'padding': 20,
            'randomize': False,
            'componentSpacing': 60,
            'nodeRepulsion': 3000,
            'edgeElasticity': 80,
            # Lower nestingFactor keeps compound (expanded) boxes tighter;
            # the default of 5 caused them to push too far from their neighbours.
            'nestingFactor': 0.1,
            'gravity': 100,
            'numIter': 1000,
            'initialTemp': 200,
            'coolingFactor': 0.95,
            'minTemp': 1.0,
        },
        style={
            'width': '100%',
            'height': '600px',
            'backgroundColor': '#f8f9fa',
            'border': '1px solid #ddd',
        },
        elements=[],
        stylesheet=get_cytoscape_stylesheet(),
        minZoom=0.01,
        maxZoom=10,
    )
