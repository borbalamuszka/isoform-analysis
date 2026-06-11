import numpy as np
import pandas as pd
from scipy.stats import rankdata
import dash_cytoscape as cyto
from dash import html

def generate_network_elements(gene_coexpression=None, gene_coexpression_idx=None, 
                              isoform_coexpression=None, isoform_coexpression_idx=None, 
                              target_gene=None, expanded_genes=None, 
                              top_k_genes=10, top_k_isoforms=10, isoforms_by_gene=None):
    """
    Generate Cytoscape elements using precomputed sparse matrices.
    If target_gene is provided, we compute the gene network centered around it.
    If expanded_genes is provided, those genes will be rendered as compound nodes containing their isoforms.
    """
    elements = []
    if expanded_genes is None:
        expanded_genes = set()
    else:
        expanded_genes = set(expanded_genes)
        
    if gene_coexpression is None or gene_coexpression_idx is None:
        return []

    # Clean indices (handle bytes from numpy serialization)
    gene_idx_clean = [g.decode('utf-8') if isinstance(g, bytes) else str(g) for g in gene_coexpression_idx]

    # Map ID to index (handling potential .version mismatches)
    gene_to_idx = {}
    for i, g in enumerate(gene_idx_clean):
        gene_to_idx[g] = i
        if '.' in g:
            gene_to_idx[g.split('.')[0]] = i
    
    gene_nodes = set()
    gene_edges = []
    
    target_idx = None
    actual_target_gene = target_gene

    if target_gene:
        if target_gene in gene_to_idx:
            target_idx = gene_to_idx[target_gene]
        elif target_gene.split('.')[0] in gene_to_idx:
            target_idx = gene_to_idx[target_gene.split('.')[0]]
            
        if target_idx is not None:
            actual_target_gene = gene_idx_clean[target_idx]
        else:
            # Fallback: Matrix doesn't have this gene, but we still render it as an isolated node
            gene_nodes.add(target_gene)

    if target_idx is not None:
        # Localized network
        
        # Get correlations for this gene
        row = gene_coexpression.getrow(target_idx).toarray()[0]
        
        # We need the top K genes correlated to the target
        abs_row = np.abs(row)
        
        # Exclude self
        abs_row[target_idx] = 0
        
        if len(abs_row) <= top_k_genes:
            top_indices = np.arange(len(abs_row))
        else:
            top_indices = np.argpartition(abs_row, -top_k_genes)[-top_k_genes:]
            
        selected_genes = [actual_target_gene]
        for idx in top_indices:
            if abs_row[idx] > 0:
                selected_genes.append(gene_idx_clean[idx])
                
        # Now get the subgraph between all selected genes
        selected_indices = [gene_to_idx[g] for g in selected_genes if g in gene_to_idx]
        
        for i in range(len(selected_indices)):
            for j in range(i + 1, len(selected_indices)):
                idx_i = selected_indices[i]
                idx_j = selected_indices[j]
                val = gene_coexpression[idx_i, idx_j]
                if val != 0 and not np.isnan(val):
                    gene_edges.append({
                        'source': gene_idx_clean[idx_i],
                        'target': gene_idx_clean[idx_j],
                        'weight': float(val),
                        'abs_weight': abs(float(val))
                    })
        
        # Sort and keep top_k_genes edges to prevent clutter
        gene_edges.sort(key=lambda x: x['abs_weight'], reverse=True)
        gene_edges = gene_edges[:top_k_genes]
        
        for edge in gene_edges:
            gene_nodes.add(edge['source'])
            gene_nodes.add(edge['target'])
                    
        gene_nodes.add(actual_target_gene) # Ensure target is always present

    # Add Gene Nodes
    for gene in gene_nodes:
        is_expanded = gene in expanded_genes
        elements.append({
            'data': {
                'id': gene,
                'label': gene,
                'type': 'gene',
                'expanded': is_expanded
            },
            'classes': 'gene-node' + (' expanded' if is_expanded else '')
        })
        
        # Isoform processing
        if is_expanded and isoform_coexpression is not None and isoform_coexpression_idx is not None and isoforms_by_gene:
            iso_idx_clean = [iso.decode('utf-8') if isinstance(iso, bytes) else str(iso) for iso in isoform_coexpression_idx]
            # Find isoforms for this gene
            gene_isoforms = isoforms_by_gene.get(gene, {})
            if gene_isoforms:
                iso_names = list(gene_isoforms.keys())
                
                iso_to_idx = {}
                for i, iso in enumerate(iso_idx_clean):
                    iso_to_idx[iso] = i
                    if '.' in iso:
                        iso_to_idx[iso.split('.')[0]] = i
                
                valid_iso_indices = []
                valid_iso_names = []
                for iso in iso_names:
                    if iso in iso_to_idx:
                        idx = iso_to_idx[iso]
                        valid_iso_indices.append(idx)
                        valid_iso_names.append(iso_idx_clean[idx])
                    elif iso.split('.')[0] in iso_to_idx:
                        idx = iso_to_idx[iso.split('.')[0]]
                        valid_iso_indices.append(idx)
                        valid_iso_names.append(iso_idx_clean[idx])
                
                # We need the top K edges within this isoform set
                iso_edges = []
                for i in range(len(valid_iso_indices)):
                    for j in range(i + 1, len(valid_iso_indices)):
                        idx_i = valid_iso_indices[i]
                        idx_j = valid_iso_indices[j]
                        val = isoform_coexpression[idx_i, idx_j]
                        if val != 0 and not np.isnan(val):
                            iso_edges.append({
                                'source': valid_iso_names[i],
                                'target': valid_iso_names[j],
                                'weight': float(val),
                                'abs_weight': abs(float(val))
                            })
                
                # Sort and keep top_k_isoforms
                iso_edges.sort(key=lambda x: x['abs_weight'], reverse=True)
                iso_edges = iso_edges[:top_k_isoforms]
                
                # Always add all valid isoforms as nodes so they appear
                for iso in valid_iso_names:
                    elements.append({
                        'data': {
                            'id': iso,
                            'label': iso,
                            'parent': gene,
                            'type': 'isoform'
                        },
                        'classes': 'isoform-node'
                    })
                    
                for edge in iso_edges:
                    elements.append({
                        'data': {
                            'source': edge['source'],
                            'target': edge['target'],
                            'weight': edge['weight'],
                            'abs_weight': edge['abs_weight'],
                            'edge_type': 'isoform_isoform'
                        },
                        'classes': 'isoform-edge' + (' positive' if edge['weight'] > 0 else ' negative')
                    })
                    
    # Add Gene Edges
    for edge in gene_edges:
        elements.append({
            'data': {
                'source': edge['source'],
                'target': edge['target'],
                'weight': edge['weight'],
                'abs_weight': edge['abs_weight'],
                'edge_type': 'gene_gene'
            },
            'classes': 'gene-edge' + (' positive' if edge['weight'] > 0 else ' negative')
        })

    return elements

def get_cytoscape_stylesheet():
    return [
        {
            'selector': 'node',
            'style': {
                'label': 'data(label)',
                'font-size': '10px',
                'text-valign': 'center',
                'text-halign': 'center',
                'color': '#fff',
                'text-outline-width': 1,
                'text-outline-color': '#888'
            }
        },
        {
            'selector': '.gene-node',
            'style': {
                'background-color': '#0074D9',
                'width': '40px',
                'height': '40px',
                'font-weight': 'bold',
                'font-size': '12px'
            }
        },
        {
            'selector': '.gene-node.expanded',
            'style': {
                'background-color': 'rgba(0, 116, 217, 0.1)',
                'border-width': 2,
                'border-color': '#0074D9',
                'color': '#333',
                'text-valign': 'top',
                'text-outline-width': 0
            }
        },
        {
            'selector': '.isoform-node',
            'style': {
                'background-color': '#2ECC40',
                'width': '20px',
                'height': '20px',
                'font-size': '8px'
            }
        },
        {
            'selector': 'edge',
            'style': {
                'width': 'mapData(abs_weight, 0, 1, 1, 5)',
                'opacity': 0.8
            }
        },
        {
            'selector': '.positive',
            'style': {
                'line-color': '#FF4136'
            }
        },
        {
            'selector': '.negative',
            'style': {
                'line-color': '#0074D9'
            }
        }
    ]

def create_coexpression_widget():
    return cyto.Cytoscape(
        id='coexpression-network',
        layout={'name': 'cose', 'idealEdgeLength': 100, 'nodeOverlap': 20, 'refresh': 20, 'fit': True, 'padding': 30, 'randomize': False, 'componentSpacing': 100, 'nodeRepulsion': 4000, 'edgeElasticity': 100, 'nestingFactor': 5, 'gravity': 80, 'numIter': 1000, 'initialTemp': 200, 'coolingFactor': 0.95, 'minTemp': 1.0},
        style={'width': '100%', 'height': '600px', 'backgroundColor': '#f9f9f9', 'border': '1px solid #ddd'},
        elements=[],
        stylesheet=get_cytoscape_stylesheet(),
        minZoom=0.01,
        maxZoom=10
    )
