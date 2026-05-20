import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import defaultdict
import os

INPUT_FILE = "data/neuro_project/expressed_isoforms.positions.txt"
OUTPUT_DIR = "data/neuro_project/output/exons"

def parse_isoform_file(filename):
    """
    Parse the expressed isoforms positions file.
    Returns a dictionary with isoform data grouped by gene.
    """
    isoforms_by_gene = defaultdict(lambda: defaultdict(list))
    
    with open(filename, 'r') as f:
        header = f.readline()  # Skip header
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 12:
                continue
            
            chrom = parts[0]
            exon_start = int(parts[1])
            exon_end = int(parts[2])
            try:
                cds_start = int(parts[5]) if parts[5] else None
                cds_end = int(parts[6]) if parts[6] else None
            except Exception:
                cds_start = None
                cds_end = None
            gene_id = parts[9]
            transcript_id = parts[10]
            
            isoforms_by_gene[gene_id][transcript_id].append({
                "exon_start": exon_start,
                "exon_end": exon_end,
                "cds_start": cds_start,
                "cds_end": cds_end
            })
    
    return isoforms_by_gene

def get_gene_boundaries(isoforms):
    """Get the start and end coordinates for all isoforms of a gene."""
    gene_start = float('inf')
    gene_end = 0
    
    for transcript_id, exons in isoforms.items():
        for exon in exons:
            gene_start = min(gene_start, exon["exon_start"])
            gene_end = max(gene_end, exon["exon_end"])
    
    return gene_start, gene_end

def visualize_isoforms(gene_id, isoforms_by_gene, output_file):
    """
    Visualize all isoforms of a gene as binary patterns with CDS regions highlighted.
    """
    if gene_id not in isoforms_by_gene:
        print(f"Error: No isoforms found for gene {gene_id}")
        return
    
    isoforms = isoforms_by_gene[gene_id]
    gene_start, gene_end = get_gene_boundaries(isoforms)
    gene_length = gene_end - gene_start + 1
    
    # Sort isoforms by name for consistent ordering
    sorted_isoforms = sorted(isoforms.items())
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(len(sorted_isoforms) * 0.35, 1.5)))
    
    # Plot each isoform
    for idx, (transcript_id, exons) in enumerate(sorted_isoforms):
        # Sort exons by start position
        exons = sorted(exons, key=lambda e: e["exon_start"])
        
        y_pos = len(sorted_isoforms) - idx - 1
        
        # Draw exons
        for exon in exons:
            rel_start = exon["exon_start"] - gene_start
            rel_end = exon["exon_end"] - gene_start
            width = rel_end - rel_start + 1
            
            # Draw exon rectangle
            rect = patches.Rectangle(
                (rel_start, y_pos), 
                width, 
                0.8,
                linewidth=0,
                edgecolor='none',
                facecolor='#2E86AB',
                alpha=0.9
            )
            ax.add_patch(rect)
            
            # Draw CDS region if available
            if exon["cds_start"] is not None and exon["cds_end"] is not None and exon["cds_end"] > exon["cds_start"]:
                # Convert CDS coordinates from exon-relative to genomic
                cds_genomic_start = exon["exon_start"] + (exon["cds_start"] - 1)
                cds_genomic_end = exon["exon_start"] + (exon["cds_end"] - 1)
                cds_rel_start = cds_genomic_start - gene_start
                cds_rel_end = cds_genomic_end - gene_start
                cds_width = cds_rel_end - cds_rel_start + 1
                
                # Only draw CDS if it falls within the exon boundaries
                if cds_genomic_start >= exon["exon_start"] and cds_genomic_end <= exon["exon_end"]:
                    cds_rect = patches.Rectangle(
                        (cds_rel_start, y_pos + 0.2),
                        cds_width,
                        0.4,
                        linewidth=0,
                        edgecolor='none',
                        facecolor='#F5B041',
                        alpha=1.0
                    )
                    ax.add_patch(cds_rect)
        
        # Add transcript label
        label_offset = gene_length * 0.02
        ax.text(
            -label_offset, 
            y_pos + 0.4, 
            transcript_id, 
            ha='right', 
            va='center', 
            fontsize=6,
            family="Arial"
        )
        
        # Add separator line between isoforms
        if idx < len(sorted_isoforms) - 1:
            ax.axhline(y=y_pos - 0.1, color='#CCCCCC', linewidth=0.5, linestyle='-', alpha=0.5)
    
    # Set limits and remove axes
    ax.set_xlim(-gene_length * 0.02, gene_length)
    ax.set_ylim(-0.5, len(sorted_isoforms))
    ax.axis('off')
    
    # Add title with gene length
    plt.title(f'Gene: {gene_id} ({gene_length:,} bp)', 
              fontsize=7, pad=20, loc='left', fontweight='bold', family="Arial")
    
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"{gene_id}: {len(isoforms)} isoforms, {gene_length:,} bp → {output_file}")
    
    plt.close()

if __name__ == "__main__":
    import sys
    
    # Parse the hardcoded input file
    isoforms_by_gene = parse_isoform_file(INPUT_FILE)
    
    if len(sys.argv) < 2:
        print("Usage: python exon_visualisation.py <gene_id>")
        print(f"\nAvailable genes ({len(isoforms_by_gene)} total):")
        for gene_id in sorted(isoforms_by_gene.keys())[:10]:
            print(f"  {gene_id} ({len(isoforms_by_gene[gene_id])} isoforms)")
        if len(isoforms_by_gene) > 10:
            print(f"  ... and {len(isoforms_by_gene) - 10} more")
        sys.exit(1)
    
    gene_id = sys.argv[1]
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Generate output file path
    output_file = os.path.join(OUTPUT_DIR, f"{gene_id}.png")
    
    # Generate visualization
    visualize_isoforms(gene_id, isoforms_by_gene, output_file)