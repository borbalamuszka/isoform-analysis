"""Configuration constants for the isoform dashboard.

Centralized place for colors, dimensions, styling, and other magic values.
Avoids hardcoding constants throughout the codebase.
"""


class Colors:
    """Color scheme for visualizations and UI."""
    
    # Exon visualization
    EXON_CDS = '#F5B041'                # Orange for coding exons
    EXON_UTR = '#2E86AB'                # Blue for UTRs
    
    # Selection/highlight
    TRANSCRIPT_HIGHLIGHT = '#00BFFF'    # Deep Sky Blue
    TRANSCRIPT_HIGHLIGHT_BG = 'rgba(0, 191, 255, 0.15)'
    TRANSCRIPT_HIGHLIGHT_LINE = 'rgba(0, 191, 255, 0)'
    
    # Charts
    SCATTER_NEGATIVE = 'red'
    SCATTER_NEUTRAL = 'grey'
    SCATTER_POSITIVE = 'green'
    
    # Domain visualization
    DOMAIN_DEFAULT = '#22C55E'          # Green for general domains
    DOMAIN_BY_TYPE = {
        'DOMAIN': '#22C55E',            # Green
        'FAMILY': '#F97316',            # Orange
        'REGION': '#3B82F6',            # Blue
        'SUPERFAMILY': '#A855F7',       # Purple
    }
    # Colors for specific domain names (cycling palette)
    DOMAIN_NAME_COLORS = [
        '#16A34A',  # Green
        '#BE185D',  # Magenta
        '#7C3AED',  # Purple
        '#0D9488',  # Teal
        '#2563EB',  # Blue
        '#DB2777',  # Pink
        '#1D4ED8',  # Deep blue
        '#6D28D9',  # Indigo
    ]
    
    # UI elements
    BUTTON_SUCCESS = '#4CAF50'
    TEXT_SECONDARY = '#999'
    TEXT_DISABLED = '#cccccc'
    BORDER_LIGHT = '#ddd'
    BG_LIGHT = '#f9f9f9'
    BG_LIGHT_BLUE = '#f0f8ff'


class Dimensions:
    """Layout and sizing constants."""
    
    # Chart heights
    CHART_HEIGHT_ISOFORMS = 450
    CHART_HEIGHT_EXONS = 200  # min height
    EXONS_HEIGHT_PER_ISOFORM = 40
    
    # Exon visualization
    EXON_SCALE_BP_PER_UNIT = 5.0      # How many bp per visual unit (smaller -> wider exons)
    INTRON_WIDTH_UNITS = 10.0           # Fixed width for introns
    EXON_BAR_HEIGHT = 0.6
    DOMAIN_BAR_HEIGHT = 0.35            # Smaller than exons for visual hierarchy
    ROW_HEIGHT = 1.0
    
    # Margins (for plots)
    MARGIN_LEFT = 60
    MARGIN_RIGHT = 30
    MARGIN_TOP = 60
    MARGIN_BOTTOM = 120
    
    # Table
    TABLE_PAGE_SIZE = 20
    TABLE_MIN_WIDTH_PX = 70
    TABLE_MAX_WIDTH_PX = 160
    
    # Data limits
    MAX_ISOFORMS_DISPLAY = 120


class Styles:
    """Common style dictionaries for Dash components."""
    
    BASE = {
        "fontFamily": "Arial, sans-serif",
        "fontSize": "14px",
        "color": "#333333",
    }
    
    HEADER = {
        **BASE,
        "fontWeight": "bold",
        "fontSize": "16px",
    }
    
    HEADER_SMALL = {
        **BASE,
        "fontWeight": "bold",
        "fontSize": "14px",
    }
    
    BUTTON_PRIMARY = {
        'backgroundColor': Colors.BUTTON_SUCCESS,
        'color': 'white',
        'padding': '10px 20px',
        'border': 'none',
        'borderRadius': '4px',
        'cursor': 'pointer',
        'fontSize': '14px',
        'fontWeight': 'bold'
    }
    
    TIP_BOX = {
        'backgroundColor': Colors.BG_LIGHT_BLUE,
        'padding': '10px',
        'borderRadius': '4px',
        'fontSize': '13px',
        'border': f'1px solid #d0e8ff',
    }
    
    BORDER_BOTTOM = {
        "marginBottom": "20px",
        "paddingBottom": "10px",
        "borderBottom": f"2px solid {Colors.BORDER_LIGHT}"
    }
