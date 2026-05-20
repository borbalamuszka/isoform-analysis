"""Scatter plot visualizations for isoform analysis."""

import plotly.express as px
import plotly.graph_objects as go


def fig_summed_vs_top_entropy_colored_by_min_spearman(results_df, selected_gene=None, gene_names=None):
    """Create scatter plot of summed vs top isoform entropy.

    Assumes min_spearman and rank are pre-computed in results_df.
    """
    if results_df.empty:
        return go.Figure()

    df_plot = results_df[["gene_id", "top_isoform_entropy", "summed_isoform_entropy",
                          "n_isoforms", "min_spearman", "rank"]].copy()
    df_plot = df_plot.dropna(subset=["min_spearman"])

    if df_plot.empty:
        return go.Figure()

    # Add selection indicator
    df_plot["is_selected"] = df_plot["gene_id"] == selected_gene

    # Map gene_id -> gene_name (if provided) so hover labels can show the human
    # readable gene name while keeping Gene ID available in customdata for
    # callbacks (clicks still return gene_id).
    gene_names = gene_names or {}
    df_plot["Gene Name"] = df_plot["gene_id"].apply(lambda g: gene_names.get(g, g))

    df_plot.rename(columns={
        "gene_id": "Gene ID",
        "n_isoforms": "Number of Isoforms",
        "min_spearman": "Min Spearman",
        "top_isoform_entropy": "Top Isoform Entropy",
        "summed_isoform_entropy": "Summed Isoform Entropy",
        "rank": "Rank"
    }, inplace=True)

    colorscale = [[0.0, "red"], [0.5, "grey"], [1.0, "green"]]

    fig = px.scatter(
        df_plot,
        x="Summed Isoform Entropy",
        y="Top Isoform Entropy",
        color="Min Spearman",
        size="Number of Isoforms",
        text="Rank",
        # Prefer showing the gene name when available; fall back to Gene ID.
        hover_name=("Gene Name" if "Gene Name" in df_plot.columns else "Gene ID"),
        hover_data=["Number of Isoforms", "Min Spearman", "Rank"],
        custom_data=["Gene ID"],
        color_continuous_scale=colorscale,
        range_color=(-1.0, 1.0),
    )
    fig.update_traces(textposition='top center', textfont=dict(size=9))
    fig.update_traces(hoverlabel=dict(bgcolor="#E8F4FF", font_color="black"))

    # Highlight selected gene
    if selected_gene:
        selected_rows = df_plot[df_plot["Gene ID"] == selected_gene]
        if not selected_rows.empty:
            fig.add_trace(go.Scatter(
                x=selected_rows["Summed Isoform Entropy"],
                y=selected_rows["Top Isoform Entropy"],
                mode='markers',
                marker=dict(
                    size=20,
                    color='rgba(255, 215, 0, 0)',  # Transparent fill
                    line=dict(color='gold', width=3)
                ),
                showlegend=False,
                hoverinfo='skip',
                name='Selected'
            ))

    fig.update_layout(margin=dict(l=40, r=20, t=40, b=40))
    return fig
