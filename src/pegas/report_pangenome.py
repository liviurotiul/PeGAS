import os

import pandas as pd
import plotly.graph_objects as go
from plotly.colors import qualitative
from tqdm import tqdm

from utils import configure_pie_chart, fix_colorbar


def create_pangenome_heatmaps(species_dict):
    """
    Generates pangenome heatmaps for each species present in the 'pangenome' directory.

    Parameters:
    - species_dict (dict): Mapping from species codes to common names.

    Returns:
    - figures (dict): Mapping from species to Plotly heatmap figures.
    """
    # Loop over each species directory in 'pangenome'.
    pangenome_dir = "pangenome"
    if not os.path.exists(pangenome_dir):
        tqdm.write(f"Directory '{pangenome_dir}' does not exist. Skipping pangenome heatmaps.")
        return ""

    figures = {}

    for sp in [x for x in os.listdir(pangenome_dir) if os.path.isdir(os.path.join(pangenome_dir, x))]:
        species_path = os.path.join(pangenome_dir, sp)

        # Check if 'output' directory exists.
        if "output" not in os.listdir(species_path):
            continue

        output_path = os.path.join(species_path, "output")

        # Find the 'gene_presence_absence.csv' file from Roary output.
        roary_csv_path = os.path.join(output_path, "gene_presence_absence.csv")
        if not os.path.exists(roary_csv_path):
            tqdm.write(f"'gene_presence_absence.csv' not found in '{output_path}'. Skipping species '{sp}'.")
            continue

        # Read the gene_presence_absence.csv file.
        df_pangenome = pd.read_csv(roary_csv_path)

        # Process the DataFrame.
        df_pangenome = df_pangenome.set_index("Gene")
        df_pangenome = df_pangenome.drop(df_pangenome.columns[:13], axis=1)
        # Replace NaN values with 0 and any other value with 1.
        df_pangenome = df_pangenome.fillna(0).applymap(lambda x: 0 if x == 0 else 1)
        # Transpose the DataFrame so that samples are rows and genes are columns.
        df_pangenome = df_pangenome.T
        # Sort columns by sum.
        df_pangenome = df_pangenome[df_pangenome.sum().sort_values(ascending=False).index]
        # Filter out rows with sum zero.
        df_pangenome = df_pangenome.loc[(df_pangenome.sum(axis=1) != 0), :]

        if df_pangenome.empty:
            tqdm.write(f"No data available for species '{sp}'. Skipping.")
            continue

        # Define a custom colorscale from white to very pastel blue.
        custom_colorscale = [
            [0, 'white'],  # Start at white.
            [1, '#9aafc7']  # End at very pastel blue.
        ]

        # Create the heatmap.
        heatmap_pangenomic = go.Figure(
            data=[
                go.Heatmap(
                    x=df_pangenome.columns,
                    y=df_pangenome.index,
                    z=df_pangenome.values,
                    colorscale=custom_colorscale,
                    hoverongaps=False,
                    hovertemplate=
                        'Value: %{z}<br>' +
                        'Gene: %{x}<br>' +
                        'Sample: %{y}<br>' +
                        '<extra></extra>'
                )
            ],
            layout=go.Layout(
                title=f"Gene Presence and Absence - {species_dict.get(sp, sp)}",
                xaxis_title="Genes",
                yaxis_title="Samples",
                showlegend=False
            )
        )

        # Update the layout of the heatmap.
        heatmap_pangenomic.update_layout(
            height=max(len(df_pangenome.index) * 20, 200),
            title_font_size=20,
            xaxis=dict(title_font_size=18),
            yaxis=dict(title_font_size=18)
        )

        # Fix colorbar (assuming fix_colorbar is defined).
        heatmap_pangenomic = fix_colorbar(heatmap_pangenomic, facet=False)

        figures[sp] = heatmap_pangenomic

    return figures


def create_pangenome_pie_charts(species_dict):

    figures = {}

    for sp in [x for x in os.listdir(f"pangenome") if os.path.isdir(f"pangenome/{x}")]:

        if "output" not in os.listdir(f"pangenome/{sp}"):
            continue

        with open(f'pangenome/{sp}/output/summary_statistics.txt', 'r') as f:
            lines = f.readlines()

        # Create a dictionary with the first word as the key and the last word as the value.
        data_dict = {}
        for line in lines:
            line = line.strip().split()
            data_dict[line[0]] = line[-1]

        # Create a DataFrame from the dictionary.
        df_pangenome = pd.DataFrame.from_dict(data_dict, orient='index', columns=['Value'])
        df_pangenome = df_pangenome.drop('Total', axis=0)

        # Create a Pie chart trace with Set3 colors.
        trace = go.Pie(
            labels=df_pangenome.index,
            values=df_pangenome['Value'],
            marker=dict(colors=qualitative.Set2)  # Apply Set3 color palette.
        )

        # Create a layout with a title and set height.
        layout = go.Layout(
            height=400,
            title_font_size=20
        )

        # Create a Figure object and add the trace and layout.
        pangenome_pie_chart = go.Figure(data=[trace], layout=layout)
        configure_pie_chart(pangenome_pie_chart)
        figures[sp] = pangenome_pie_chart

    return figures
