import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from tqdm import tqdm


def create_subtype_polar_charts(df, species_dict):

    figures = {}

    # Define a color palette for resistance types.
    color_palette = px.colors.qualitative.Set3
    unique_resistances = sorted(df['RESISTANCE'].dropna().unique())
    resistance_to_color = {resistance: color_palette[i % len(color_palette)] for i, resistance in enumerate(unique_resistances)}
    resistance_to_color['000'] = 'rgba(255,255,255,0)'  # Transparent color for empty space.

    for sp in df["SPECIES"].unique().tolist():

        tqdm.write(f"[pegas] Creating polar chart for species: {sp}")

        # Filter for the specific species and drop rows with NaN in resistance.
        df_species = df[df["SPECIES"] == sp].dropna(subset=["RESISTANCE"])

        # Prepare the data.
        counts = df_species.groupby(['SUBTYPE', 'RESISTANCE']).size().reset_index(name='gene_count')

        # Calculate the average number of genes per sample for each subtype.
        avg_genes_per_subtype = {subtype: len(df_species[df_species['SUBTYPE'] == subtype]) / len(df_species[df_species['SUBTYPE'] == subtype]['SAMPLE'].unique()) for subtype in df_species['SUBTYPE'].unique()}

        # Sort subtypes by average number of genes per sample.
        subtypes = sorted(avg_genes_per_subtype, key=avg_genes_per_subtype.get, reverse=True)

        resistances = sorted(df_species['RESISTANCE'].unique())

        if len(resistances) == 0:
            tqdm.write(f"[pegas] No resistance data for species '{sp}'. Skipping.")
            continue

        # Count the number of samples in each subtype.
        subtype_counts = df_species.groupby('SUBTYPE')['SAMPLE'].nunique()
        subtype_counts_dict = subtype_counts.to_dict()

        # Calculate maximum gene count to scale the empty section.
        max_gene_count = max(avg_genes_per_subtype.values())
        empty_space_count = max_gene_count * 2 / 8  # Represents 30% of the max gene count average.

        # Add the new "000" resistance type with 30% of the max gene count.
        empty_rows = pd.DataFrame([
            {'SUBTYPE': subtype, 'RESISTANCE': '000', 'gene_count': empty_space_count * subtype_counts_dict[subtype]}
            for subtype in subtypes
        ])
        counts = pd.concat([counts, empty_rows], ignore_index=True)

        # Update the subtypes and resistances lists.
        resistances = sorted(counts['RESISTANCE'].unique())
        bar_width = min(360 / len(subtypes), 10)  # Degrees.
        subtype_to_angle = {subtype: index * bar_width for index, subtype in enumerate(subtypes)}

        # Create the windrose chart.
        data = []

        for resistance in resistances:
            df_res = counts[counts['RESISTANCE'] == resistance]
            df_res_full = pd.DataFrame({'SUBTYPE': subtypes})
            df_res_full = df_res_full.merge(df_res[['SUBTYPE', 'gene_count']], on='SUBTYPE', how='left')
            df_res_full['gene_count'] = df_res_full['gene_count'].fillna(0)
            df_res_full['theta'] = df_res_full['SUBTYPE'].map(subtype_to_angle)
            df_res_full['sample_count'] = df_res_full['SUBTYPE'].map(subtype_counts_dict)
            df_res_full['value'] = df_res_full['gene_count'] / [subtype_counts_dict[x] for x in df_res_full['SUBTYPE']]
            color = resistance_to_color[resistance]

            trace = go.Barpolar(
                r=df_res_full['value'] * 10,
                theta=df_res_full['theta'],
                width=[bar_width] * len(df_res_full),
                name=resistance,
                marker_color=color,
                text=[
                    f"Subtype: {subtype}<br>Resistance: {resistance}<br>Genes: {count}"
                    for subtype, count, value in zip(df_res_full['SUBTYPE'], df_res_full['gene_count'], df_res_full['value'])
                ],
                hoverinfo='text'
            )
            data.append(trace)

        # Calculate cumulative gene counts for positioning labels.
        cumulative_counts = counts.groupby(['SUBTYPE', 'RESISTANCE'])['gene_count'].sum().unstack(fill_value=0)
        cumulative_counts['total'] = cumulative_counts.sum(axis=1)
        max_r_values = cumulative_counts['total']

        max_r_value = max(max_r_values)

        # Prepare the data for the labels.
        label_angles = [subtype_to_angle[subtype] for subtype in subtypes]
        label_radii = [max_r_values[subtype] / subtype_counts_dict[subtype] * 10 + np.log(max_r_value) * 1.2 + 1 for subtype in subtypes]  # Slightly beyond the bar.

        # Add a Scatterpolar trace for the subtype labels.
        label_trace = go.Scatterpolar(
            r=label_radii,
            theta=label_angles,
            mode='text',
            text=[f"<b>ST{str(subtype)}</b>" for subtype in subtypes],
            textfont=dict(
                size=14, color='grey',
                weight='bold'
            ),
            textposition='middle center',

            hoverinfo='none',
            showlegend=False
        )
        data.append(label_trace)

        # Calculate the number of samples per subtype.
        samples_per_subtype = df_species.groupby('SUBTYPE')['SAMPLE'].nunique()

        # Prepare data for the circular line plot.
        line_radii = [samples_per_subtype.get(subtype, 0) for subtype in subtypes]
        line_angles = [subtype_to_angle[subtype] for subtype in subtypes]

        max_avg_genes = max(avg_genes_per_subtype.values())
        max_sample_count = max(samples_per_subtype.values) if len(samples_per_subtype.values) > 0 else 0
        max_bar_radius = (max_avg_genes + empty_space_count) * 10

        # Normalize sample counts so the line stays inside the wedge stack.
        line_min = max_bar_radius * 0.35
        line_max = max_bar_radius * 0.75
        if max_sample_count > 0:
            scaled_line_radii = [
                line_min + (count / max_sample_count) * (line_max - line_min)
                for count in line_radii
            ]
        else:
            scaled_line_radii = [line_min] * len(line_radii)

        remaining_line_angles = [angle for angle in range(0, 360, int(bar_width)) if angle not in line_angles]
        scaled_line_radii += [max_avg_genes * 1.8] * len(remaining_line_angles)
        line_angles += remaining_line_angles

        # Close the circle by adding the first point at the end.
        scaled_line_radii[-1] = scaled_line_radii[0]
        line_angles[-1] = line_angles[0]

        line_trace = go.Scatterpolar(
            r=scaled_line_radii,
            theta=line_angles,
            mode='lines+markers',
            line=dict(color='#B3D9EE', width=2),
            marker=dict(size=2, color='#B3D9EE'),
            # fill='toself',  # This fills the area under the line.
            # fillcolor='#f5e7e4',  # Set the fill color and opacity.
            name='Sample Count',
            hoverinfo='text',
            line_shape='spline',
            text=[f"Subtype: {subtype}<br>Samples: {samples_per_subtype.get(subtype, 0)}" for subtype in subtypes] + [f"Subtype: {subtypes[0]}<br>Samples: {samples_per_subtype.get(subtypes[0], 0)}"],
            showlegend=False
        )
        data.append(line_trace)

        # Add the species title in the middle with enhanced styling.
        title_trace = go.Scatterpolar(
            r=[0],
            theta=[0],
            mode='text',
            text=[species_dict.get(sp, sp).replace(" ", '<br>')],
            textfont=dict(
                size=24,  # Slightly larger font size for emphasis.
                color="violet",  # Font color.
                family='Arial Narrow Italic, Arial Italic, sans-serif',  # Italic font family.
                weight='bold'  # Bold font weight for emphasis.
            ),
            textposition='middle center',
            hoverinfo='none',
            showlegend=False
        )
        data.append(title_trace)

        # Customize the chart layout.
        fig = go.Figure(data=data)
        fig.update_layout(
            template='plotly_white',
            polar=dict(
                radialaxis=dict(
                    ticks='',
                    showgrid=False,
                    showticklabels=False,
                    showline=False
                ),
                angularaxis=dict(
                    ticks='',
                    rotation=90,
                    direction='clockwise',
                    showticklabels=False,
                    showgrid=False,
                    showline=False
                )
            ),
            bargap=0,
            height=1800,
            width=1800,
            margin=dict(l=0, r=0, t=0, b=0),
            showlegend=False,
            barmode='stack'
        )

        figures[sp] = fig

    return figures


def generate_virulence_factor_bar_chart(df, species_dict):
    """
    Generates bar charts for virulence factors for each species.

    Parameters:
    - df (DataFrame): The DataFrame containing the data.
    - species_dict (dict): Dictionary mapping species codes to species names.

    Returns:
    - virulence_factor_plots (dict): Dictionary with species names as keys and bar chart figures as values.
    """
    virulence_factor_plots = {}

    # Determine a consistent bar height (in pixels).
    bar_height = 20  # Adjust as needed for desired bar thickness.

    # Determine the maximum number of virulence factors across all species.
    max_virulence_factors = 0
    species_data = {}

    species = df['SPECIES'].unique()
    chart_heights = {}

    for sp in species:
        species_name = species_dict.get(sp, sp)
        df_species = df[(df['SPECIES'] == sp) & (df['PREDICTION_SOURCE'] == "VFDB")]

        # Count the appearance of each virulence factor.
        virulence_counts = df_species['GENE'].value_counts()

        if len(virulence_counts) == 0:
            continue

        # Sort the virulence factors by count.
        virulence_counts = virulence_counts.sort_values(ascending=False)

        # Update the maximum number of virulence factors.
        max_virulence_factors = max(max_virulence_factors, len(virulence_counts))

        # Store the data for later use.
        species_data[species_name] = virulence_counts

        chart_heights[species_name] = bar_height * len(df_species['GENE'].unique())

    # Now generate the figures with consistent bar widths.
    for species_name, virulence_counts in species_data.items():
        # Calculate figure height based on the number of bars.
        chart_height = chart_heights[species_name]

        # Create a bar chart.
        fig = go.Figure()
        fig.add_trace(go.Bar(
            y=virulence_counts.index,  # Horizontal bars.
            x=virulence_counts.values,  # Corresponding values.
            orientation='h',  # Horizontal orientation.
            marker_color='rgb(255, 204, 153)',  # Pastel yellow-orange.
            hoverinfo='x',
            text=virulence_counts.values,  # Show count values.
            textposition='outside',  # Place text outside the bar.
        ))

        # Customize layout.
        fig.update_layout(
            title=f"{species_name}",
            title_font_size=30,
            xaxis_title="Counts",
            yaxis_title="Virulence Factors",
            plot_bgcolor="white",  # White plot background.
            paper_bgcolor="white",  # White paper background.
            font=dict(
                family="Arial, sans-serif",
                size=14
            ),
            margin=dict(l=150, r=50, t=80, b=50),  # Adjust margins for readability.
            height=chart_height,  # Set chart height based on max bars.
            yaxis=dict(
                automargin=True,
                tickfont=dict(size=12),
                categoryorder='total ascending'  # Ensure consistent order.
            ),
            xaxis=dict(
                automargin=True,
                tickfont=dict(size=12)
            )
        )

        # Update bar thickness.
        fig.update_traces(
            marker_line_width=1,
            marker_line_color='black',
            selector=dict(type='bar')
            # width=0.8  # Set a fixed bar width relative to the plot.
        )

        # Store the figure in the dictionary.
        virulence_factor_plots[species_name] = fig

    return virulence_factor_plots
