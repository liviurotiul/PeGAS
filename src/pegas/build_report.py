import pandas as pd
import itertools
import re
import networkx as nx
import plotly.graph_objects as go
import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt
import plotly
import matplotlib
import math
import warnings

from plotly.subplots import make_subplots
from itertools import permutations
from bs4 import BeautifulSoup as bs
from collections import Counter
from jinja2 import Template
from jinja2 import Environment, FileSystemLoader
from tqdm import tqdm
from plotly.colors import qualitative

from utils import *
from itertools import chain

import warnings
import pandas as pd
import plotly.graph_objs as go
import plotly.tools
import matplotlib


def create_sunburst_chart(df, species_dict, species_color_mapping):
    """
    Creates a sunburst chart based on sample, subtype, and species data.
    
    Parameters:
    - df (DataFrame): The dataframe containing the data.
    - species_dict (dict): Mapping from species codes to common names.
    - species_color_mapping (dict): Mapping from species common names to colors.
    
    Returns:
    - sunburst_html (str): HTML code for the generated sunburst chart.
    """
    label = []
    parent = []
    value = []
    node_info = []
    colors = []

    # Grouping and preparing data
    df_samples = df.groupby(["SAMPLE", "SUBTYPE", "SPECIES"]).size().reset_index(name="count")
    df_samples["SPECIES_COMMON_NAME"] = df_samples["SPECIES"].map(species_dict)
    sample_list = df_samples["SAMPLE"].unique().tolist()

    for _, sample in df_samples.iterrows():
        # Use SPECIES_COMMON_NAME for species
        subtype_species = f"ST{sample['SUBTYPE']}({sample['SPECIES_COMMON_NAME']})"
        
        # Append the sample label
        label.append(sample["SAMPLE"])
        parent.append(subtype_species)
        value.append(1)
        
        # Check if the subtype-species combination has been added as a label
        if subtype_species not in label:
            label.append(subtype_species)
            parent.append(sample["SPECIES_COMMON_NAME"])
            value.append(0)
        
        # Check if the species has been added as a label
        if sample["SPECIES_COMMON_NAME"] not in label:
            label.append(sample["SPECIES_COMMON_NAME"])
            parent.append("")
            value.append(0)

    for par, item in zip(parent, label):
        # Condition for species-level nodes
        if item in df_samples["SPECIES_COMMON_NAME"].unique():
            no_of_samples = df_samples[df_samples["SPECIES_COMMON_NAME"] == item].shape[0]
            node_info.append(
                f"Number of isolates: {no_of_samples} <br> {round(no_of_samples / len(sample_list) * 100, 2)}% of total isolates"
            )
            colors.append(species_color_mapping[item])

        # Condition for subtype-level nodes
        elif item.split("(")[0].replace("ST", "") in df_samples["SUBTYPE"].astype(str).tolist():
            subtype = item.split("(")[0].replace("ST", "")
            species_common_name = item.split("(")[1].replace(")", "")
            no_of_samples_species = df_samples[df_samples["SPECIES_COMMON_NAME"] == species_common_name]
            no_of_samples_subtype = no_of_samples_species[no_of_samples_species["SUBTYPE"] == subtype].shape[0]
            node_info.append(
                f"Number of isolates: {no_of_samples_subtype} <br> {round(no_of_samples_subtype / len(sample_list) * 100, 2)}% of total isolates"
            )
            diluted_color = dilute_hex_color(species_color_mapping[species_common_name], 0.05)
            colors.append(diluted_color)

        # Condition for sample-level nodes
        elif item in df_samples["SAMPLE"].tolist():
            node_info.append("")
            species_common_name = par.split("(")[1].replace(")", "")
            diluted_color = dilute_hex_color(species_color_mapping[species_common_name], 0.1)
            colors.append(diluted_color)

    sunburst_data = {
        "label": label,
        'parent': parent,
        'value': value
    }

    # Define the trace for the sunburst chart
    sunburst_trace = go.Sunburst(
        labels=sunburst_data['label'],
        parents=sunburst_data['parent'],
        values=sunburst_data['value'],
        customdata=node_info,
        hovertemplate='<b>%{label}</b><br>%{customdata}<extra></extra>',
        marker=dict(colors=colors)
    )

    sunburst_figure = go.Figure(data=sunburst_trace)
    sunburst_figure.update_layout(
        font=dict(
            family="Arial, sans-serif",  # Use a professional sans-serif font like Arial
            size=14,                      # Adjust font size as needed
            color="black"                 # Set font color
        )
    )
    sunburst_figure.update_layout(height=1500)
    sunburst_figure.update_layout(
        hoverlabel_namelength=0,
        title="Sequence typing: sunburst chart"
    )

    # Generate HTML code for the sunburst chart (assuming create_html_element is defined)
    sunburst_html = create_html_element(sunburst_figure, "")

    return sunburst_figure

def load_and_preprocess_data():
    """
    Loads the results.csv file and preprocesses the dataframe.
    
    Returns:
    - df (DataFrame): Preprocessed dataframe.
    """
    warnings.filterwarnings("ignore")
    matplotlib.use('Agg')

    # Read results.csv
    df = pd.read_csv("dataframe/results.csv", dtype={'SAMPLE': str})

    # Set the columns to the correct data types
    df = df.astype({
        'SAMPLE': str,
        'SPECIES': str,
        'SUBTYPE': str,
        'GENE': str,
        'ACCESSION': str,
        'PRODUCT': str,
        '%IDENTITY': float,
        '%COVERAGE': float,
        'PREDICTION_SOURCE': str,
        'CONTIG_LENGTH': int,
        'CONTIG_COVERAGE': float,
        'CONTIG_NUMBER': int
    })

    # Assign colors to species (assuming species_dict is defined)
    species_color_mapping = assign_colors(species_dict)

    df["SEQUENCE"] = df["SEQUENCE"].str.replace("contig", "")

    # Convert to int
    df["CONTIG"] = df["SEQUENCE"].astype(int)

    return df.drop(columns=["SEQUENCE"])

def create_heatmap_figure(df, prediction_source, title):
    """
    Creates identity and coverage heatmaps for the specified prediction source.

    Parameters:
    - df (DataFrame): The dataframe containing the data.
    - prediction_source (str): The prediction source to filter by (e.g., 'VFDB').
    - title (str): The title for the heatmap.

    Returns:
    - heatmap_html (str): HTML code for the generated heatmap.
    """
    # Filter the dataframe based on the prediction source
    df_filtered = df[df["PREDICTION_SOURCE"] == prediction_source]

    if df_filtered.empty:
        print(f"No data found for prediction source '{prediction_source}'.")
        return ""

    # Add an ID column
    df_filtered = df_filtered.copy()
    df_filtered.reset_index(inplace=True)
    df_filtered.rename(columns={'index': 'ID'}, inplace=True)

    # Create a mapping from ID to Gene
    id_to_gene = dict(zip(df_filtered['ID'], df_filtered['GENE']))

    # Pivot the dataframe for identity and coverage
    df_identity = df_filtered.pivot_table(index='SAMPLE', columns='ID', values='%IDENTITY', aggfunc='max').fillna(0)
    df_coverage = df_filtered.pivot_table(index='SAMPLE', columns='ID', values='%COVERAGE', aggfunc='max').fillna(0)

    colorscale = "Blues"

    # Create Heatmap for Identity
    data_identity = go.Heatmap(
        z=df_identity.values,
        x=[id_to_gene[item] for item in df_identity.columns],
        y=df_identity.index,
        colorscale=colorscale,
        hoverongaps=False,
        hovertemplate=(
            'Isolate: %{y}<br>' +
            'Gene: %{x}<br>' +
            'Identity: %{z}%<br><extra></extra>'
        )
    )

    # Create Heatmap for Coverage
    data_coverage = go.Heatmap(
        z=df_coverage.values,
        x=[id_to_gene[item] for item in df_coverage.columns],
        y=df_coverage.index,
        colorscale=colorscale,
        hoverongaps=False,
        hovertemplate=(
            'Isolate: %{y}<br>' +
            'Gene: %{x}<br>' +
            'Coverage: %{z}%<br><extra></extra>'
        )
    )

    # Combine the heatmaps into a subplot
    heatmap_figure = plotly.tools.make_subplots(rows=2, cols=1, shared_xaxes=True)
    heatmap_figure.append_trace(data_identity, 1, 1)
    heatmap_figure.append_trace(data_coverage, 2, 1)

    heatmap_figure.update_layout(
        title={
            'text': title,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 20}
        },
        height=max(300 + len(df_identity.index) * 10, 600)
    )

    # Fix the colorbar (assuming fix_colorbar is defined)
    heatmap_figure = fix_colorbar(heatmap_figure, facet=True, dtick=10)

    # Generate HTML code for the heatmap (assuming create_html_element is defined)
    heatmap_html = create_html_element(heatmap_figure, f"{prediction_source} output - identity & coverage heatmap")

    return heatmap_figure

def create_annotation_tables(df, species_dict):
    """
    Creates three HTML tables summarizing annotations from NCBI, VFDB, and PlasmidFinder,
    with custom CSS to set the table width to 100%.

    Parameters:
    - df (DataFrame): The dataframe containing the data.
    - species_dict (dict): Mapping from species codes to common names.

    Returns:
    - tables (dict): A dictionary with keys 'Resistance', 'Virulence', 'Plasmid' and values as HTML code for the tables.
    """
    from bs4 import BeautifulSoup

    # Map species codes to common names in the DataFrame
    df['SPECIES_NAME'] = df['SPECIES'].map(species_dict)

    # Define predictors and their corresponding table names
    predictors = {
        'NCBI': 'Resistance',
        'VFDB': 'Virulence',
        'PlasmidFinder': 'Plasmid'
    }

    # Initialize a dictionary to hold the HTML tables
    tables = {}

    for predictor, table_name in predictors.items():
        # Filter the DataFrame for the current predictor
        df_predictor = df[df['PREDICTION_SOURCE'] == predictor]

        # Check if the filtered DataFrame is not empty
        if not df_predictor.empty:
            # Select relevant columns (exclude unnecessary ones)
            df_table = df_predictor.drop(columns=["CONTIG_LENGTH", "CONTIG_COVERAGE", "CONTIG_NUMBER"], errors='ignore')

            # Remove duplicates if any
            df_table = df_table.drop_duplicates()

            # Sort the table by SAMPLE for better readability
            df_table = df_table.sort_values(by=['SAMPLE'])

            # Generate HTML code for the table
            table_html = df_table.to_html(
                classes='table table-striped', # Might have to comment this; it s showing search and pagination 2 times
                index=False,
                border=0,
                table_id=f'{table_name}_table'
            )

            # Parse the HTML to modify the table tag
            soup = BeautifulSoup(table_html, 'html.parser')

            # Set the style attribute to make the table width 100%
            table_tag = soup.find('table')
            if table_tag:
                table_tag['style'] = 'width: 100%;'

            # Convert back to string
            table_html = str(soup)

            # Store the HTML table in the dictionary
            tables[table_name] = create_html_widget(table_html, f"{table_name}")
        else:
            # If there is no data for this predictor, return a message
            tables[table_name] = f'<p>No data available for {table_name}</p>'

    return tables

def create_fastqc_table(df, species_dict, gc_content_dict):
    """
    Creates an HTML table summarizing FastQC outputs for all samples.

    Parameters:
    - df_samples (DataFrame): DataFrame containing sample information.
    - species_dict (dict): Mapping from species codes to common names.
    - gc_content_dict (dict): Mapping from species common names to expected GC content.

    Returns:
    - fastqc_table_html (str): HTML code for the generated FastQC table.
    """

    df_samples = df.groupby("SAMPLE").agg({
        "CONTIG_LENGTH": "first",
        "CONTIG_COVERAGE": "first",
        "SPECIES": "first",
        "CONTIG_NUMBER": "first"
    }).reset_index()

    # Find all FastQC HTML files
    html_files = []
    for root, dirs, files in os.walk("fastqc"):
        for file in files:
            if file.endswith(".html"):
                html_files.append(os.path.join(root, file))

    # Define the strings we're looking for in the FastQC reports
    search_strings = ["Total Sequences", "Sequences flagged as poor quality", "Sequence length", "%GC"]

    # Function to extract sample name
    def get_core_sample_name(filename):
        """Extracts the core sample name by removing _R1 or _R2 and other suffixes."""
        return os.path.basename(filename).replace("_R1", "").replace("_R2", "").replace(".fastq.gz", "").replace("_fastqc.html", "")

    # Function to extract read direction from file path
    def get_read_direction_from_file(file_path):
        base_name = os.path.basename(file_path)
        match = re.search(r'_R([12])', base_name)
        if match:
            return 'R' + match.group(1)
        else:
            # Try to infer from the directory structure if not in filename
            if 'R1' in file_path or 'read1' in file_path.lower():
                return 'R1'
            elif 'R2' in file_path or 'read2' in file_path.lower():
                return 'R2'
            else:
                return None

    # Build a set of unique sample names
    sample_names = set()
    for file_path in html_files:
        sample_name = get_core_sample_name(file_path)
        sample_names.add(sample_name)

    # Create the DataFrame with the corrected index
    df_qc = pd.DataFrame(
        index=sorted(sample_names),
        columns=[f"{x} R1" for x in search_strings] + [f"{x} R2" for x in search_strings] + ['SPECIES']
    )

    # Loop over each file path and extract data
    for file_path in html_files:
        # Load the HTML file into BeautifulSoup
        with open(file_path, "r", encoding='utf-8') as f:
            soup = bs(f, "html.parser")

        sample_name = get_core_sample_name(file_path)
        read_direction = get_read_direction_from_file(file_path)

        if read_direction is None:
            tqdm.write(f"Warning: Could not determine read direction for file '{file_path}'. Skipping.")
            continue  # Skip files without read direction

        # Loop over each search string
        for search_string in search_strings:
            # Find the <td> element with the search string as its content
            td = soup.find("td", string=search_string)

            # If the <td> element exists, extract the content of the next <td> element
            if td is not None:
                next_td = td.find_next_sibling("td")
                if next_td is not None:
                    column_name = f"{search_string} {read_direction}"
                    df_qc.loc[sample_name, column_name] = next_td.get_text().strip()
                else:
                    tqdm.write(f"Warning: Next <td> element not found for '{search_string}' in file '{file_path}'.")
            else:
                tqdm.write(f"Warning: '{search_string}' not found in file '{file_path}'.")

    # Assign 'SPECIES' to each sample in df_qc based on df_samples and species_dict
    for sample_name in df_qc.index:
        sample_row = df_samples[df_samples["SAMPLE"] == sample_name]

        if not sample_row.empty:
            species_code = sample_row["SPECIES"].values[0]
            species_common_name = species_dict.get(species_code, 'Unknown')
            df_qc.loc[sample_name, "SPECIES"] = species_common_name
        else:
            tqdm.write(f"Warning: Sample '{sample_name}' not found in df_samples. Assigning 'Unknown' to SPECIES.")
            df_qc.loc[sample_name, "SPECIES"] = 'Unknown'

    # Convert relevant columns to numeric, handling non-numeric values
    numeric_columns = [col for col in df_qc.columns if any(x in col for x in ["Total Sequences", "Sequences flagged as poor quality"])]
    for col in numeric_columns:
        df_qc[col] = pd.to_numeric(df_qc[col], errors='coerce')

    # Calculate percentages and format the 'Sequences flagged as poor quality' columns
    for read_direction in ['R1', 'R2']:
        total_sequences_col = f"Total Sequences {read_direction}"
        poor_quality_col = f"Sequences flagged as poor quality {read_direction}"
        if total_sequences_col in df_qc.columns and poor_quality_col in df_qc.columns:
            # Avoid division by zero and handle missing data
            total_sequences = df_qc[total_sequences_col].replace(0, pd.NA)
            poor_quality_sequences = df_qc[poor_quality_col]
            with pd.option_context('mode.use_inf_as_na', True):
                percentage = (poor_quality_sequences / total_sequences * 100).round(2)
            # Combine counts and percentages
            df_qc[poor_quality_col] = poor_quality_sequences.fillna(0).astype(int).astype(str) + ' (' + percentage.fillna(0).astype(str) + '%)'

    # Compute %GC intervals and annotate with expected values
    for read_direction in ['R1', 'R2']:
        gc_col = f"%GC {read_direction}"
        if gc_col in df_qc.columns:
            def compute_GC_intervals(row):
                sp = row["SPECIES"]
                actual_gc = row[gc_col]
                if sp not in gc_content_dict:
                    return f"{actual_gc}% (expected: N/A)"
                expected_gc = gc_content_dict[sp]
                return f"{actual_gc}% (expected: {expected_gc}%)"
            df_qc[gc_col] = df_qc.apply(compute_GC_intervals, axis=1)

    # Generate HTML table code
    df_reads_table_code = df_qc.to_html(classes='table table-striped', table_id='full_qa_table', index=True)

    # Parse the HTML code using BeautifulSoup
    soup = bs(df_reads_table_code, 'html.parser')

    # Style the table headers
    for th in soup.find_all('th'):
        th['style'] = 'text-align: left'

    # Apply color coding to cells based on conditions
    # Define a function to apply colors based on the difference in %GC
    def apply_gc_color(td):
        text = td.get_text().strip()
        if "N/A" in text:
            td['style'] = 'background-color: #f8d7da;'  # Light red
            return

        match = re.search(r'([\d\.]+)% \(expected: ([\d\.]+)%\)', text)
        if match:
            actual_gc = float(match.group(1))
            expected_gc = float(match.group(2))
            difference = abs(actual_gc - expected_gc)
            if difference > 6:
                color = '#f8d7da'  # Light red
            else:
                color = '#d4edda'  # Light green
            td['style'] = f'background-color: {color};'
        else:
            td['style'] = 'background-color: #ffffff;'

    for row in soup.find('table').find_all('tr'):
        for td in row.find_all('td'):
            # Assuming the %GC columns are at specific positions
            for gc_col in [f"%GC R1", f"%GC R2"]:
                if gc_col in df_qc.columns:
                    apply_gc_color(td)
     
    # Convert the modified table back to HTML
    df_reads_table_code = str(soup)

    # Wrap the table in a widget
    fastqc_table_html = create_html_widget(df_reads_table_code, "FastQC output - table")

    return fastqc_table_html

def create_subtype_polar_charts(df, species_dict):
    
    figures = {}

    # Define a color palette for resistance types
    color_palette = px.colors.qualitative.Set3
    unique_resistances = sorted(df['RESISTANCE'].dropna().unique())
    resistance_to_color = {resistance: color_palette[i % len(color_palette)] for i, resistance in enumerate(unique_resistances)}
    resistance_to_color['000'] = 'rgba(255,255,255,0)'  # Transparent color for empty space

    for sp in df["SPECIES"].unique().tolist():

        # Filter for the specific species and drop rows with NaN in resistance
        df_species = df[df["SPECIES"] == sp].dropna(subset=["RESISTANCE"])

        # Prepare the data
        counts = df_species.groupby(['SUBTYPE', 'RESISTANCE']).size().reset_index(name='gene_count')

        # Calculate the average number of genes per sample for each subtype
        avg_genes_per_subtype = {subtype: len(df_species[df_species['SUBTYPE'] == subtype]) / len(df_species[df_species['SUBTYPE'] == subtype]['SAMPLE'].unique()) for subtype in df_species['SUBTYPE'].unique()}

        # Sort subtypes by average number of genes per sample
        subtypes = sorted(avg_genes_per_subtype, key=avg_genes_per_subtype.get, reverse=True)

        resistances = sorted(df_species['RESISTANCE'].unique())

        if len(resistances) == 0 or len(subtypes) < 2:
            continue

        # Count the number of samples in each subtype
        subtype_counts = df_species.groupby('SUBTYPE')['SAMPLE'].nunique()
        subtype_counts_dict = subtype_counts.to_dict()

        # Calculate maximum gene count to scale the empty section
        total_gene_counts = counts.groupby('SUBTYPE')['gene_count'].sum()
        max_gene_count = max(avg_genes_per_subtype.values())
        empty_space_count = max_gene_count * 2 / 8 # Represents 30% of the max gene count average
    
        # Add the new "000" resistance type with 30% of the max gene count
        empty_rows = pd.DataFrame([
            {'SUBTYPE': subtype, 'RESISTANCE': '000', 'gene_count': empty_space_count*subtype_counts_dict[subtype]}
            for subtype in subtypes
        ])
        counts = pd.concat([counts, empty_rows], ignore_index=True)

        # Update the subtypes and resistances lists
        resistances = sorted(counts['RESISTANCE'].unique())
        num_subtypes = len(subtypes)
        bar_width = min(360 / num_subtypes, 10)  # Degrees
        subtype_to_angle = {subtype: index * bar_width for index, subtype in enumerate(subtypes)}

        # Create the windrose chart
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
                r=df_res_full['value']*10,
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

        # Calculate cumulative gene counts for positioning labels
        cumulative_counts = counts.groupby(['SUBTYPE', 'RESISTANCE'])['gene_count'].sum().unstack(fill_value=0)
        cumulative_counts['total'] = cumulative_counts.sum(axis=1)
        max_r_values = cumulative_counts['total']

        max_r_value = max(max_r_values)

        # Prepare the data for the labels
        label_angles = [subtype_to_angle[subtype] for subtype in subtypes]
        label_radii = [max_r_values[subtype] / subtype_counts_dict[subtype] * 10 + np.log(max_r_value)*1.2 + 1 for subtype in subtypes]  # Slightly beyond the bar

        # Add a Scatterpolar trace for the subtype labels
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

        # Calculate the number of samples per subtype
        samples_per_subtype = df_species.groupby('SUBTYPE')['SAMPLE'].nunique()

        # Prepare data for the circular line plot
        line_radii = [samples_per_subtype.get(subtype, 0) for subtype in subtypes]
        line_angles = [subtype_to_angle[subtype] for subtype in subtypes]

        # Calculate scaling factor
        max_sample_count = max(samples_per_subtype.values)
        scaling_factor = 1 / max_r_value

        max_avg_genes = max(avg_genes_per_subtype.values())

        # Scale the sample counts
        scaled_line_radii = [max_avg_genes*1.8 + count for count in line_radii]

        # Count the number of polar bars is still space for in the plot
        num_bars = (360 - len(resistances) * bar_width) / bar_width

        remaining_line_angles = [angle for angle in range(0, 360, int(bar_width)) if angle not in line_angles]
        scaled_line_radii += [max_avg_genes*1.8] * len(remaining_line_angles)    
        line_angles += remaining_line_angles

        # Close the circle by adding the first point at the end
        scaled_line_radii[-1] = scaled_line_radii[0]
        line_angles[-1] = line_angles[0]

        line_trace = go.Scatterpolar(
            r=scaled_line_radii,
            theta=line_angles,
            mode='lines+markers',
            line=dict(color='#B3D9EE', width=2),
            marker=dict(size=2, color='#B3D9EE'),
            # fill='toself',  # This fills the area under the line
            # fillcolor='#f5e7e4',  # Set the fill color and opacity
            name='Sample Count',
            hoverinfo='text',
            line_shape='spline',
            text=[f"Subtype: {subtype}<br>Samples: {samples_per_subtype.get(subtype, 0)}" for subtype in subtypes] + [f"Subtype: {subtypes[0]}<br>Samples: {samples_per_subtype.get(subtypes[0], 0)}"],
            showlegend=False
        )
        data.append(line_trace)

        # Add the species title in the middle with enhanced styling
        title_trace = go.Scatterpolar(
            r=[0],
            theta=[0],
            mode='text',
            text=[species_dict.get(sp, sp).replace(" ", '<br>')],
            textfont=dict(
            size=24,                        # Slightly larger font size for emphasis
            color="violet",                 # Font color
            family='Arial Narrow Italic, Arial Italic, sans-serif',  # Italic font family
            weight='bold'                   # Bold font weight for emphasis
            ),
            textposition='middle center',
            hoverinfo='none',
            showlegend=False
        )
        data.append(title_trace)

        # Customize the chart layout
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

def create_pangenome_heatmaps(species_dict):
    """
    Generates pangenome heatmaps for each species present in the 'pangenome' directory.

    Parameters:
    - species_dict (dict): Mapping from species codes to common names.

    Returns:
    - pangenome_heatmaps_html (str): HTML code containing all the pangenome heatmaps.
    """
    import os
    import pandas as pd
    import plotly.graph_objs as go

    pangenome_heatmaps_html = ""

    # Loop over each species directory in 'pangenome'
    pangenome_dir = "pangenome"
    if not os.path.exists(pangenome_dir):
        tqdm.write(f"Directory '{pangenome_dir}' does not exist. Skipping pangenome heatmaps.")
        return ""

    figures = {}

    for sp in [x for x in os.listdir(pangenome_dir) if os.path.isdir(os.path.join(pangenome_dir, x))]:
        species_path = os.path.join(pangenome_dir, sp)

        # Check if 'output' directory exists
        if "output" not in os.listdir(species_path):
            continue

        output_path = os.path.join(species_path, "output")

        # Find the 'gene_presence_absence.csv' file from Roary output
        roary_csv_path = os.path.join(output_path, "gene_presence_absence.csv")
        if not os.path.exists(roary_csv_path):
            tqdm.write(f"'gene_presence_absence.csv' not found in '{output_path}'. Skipping species '{sp}'.")
            continue

        # Read the gene_presence_absence.csv file
        df_pangenome = pd.read_csv(roary_csv_path)

        # Process the DataFrame
        df_pangenome = df_pangenome.set_index("Gene")
        df_pangenome = df_pangenome.drop(df_pangenome.columns[:13], axis=1)
        # Replace NaN values with 0 and any other value with 1
        df_pangenome = df_pangenome.fillna(0).applymap(lambda x: 0 if x == 0 else 1)
        # Transpose the DataFrame so that samples are rows and genes are columns
        df_pangenome = df_pangenome.T
        # Sort columns by sum
        df_pangenome = df_pangenome[df_pangenome.sum().sort_values(ascending=False).index]
        # Filter out rows with sum zero
        df_pangenome = df_pangenome.loc[(df_pangenome.sum(axis=1) != 0), :]

        if df_pangenome.empty:
            tqdm.write(f"No data available for species '{sp}'. Skipping.")
            continue

        # Define a custom colorscale from white to very pastel blue
        custom_colorscale = [
            [0, 'white'],        # Start at white
            [1, '#9aafc7']       # End at very pastel blue
        ]

        # Create the heatmap
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

        # Update the layout of the heatmap
        heatmap_pangenomic.update_layout(
            height=max(len(df_pangenome.index) * 20, 200),
            title_font_size=20,
            xaxis=dict(title_font_size=18),
            yaxis=dict(title_font_size=18)
        )

        # Fix colorbar (assuming fix_colorbar is defined)
        heatmap_pangenomic = fix_colorbar(heatmap_pangenomic, facet=False)

        figures[sp] = heatmap_pangenomic

        # Append to the cumulative HTML code
        # pangenome_heatmaps_html += heatmap_pangenomic_code

    return figures

def create_pangenome_pie_charts(species_dict):

    pangenome_pie_chart_codes = ""
    figures = {}

    for sp in [x for x in os.listdir(f"pangenome") if os.path.isdir(f"pangenome/{x}")]:

        if "output" not in os.listdir(f"pangenome/{sp}"):
            continue

        with open(f'pangenome/{sp}/output/summary_statistics.txt', 'r') as f:
            lines = f.readlines()

        # Create a dictionary with the first word as the key and the last word as the value
        data_dict = {}
        for line in lines:
            line = line.strip().split()
            data_dict[line[0]] = line[-1]

        # Create a dataframe from the dictionary
        df_pangenome = pd.DataFrame.from_dict(data_dict, orient='index', columns=['Value'])
        df_pangenome = df_pangenome.drop('Total', axis=0)

        # Create a Pie chart trace with Set3 colors
        trace = go.Pie(
            labels=df_pangenome.index,
            values=df_pangenome['Value'],
            marker=dict(colors=qualitative.Set2)  # Apply Set3 color palette
        )

        # Create a layout with a title and set height
        layout = go.Layout(
            height=400,
            title_font_size=20
        )

        # Combine the trace and layout into a figure
        fig = go.Figure(data=[trace], layout=layout)


        # Create a Figure object and add the trace and layout
        pangenome_pie_chart = go.Figure(data=[trace], layout=layout)
        configure_pie_chart(pangenome_pie_chart)
        figures[sp] = pangenome_pie_chart
    
    return figures

def generate_table_html_for_species(sp, df):
    """
    Generates a grouped HTML table for the given species, organized by SUBTYPE, SAMPLE, and RESISTANCE.

    Parameters:
    - sp (str): The species identifier.
    - df (DataFrame): The dataframe containing the data.

    Returns:
    - str: The HTML string of the grouped table.
    """
    # Filter for the specified species and drop rows with missing RESISTANCE
    df_species = df[df['SPECIES'] == sp].dropna(subset=['RESISTANCE'])

    # Define the columns you want to display
    columns_to_display = ['SUBTYPE', 'SAMPLE', 'RESISTANCE', 'GENE', '%IDENTITY']
    df_species = df_species[columns_to_display]

    # Sort the data for consistent grouping
    df_species.sort_values(['SUBTYPE', 'SAMPLE', 'RESISTANCE'], inplace=True)

    # Initialize the HTML table structure with Bootstrap classes
    table_html = """
    <table class='table table-sm' style='font-family: Arial, sans-serif;'>
        <thead class='thead-light'>
            <tr>
                <th>SUBTYPE</th>
                <th>SAMPLE</th>
                <th>RESISTANCE</th>
                <th>GENES</th>
            </tr>
        </thead>
        <tbody>
    """

    # Group by SUBTYPE, SAMPLE, and RESISTANCE
    grouped = df_species.groupby(['SUBTYPE', 'SAMPLE', 'RESISTANCE'])

    # Initialize variables to keep track of rowspan and previous values
    prev_subtype = None
    prev_sample = None
    subtype_rowspan = {}
    sample_rowspan = {}

    # Compute rowspan for SUBTYPE and SAMPLE
    for (subtype, sample), group in df_species.groupby(['SUBTYPE', 'SAMPLE']):
        subtype_rowspan[subtype] = subtype_rowspan.get(subtype, 0) + len(group['RESISTANCE'].unique())
        sample_rowspan[(subtype, sample)] = sample_rowspan.get((subtype, sample), 0) + len(group['RESISTANCE'].unique())

    # Iterate over groups and build table rows
    for (subtype, sample, resistance), group in grouped:
        genes_html = ""
        for _, row in group.iterrows():
            genes_html += f"""
            <span style='display: inline-block; margin-right: 5px;' title='%IDENTITY: {row['%IDENTITY']}%'>
                {row['GENE']}
            </span>
            """

        # Start table row
        table_html += "<tr>"

        # SUBTYPE cell
        if subtype != prev_subtype:
            rowspan = subtype_rowspan[subtype]
            table_html += f"<td rowspan='{rowspan}' class='subtype-cell'>{subtype}</td>"
            prev_subtype = subtype
            prev_sample = None  # Reset prev_sample when subtype changes

        else:
            # If not first occurrence, skip the subtype cell
            pass

        # SAMPLE cell
        if sample != prev_sample:
            rowspan = sample_rowspan[(subtype, sample)]
            table_html += f"<td rowspan='{rowspan}' class='sample-cell'>{sample}</td>"
            prev_sample = sample
        else:
            # If not first occurrence, skip the sample cell
            pass

        # RESISTANCE cell
        table_html += f"<td class='resistance-cell'>{resistance}</td>"

        # GENES cell
        table_html += f"<td class='genes-cell'>{genes_html}</td>"

        # End table row
        table_html += "</tr>"

    # Close the table
    table_html += "</tbody></table>"

    # Include custom CSS
    custom_css = """
    <style>
    .table {
        border-collapse: collapse;
    }
    .table th, .table td {
        vertical-align: middle;
    }
    /* Remove borders from all cells by default */
    .table td, .table th {
        border: none;
    }
    /* Subtype styling */
    .table td.subtype-cell {
        font-weight: bold;
        font-size: 18px; /* Larger font for SUBTYPE */
        border-top: 2px solid #000; /* Line between subtypes */
    }
    /* Sample styling */
    .table td.sample-cell {
        font-weight: bold;
        font-size: 16px; /* Smaller font for SAMPLE */
        border-top: 1px solid #000; /* Line between samples */
    }
    /* Resistance styling */
    .table td.resistance-cell {
        font-weight: bold; /* Resistance types should be bold */
        font-size: 14px; /* Font size same as GENES */
    }
    /* Genes styling */
    .table td.genes-cell {
        font-size: 14px; /* Font size same as RESISTANCE */
    }
    </style>
    """

    return custom_css + table_html

def create_contig_plot(df, species_dict):
    """
    Creates cumulative contig length plots for each species, highlighting N50 segments.

    Parameters:
    - df (DataFrame): The dataframe containing species and sample information.
    - species_dict (dict): A dictionary mapping species codes to species names.

    Returns:
    - species_contig_plots (dict): Dictionary of cumulative contig length plots per species.
    - n50_contig_length_box_plot (Figure): Box plot of N50 contig lengths for all samples, colored by species.
    - contig_coverage_box_plot (Figure): Box plot of average contig coverage for all samples, colored by species.
    """
    import numpy as np

    # Initialize dictionaries and lists to store plots and data
    species_contig_plots = {}
    n50_lengths = []
    average_coverages = []

    species_list = df['SPECIES'].unique()

    # Use pastel colors from Plotly's qualitative color scales
    pastel_colors = qualitative.Pastel1 + qualitative.Pastel2

    for sp in species_list:
        species_name = species_dict.get(sp, sp)
        samples = df[df['SPECIES'] == sp]['SAMPLE'].unique()

        # Initialize a figure for the species
        fig = go.Figure()

        for idx, sample in enumerate(samples):
            # Get the shovill/contigs.fa file
            contigs_file = f"results/{sample}/shovill/contigs.fa"

            # Check if the contigs file exists
            if not os.path.exists(contigs_file):
                tqdm.write(f"Contigs file for sample {sample} not found.")
                continue

            # Read and parse the contig headers to extract lengths and coverages
            contig_lengths = []
            contig_coverages = []
            with open(contigs_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        header = line.strip()
                        # Extract length and coverage
                        header_fields = header.split()
                        len_field = [x for x in header_fields if x.startswith('len=')]
                        cov_field = [x for x in header_fields if x.startswith('cov=')]
                        if len(len_field) > 0 and len(cov_field) > 0:
                            contig_length = int(len_field[0].split('=')[1])
                            contig_coverage = float(cov_field[0].split('=')[1])
                            contig_lengths.append(contig_length)
                            contig_coverages.append(contig_coverage)

            # Check if we have contig data
            if not contig_lengths:
                tqdm.write(f"No contig data found for sample {sample}.")
                continue

            # Sort the contig lengths in descending order
            sorted_contig_lengths = sorted(contig_lengths, reverse=True)

            # Calculate cumulative contig lengths
            cumulative_lengths = np.cumsum(sorted_contig_lengths)
            cumulative_lengths = np.insert(cumulative_lengths, 0, 0)

            # Prepare x-values (contig indices), starting from zero
            x_values = list(range(len(cumulative_lengths)))

            # Calculate N50 and find the index where N50 is reached
            n50_value, n50_index = calculate_n50_and_index(sorted_contig_lengths)

            n50_lengths.append({
                'Sample': sample,
                'N50': n50_value,
                'Species': species_name
            })

            # Calculate average coverage
            n50_coverage = contig_coverages[n50_index]
            average_coverages.append({
                'Sample': sample,
                'N50_Coverage': n50_coverage,
                'Species': species_name
            })

            # Adjust indices when splitting data at N50 index
            x_n50 = x_values[:n50_index + 2]
            y_n50 = cumulative_lengths[:n50_index + 2]
            x_rest = x_values[n50_index + 1:]
            y_rest = cumulative_lengths[n50_index + 1:]

            # Use pastel color for the sample line
            color_idx = idx % len(pastel_colors)
            sample_color = pastel_colors[color_idx]

            # Plot segment up to N50 in red
            fig.add_trace(go.Scatter(
                x=x_n50,
                y=y_n50,
                mode='lines',
                name=f'{sample} (up to N50)',
                line=dict(color='red', width=3),
                hovertemplate=f'Sample: {sample}<br>Contig Index: %{{x}}<br>Cumulative Length: %{{y}} bp<extra></extra>'
            ))

            # Plot the rest of the line in pastel color
            fig.add_trace(go.Scatter(
                x=x_rest,
                y=y_rest,
                mode='lines',
                name=f'{sample}',
                line=dict(color=sample_color, width=3),
                hovertemplate=f'Sample: {sample}<br>Contig Index: %{{x}}<br>Cumulative Length: %{{y}} bp<extra></extra>'
            ))

        # Update figure layout
        fig.update_layout(
            title=f'Cumulative Contig Lengths for Species {species_name}',
            xaxis_title='Contig Index',
            yaxis_title='Cumulative Contig Length (bp)',
            legend_title='Samples',
            template='plotly_white',
            height=600,
            xaxis=dict(range=[0, None]),  # Ensure x-axis starts from 0
            yaxis=dict(range=[0, None])   # Ensure y-axis starts from 0
        )

        # Store the figure in the dictionary
        species_contig_plots[species_name] = fig


    # Create DataFrames for N50 values and average coverages
    n50_df = pd.DataFrame(n50_lengths)
    coverage_df = pd.DataFrame(average_coverages)

    # Use Set3 color palette for box plots
    colors = qualitative.Set2

    # Create box plot with N50 contig length for all samples, colored by species
    n50_box_fig = px.box(
        n50_df,
        x='Species',
        y='N50',
        color='Species',
        points='all',
        title='N50 Contig Length Distribution Across Species',
        template='plotly_white',
        color_discrete_sequence=colors
    )
    n50_box_fig.update_traces(
        jitter=0.3,
        pointpos=-1.8,
        marker_size=6,
        hovertemplate='Sample: %{customdata[0]}<br>Species: %{x}<br>N50: %{y} bp<extra></extra>',
        customdata=n50_df[['Sample']]
    )
    n50_box_fig.update_layout(
        yaxis_title='N50 Contig Length (bp)',
        xaxis_title='Species',
        height=400  # Adjusted height
    )

    # Create box plot with average contig coverage for all samples, colored by species
    coverage_box_fig = px.box(
        coverage_df,
        x='Species',
        y='N50_Coverage',
        color='Species',
        points='all',
        title='N50 Contig Coverage Across Species',
        template='plotly_white',
        color_discrete_sequence=colors
    )
    coverage_box_fig.update_traces(
        jitter=0.3,
        pointpos=-1.8,
        marker_size=6,
        hovertemplate='Sample: %{customdata[0]}<br>Species: %{x}<br>Coverage: %{y:.2f}X<extra></extra>',
        customdata=coverage_df[['Sample']]
    )
    coverage_box_fig.update_layout(
        yaxis_title='N50 Contig Coverage (X)',
        xaxis_title='Species',
        height=400  # Adjusted height
    )

    return species_contig_plots, n50_box_fig, coverage_box_fig

def calculate_n50_and_index(sorted_lengths):
    """
    Calculates the N50 value and the index where N50 is reached from a list of contig lengths sorted in descending order.

    Parameters:
    - sorted_lengths (list): List of contig lengths sorted in descending order.

    Returns:
    - n50 (int): The N50 contig length.
    - index (int): The index at which N50 is reached.
    """
    total_length = sum(sorted_lengths)
    half_total = total_length / 2
    running_total = 0
    for idx, length in enumerate(sorted_lengths):
        running_total += length
        if running_total >= half_total:
            return length, idx
    return 0, 0

def calculate_n50(sorted_lengths):
    """
    Calculates the N50 value from a list of contig lengths sorted in descending order.

    Parameters:
    - sorted_lengths (list): List of contig lengths sorted in descending order.

    Returns:
    - n50 (int): The N50 contig length.
    """
    total_length = sum(sorted_lengths)
    half_total = total_length / 2
    running_total = 0
    for length in sorted_lengths:
        running_total += length
        if running_total >= half_total:
            return length
    return 0

def generate_general_info_table(df, input_dir, output_dir):
    """
    Generates an HTML table containing general information about the run.

    Parameters:
    - df (DataFrame): The dataframe containing the data.
    - input_dir (str): The path of the input directory.
    - output_dir (str): The path of the output directory.

    Returns:
    - str: The HTML string of the table.
    """
    import pandas as pd
    from datetime import datetime

    # Calculate the required values
    report_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    num_samples = len(df['SAMPLE'].unique())
    num_genes = len(df['GENE'].unique())

    # Build the table data
    table_data = [
        ['Date', report_date],
        ['Number of Samples', num_samples],
        ['Number of Unique Genes', num_genes],
        ['Input Directory', input_dir],
        ['Output Directory', output_dir]
    ]

    # Initialize the HTML table
    html = '''
    <table style="margin: 0 auto; font-family: Arial, sans-serif; border-collapse: collapse;">
    '''

    # Build the table rows
    for label, value in table_data:
        html += f'''
        <tr>
            <td style="font-weight: bold; text-align: right; padding: 5px;">{label}</td>
            <td style="text-align: left; padding: 5px;">{value}</td>
        </tr>
        '''

    html += '''
    </table>
    '''

    return html

def generate_virulence_factor_bar_chart(df, species_dict):
    """
    Generates bar charts for virulence factors for each species.
    
    Parameters:
    - df (DataFrame): The dataframe containing the data.
    - species_dict (dict): Dictionary mapping species codes to species names.
    
    Returns:
    - virulence_factor_plots (dict): Dictionary with species names as keys and bar chart figures as values.
    """
    virulence_factor_plots = {}
    
    # Determine a consistent bar height (in pixels)
    bar_height = 20  # Adjust as needed for desired bar thickness
    
    # Determine the maximum number of virulence factors across all species
    max_virulence_factors = 0
    species_data = {}
    
    species = df['SPECIES'].unique()
    chart_heights = {}

    for sp in species:
        species_name = species_dict.get(sp, sp)
        df_species = df[(df['SPECIES'] == sp) & (df['PREDICTION_SOURCE'] == "VFDB")]
        
        # Count the appearance of each virulence factor
        virulence_counts = df_species['GENE'].value_counts()
        
        # Sort the virulence factors by count
        virulence_counts = virulence_counts.sort_values(ascending=False)
        
        # Update the maximum number of virulence factors
        max_virulence_factors = max(max_virulence_factors, len(virulence_counts))
        
        # Store the data for later use
        species_data[species_name] = virulence_counts

        chart_heights[species_name] = bar_height * len(df_species['GENE'].unique())

    # Now generate the figures with consistent bar widths
    for species_name, virulence_counts in species_data.items():
        # Calculate figure height based on the number of bars
        chart_height = chart_heights[species_name]

        # Create a bar chart
        fig = go.Figure()
        fig.add_trace(go.Bar(
            y=virulence_counts.index,  # Horizontal bars
            x=virulence_counts.values,  # Corresponding values
            orientation='h',            # Horizontal orientation
            marker_color='rgb(255, 204, 153)',  # Pastel yellow-orange
            hoverinfo='x',
            text=virulence_counts.values,       # Show count values
            textposition='outside',             # Place text outside the bar
        ))
        
        # Customize layout
        fig.update_layout(
            title=f"{species_name}",
            title_font_size=30,
            xaxis_title="Counts",
            yaxis_title="Virulence Factors",
            plot_bgcolor="white",        # White plot background
            paper_bgcolor="white",       # White paper background
            font=dict(
                family="Arial, sans-serif",
                size=14
            ),
            margin=dict(l=150, r=50, t=80, b=50),  # Adjust margins for readability
            height=chart_height,         # Set chart height based on max bars
            yaxis=dict(
                automargin=True,
                tickfont=dict(size=12),
                categoryorder='total ascending'  # Ensure consistent order
            ),
            xaxis=dict(
                automargin=True,
                tickfont=dict(size=12)
            )
        )
        
        # Update bar thickness
        fig.update_traces(
            marker_line_width=1,
            marker_line_color='black',
            selector=dict(type='bar')
            # width=0.8  # Set a fixed bar width relative to the plot
        )
        
        # Store the figure in the dictionary
        virulence_factor_plots[species_name] = fig
    
    return virulence_factor_plots

def generate_virulence_factor_table(df, species_dict):
    """
    Generates HTML tables of virulence factors for each species.

    Parameters:
    - df (DataFrame): The dataframe containing the data.
    - species_dict (dict): Dictionary mapping species codes to species names.

    Returns:
    - tables (dict): Dictionary with species names as keys and HTML code for tables as values.
    """
    # Filter the DataFrame for virulence factors
    df_virulence = df[df['PREDICTION_SOURCE'] == 'VFDB']

    species = df_virulence['SPECIES'].unique()

    tables = {}

    for sp in species:
        species_name = species_dict.get(sp, sp)

        df_species = df_virulence[df_virulence['SPECIES'] == sp]
        df_species = df_species.drop(columns=['SPECIES', 'CONTIG_LENGTH', 'CONTIG_COVERAGE', 'PREDICTION_SOURCE', 'RESISTANCE', 'CONTIG_NUMBER'])
        df_species.reset_index(drop=True, inplace=True)

        # # Group by sample
        # df_species = df_species.groupby('SAMPLE').agg({
        #     'GENE': lambda x: ', '.join(x)
        #     'SUBTYPE': 
        # }).reset_index()

        # Generate HTML code for the table
        html_table = df_species.to_html(
            classes='table table-striped table-bordered',
            table_id=f'virulence_table_{sp}',
            index=False,
            border=0,
            justify='left'
        )

        # Store in the dictionary
        tables[species_name] = html_table

    return tables


def build_report(html_template_path, input_dir, output_dir):
    """
    Builds the report by generating visualizations and integrating them into the HTML template.

    Parameters:
    - html_template_path (str): Path to the directory containing the HTML template.

    Returns:
    - report_html (str): The final report as an HTML string.
    """
    # Load and preprocess data
    df = load_and_preprocess_data()

    # Check if figures directory exists
    if not os.path.exists("report/figures"):
        os.makedirs("report/figures")

    species_color_mapping = assign_colors(species_dict)

    # Create Plotly figures and collect their JSON data
    figures = []

    info_table = generate_general_info_table(df, input_dir, output_dir)

    sunburst_chart_fig = create_sunburst_chart(df, species_dict, species_color_mapping)
    sunburst_chart = create_figure_json(sunburst_chart_fig, "Sequence typing: sunburst chart", "sunburst_chart")
    sunburst_table = create_figure_with_table_json(
        sunburst_chart_fig,
        info_table,
        title="Sequence typing: sunburst chart",
        fig_id="sunburst_chart"
    )
    figures.append(sunburst_table)

    # Subtype Pie Charts
    subtype_pie_charts = create_subtype_polar_charts(df, species_dict)
    subtype_figures = []
    figure_table_pairs = []

    for sp, fig in subtype_pie_charts.items():
        # Generate the HTML for the table (you'll need to implement this)
        table_html = generate_table_html_for_species(sp, df)  # Implement this function as needed

        # Create the figure and table JSON
        figure_table = create_figure_with_table_json(
            fig,
            table_html,
            title=f"{species_dict[sp]}",
            fig_id=f"{sp}"
        )

        figure_table_pairs.append(figure_table)
        figures.append(figure_table)  # Add to the figures list

    # Pangenome Heatmaps
    pangenome_heatmaps_json = []
    pangenome_heatmaps = create_pangenome_heatmaps(species_dict)
    for sp, fig in pangenome_heatmaps.items():
        pangenome_heatmap = create_figure_json(fig, "", f"pangenome_{sp}")
        pangenome_heatmaps_json.append(pangenome_heatmap)
        figures.append(pangenome_heatmap)

    # Pangenome Pie Charts
    pangenome_pie_charts_json = []
    pangenome_pie_charts = create_pangenome_pie_charts(species_dict)
    for sp, fig in pangenome_pie_charts.items():
        pangenome_pie_chart = create_figure_json(fig, "", f"pangenome_pie_{sp}")
        pangenome_pie_charts_json.append(pangenome_pie_chart)
        figures.append(pangenome_pie_chart)

    # Combine pangenome elements. intercalate pangenome_heatmaps and pangenome_pie_charts
    pangenome_elements = [elem for pair in zip(pangenome_heatmaps_json, pangenome_pie_charts_json) for elem in pair]

    # Assuming you have your DataFrame 'df' and 'species_dict'
    tables = create_annotation_tables(df, species_dict)
    resistance_table_html = tables['Resistance']
    virulence_table_html = tables['Virulence']
    plasmid_table_html = tables['Plasmid']

    contig_line_plots, n50_contig_length_box_plot, contig_coverage_box_plot = create_contig_plot(df, species_dict)

    qc_elements = []

    # Create JSON data for the N50 box plot
    n50_contig_length_box_plot_json = create_figure_json(n50_contig_length_box_plot, "N50 Contig Length Distribution", "n50_contig_length_box_plot")
    figures.append(n50_contig_length_box_plot_json)
    qc_elements.append(n50_contig_length_box_plot_json)

    # Create JSON data for the contig coverage box plot
    contig_coverage_box_plot_json = create_figure_json(contig_coverage_box_plot, "Average Contig Coverage Distribution", "contig_coverage_box_plot")
    figures.append(contig_coverage_box_plot_json)
    qc_elements.append(contig_coverage_box_plot_json)

    # Create JSON data for the contig plots
    contig_line_plots_json = []
    for sample, fig in contig_line_plots.items():
        contig_line_plot = create_figure_json(fig, f"Contig Lengths for {sample}", f"contig_line_plot_{sample}")
        contig_line_plots_json.append(contig_line_plot)
        qc_elements.append(contig_line_plot)
        figures.append(contig_line_plot)

    fastqc_table_html = create_fastqc_table(df, species_dict, gc_content_dict)

    virulence_figures = generate_virulence_factor_bar_chart(df, species_dict)
    virulence_figure = {}
    for sp, fig in virulence_figures.items():
        virulence_figure = create_figure_json(fig, f"Virulence Factors for {sp}", f"virulence_{sp}")
        figures.append(virulence_figure)
    
    virulence_figure_tables = []
    virulence_tables = generate_virulence_factor_table(df, species_dict)
    for sp, table in virulence_tables.items():
        virulence_table = create_figure_with_table_json(
            virulence_figures[sp],
            table,
            title=f"Virulence Factors for {sp}",
            fig_id=f"virulence_{sp}"
        )
        virulence_figure_tables.append(virulence_table)

        figures.append(virulence_table)

    # Prepare the list of tabs
    tabs = [
        {'id': 'sunburst', 'title': 'MLST', 'figure_tables': [sunburst_table], 'title_text': 'MLST', 'help': sunburst_help_text_html},
        {'id': 'fastqc-table', 'title': 'FastQC', 'content': fastqc_table_html , 'title_text': 'FastQC', 'help': fastqc_help_text_html},
        {'id': 'assembly_qc', 'title': 'Assembly QC', 'fig_ids': [fig['fig_id'] for fig in qc_elements], 'titles': [fig['title'] for fig in qc_elements], 'title_text': 'Assembly QC', 'help': assembly_qc_help_text_html},
        {'id': 'resistance-table', 'title': 'Resistance annotations', 'content': resistance_table_html, 'title_text': 'Resistance annotations'},
        {'id': 'virulence-table', 'title': 'Virulence annotations', 'content': virulence_table_html, 'title_text': 'Virulence annotations'},
        {'id': 'plasmid-table', 'title': 'Plasmid annotations', 'content': plasmid_table_html, 'title_text': 'Plasmid annotations'},
        {
            'id': 'virulence',
            'title': 'Virulome',
            'figure_tables': virulence_figure_tables,
            'title_text': 'Virulome'
        },
        {
            'id': 'resistome',
            'title': 'Resistome',
            'figure_tables': figure_table_pairs,
            'title_text': 'Resistome',
            'help': resistome_help_text_html
        },
        {'id': 'pangenome', 'title': 'Pangenome', 'fig_ids': [fig['fig_id'] for fig in pangenome_elements], 'titles': [fig['title'] for fig in pangenome_elements], 'title_text': 'Pangenome'}
    ]   

    # Load the HTML template using Jinja2
    env = Environment(loader=FileSystemLoader(html_template_path))
    env.globals.update(zip=zip)
    template = env.get_template('layout.html')

    # Render the template with the variables
    report_html = template.render(
        tabs=tabs,
        figures=figures
    )

    # Write the report to a file
    with open(f"report/report.html", 'w') as file:
        file.write(report_html)
    os.system("touch flags/.report")

    return report_html
