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

from itertools import permutations
from bs4 import BeautifulSoup as bs
from collections import Counter
from jinja2 import Template

from utils import *
from itertools import chain


def build_report(html_string):

    matplotlib.use('Agg')

    # ========================= OVERVIEW =============================================================

    # Read results.csv
    df = pd.read_csv("dataframe/results.csv")

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

    # Assigning colors to species
    species_color_mapping = assign_colors(species_dict)

    # Create dataframe that has sample names as rows and genes as columns
    # At the moment the dataframe has the following structure:
    # [GENE, SAMPLE, ID, COVERAGE]

    # First, select only the plasmids rows
    df_virulence = df[df["PREDICTION_SOURCE"] == "VFDB"]

    # Only use the gene, sample, identity and coverage columns
    df_virulence_id = df_virulence[["GENE", "SAMPLE", "%IDENTITY"]]

    # Add the index as a separate column
    df_virulence["ID"] = df_virulence_id.index

    # Make a dictionary mapping each id to the gene
    id_to_gene = dict(zip(df_virulence_id.index, df_virulence_id["GENE"]))

    # Pivot the dataframe on the ID
    df_virulence_id = df_virulence.pivot_table(index='SAMPLE', columns='ID', values='%IDENTITY', aggfunc='max')
    df_virulence_cov = df_virulence.pivot_table(index='SAMPLE', columns='ID', values='%COVERAGE', aggfunc='max')

    # Delete rows tha only hae NaN values
    df_virulence_id = df_virulence_id.fillna(0)
    df_virulence_cov = df_virulence_cov.fillna(0)

    colorscale = "Blues"

    data_identity = go.Heatmap(
        z=df_virulence_id.values,
        x=[id_to_gene[item] for item in df_virulence_id.columns],
        y=df_virulence_id.index,
        colorscale=colorscale,
        hoverongaps = False,
        hovertemplate='Isolate: %{y}<br>' +
                    'Gene: %{x}<br>' +
                    'Identity: %{z}%<br><extra></extra>'
    )

    data_coverage = go.Heatmap(
        z=df_virulence_cov.values,
        x=[id_to_gene[item] for item in df_virulence_cov.columns],
        y=df_virulence_cov.index,
        colorscale=colorscale,
        hoverongaps = False,
        hovertemplate='Isolate: %{y}<br>' +
                    'Gene: %{x}<br>' +
                    'Coverage: %{z}%<br><extra></extra>'
    )

    heatmap_virulence_full_figure_coverage = plotly.tools.make_subplots(rows=2, cols=1, shared_xaxes=False)
    heatmap_virulence_full_figure_coverage.append_trace(data_identity, 1, 1)
    heatmap_virulence_full_figure_coverage.append_trace(data_coverage, 2, 1)

    heatmap_virulence_full_figure_coverage.update_layout(
        title={
            'text': "Virulence factors heatmap: identity & coverage",
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 20}
        },
        height=max(300 + len(df_virulence_id.index)*10, 600)
    )

    heatmap_virulence_full_figure_coverage = fix_colorbar(heatmap_virulence_full_figure_coverage, facet=True, dtick=10)

    heatmap_virulence_full_figure_coverage_code = create_html_element(heatmap_virulence_full_figure_coverage, "VFDB output - coverage heatmap")

    # _____________________ Plasmids quality heatmaps ________________________

    # First, select only the plasmids rows
    df_resistance = df[df["PREDICTION_SOURCE"] == "PlasmidFinder"]

    # Only use the gene, sample, identity and coverage columns
    df_resistance_id = df_resistance[["GENE", "SAMPLE", "%IDENTITY"]]

    # Add the index as a separate column
    df_resistance["ID"] = df_resistance_id.index

    # Make a dictionary mapping each id to the gene
    id_to_gene = dict(zip(df_resistance_id.index, df_resistance_id["GENE"]))

    # Pivot the dataframe on the ID
    df_plasmids_id = df_resistance.pivot_table(index='SAMPLE', columns='ID', values='%IDENTITY', aggfunc='max')
    df_plasmids_cov = df_resistance.pivot_table(index='SAMPLE', columns='ID', values='%COVERAGE', aggfunc='max')

    df_plasmids_id = df_plasmids_id.fillna(0)
    df_plasmids_cov = df_plasmids_cov.fillna(0)

    data_identity = go.Heatmap(
        z=df_plasmids_id.values,
        x=[id_to_gene[item] for item in df_plasmids_id.columns],
        y=df_plasmids_id.index,
        colorscale=colorscale,
        hoverongaps = False,
        hovertemplate='Isolate: %{y}<br>' +
                    'Gene: %{x}<br>' +
                    'Identity: %{z}%<br><extra></extra>'
    )

    # Don't show tiles for the second heatmap
    data_coverage = go.Heatmap(
        z=df_plasmids_cov.values,
        x=[id_to_gene[item] for item in df_plasmids_cov.columns],
        y=df_plasmids_cov.index,
        colorscale=colorscale,
        hoverongaps = False,
        hovertemplate='Isolate: %{y}<br>' +
                    'Gene: %{x}<br>' +
                    'Coverage: %{z}%<br><extra></extra>'
    )

    heatmap_plasmids_full_figure_coverage = plotly.tools.make_subplots(rows=2, cols=1, shared_xaxes=True)
    heatmap_plasmids_full_figure_coverage.append_trace(data_identity, 1, 1)
    heatmap_plasmids_full_figure_coverage.append_trace(data_coverage, 2, 1)

    heatmap_plasmids_full_figure_coverage.update_layout(
        title={
            'text': "Plasmids heatmap: identity & coverage",
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 20}
        },
        height=max(300 + len(df_virulence_id.index)*10, 600)
    )

    heatmap_plasmids_full_figure_coverage = fix_colorbar(heatmap_plasmids_full_figure_coverage, facet=True, dtick=10)

    heatmap_plasmids_full_figure_coverage_code = create_html_element(heatmap_plasmids_full_figure_coverage, "PlasmidFinder output - coverage heatmap")


    # _____________________ Resistance quality heatmaps ________________________

    # First, select only the plasmids rows
    df_resistance = df[df["PREDICTION_SOURCE"] == "NCBI"]

    # Only use the gene, sample, identity and coverage columns
    df_resistance_id = df_resistance[["GENE", "SAMPLE", "%IDENTITY"]]

    # Add the index as a separate column
    df_resistance["ID"] = df_resistance_id.index

    # Make a dictionary mapping each id to the gene
    id_to_gene = dict(zip(df_resistance_id.index, df_resistance_id["GENE"]))

    # Pivot the dataframe on the ID
    df_resistance_id = df_resistance.pivot_table(index='SAMPLE', columns='ID', values='%IDENTITY', aggfunc='max')
    df_resistance_cov = df_resistance.pivot_table(index='SAMPLE', columns='ID', values='%COVERAGE', aggfunc='max')

    df_resistance_id = df_resistance_id.fillna(0)
    df_resistance_cov = df_resistance_cov.fillna(0)

    data_identity = go.Heatmap(
        z=df_resistance_id.values,
        x=[id_to_gene[item] for item in df_resistance_id.columns],
        y=df_resistance_id.index.to_list(),
        colorscale=colorscale,
        hovertemplate='Isolate: %{y}<br>' +
                    'Gene: %{x}<br>' +
                    'Identity: %{z}%<br><extra></extra>'
    )

    data_coverage = go.Heatmap(
        z=df_resistance_cov.values,
        x=[id_to_gene[item] for item in df_resistance_cov.columns],
        y=df_resistance_cov.index,
        colorscale=colorscale,
        hovertemplate='Isolate: %{y}<br>' +
                    'Gene: %{x}<br>' +
                    'Coverage: %{z}%<br><extra></extra>'
    )

    heatmap_resistance_full_figure = plotly.tools.make_subplots(rows=2, cols=1, shared_xaxes=True)
    heatmap_resistance_full_figure.append_trace(data_identity, 1, 1)
    heatmap_resistance_full_figure.append_trace(data_coverage, 2, 1)

    # heatmap_plasmids_full_figure_coverage = go.Figure(data=[data])

    heatmap_resistance_full_figure.update_layout(
        title={
            'text': "Resistance factors heatmap: identity & coverage",
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 20}
        },
        height=max(300 + len(df_virulence_id.index)*10, 600)
    )

    heatmap_resistance_full_figure = fix_colorbar(heatmap_resistance_full_figure, facet=True, dtick=10)

    heatmap_resistance_full_figure_code = create_html_element(heatmap_resistance_full_figure, "Resistance output - identity & coverage heatmap")

    # if len(species_list) > 1:

    #     #____________________ Virulence heatmap __________________

    #____________Sunburst Chart________________

    label = []
    parent = []
    value = []
    node_info = []
    colors = []

    df_samples = df.groupby(["SAMPLE", "SUBTYPE", "SPECIES"]).size().reset_index(name="count")
    sample_list = df_samples["SAMPLE"].unique().tolist()

    for _, sample in df_samples.iterrows():
        label.append(sample["SAMPLE"])
        parent.append(sample["SUBTYPE"] + "~" + sample["SPECIES"])
        value.append(1)
        if sample["SUBTYPE"] + "~" + sample["SPECIES"] not in label:
            label.append(sample["SUBTYPE"] + "~" + sample["SPECIES"])
            parent.append(sample["SPECIES"])
            value.append(0)
        if sample["SPECIES"] not in label:
            label.append(sample["SPECIES"])
            parent.append("")
            value.append(0)

    for par, item in zip(parent, label):

        if item in df_samples['SPECIES'].to_list():
            no_of_samples = df_samples['SPECIES'].to_list().count(item)
            node_info.append(f"Number of isolates: {no_of_samples} <br> {round(no_of_samples/len(sample_list), 2)*100}% from total number of isolates")
            colors.append(species_color_mapping[species_dict[item]])

        elif item.split("~")[0] in df_samples['SUBTYPE'].to_list():
            subtype = item.split("~")[0]
            species = item.split("~")[1]
            no_of_samples = df_samples[df_samples['SPECIES'] == par]['SUBTYPE'].to_list().count(subtype)
            node_info.append(f"Number of isolates: {no_of_samples} <br> {round(no_of_samples/len(sample_list), 2)*100}% from total from total number of isolates")
            # species_corresponding_subtype = df_samples[df_samples["subtype"]==item.split("_")[0]]["species"].to_list()[0]
            diluted_color = dilute_hex_color(species_color_mapping[species_dict[species]], 0.05)
            colors.append(diluted_color)

        elif item in df_samples['SAMPLE'].to_list():
            node_info.append("")
            species = par.split("~")[1]
            # species_corresponding_subtype = df_samples[df_samples["name"]==item]["species"].to_list()[0]
            diluted_color = dilute_hex_color(species_color_mapping[species_dict[species]], 0.1)
            colors.append(diluted_color)

    import pdb; pdb.set_trace()

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
        customdata=node_info,  # pass the customdata
        hovertemplate='<b>%{label}</b><br>%{customdata}<extra></extra>',
        marker=dict(
            colors=colors
        )
    )

    sunburst_figure = go.Figure(data=sunburst_trace)
    sunburst_figure.update_layout(height=800)
    sunburst_figure.update_layout(
        hoverlabel_namelength=0,
        title="Sequence typing: sunburst chart"
    )
    sunburst_figure_code = create_html_element(sunburst_figure, "MLST output - Sunburst chart")

    # # ========================= SAMPLES INFO =========================================================

    # ______________ Table _____________________

    accession_product = dict(zip(df['ACCESSION'], df['GENE']))
    sample_species = dict(zip(df['SAMPLE'], df['SPECIES']))
    sample_subtype = dict(zip(df['SAMPLE'], df['SUBTYPE']))

    # Make a dictionary with the aggregate functions to be used in the groupby
    agg_funcs = {
        'GENE': lambda x: list(x),
        'SUBTYPE': 'first',
        'SPECIES': 'first',
        'CONTIG_LENGTH': 'first',
        'CONTIG_COVERAGE': 'first',
        'CONTIG_NUMBER': 'first'
    }

    # Pivot table on teh values of the PREDICTION_SOURCE column
    df_table = df_table = df.pivot_table(index='SAMPLE',columns='PREDICTION_SOURCE',values='GENE',aggfunc=lambda x: x.tolist())
    df_table["SPECIES"] = df_table.index.map(sample_species)
    df_table["SUBTYPE"] = df_table.index.map(sample_subtype)
    df_table["SPECIES"] = df_table["SPECIES"].map(species_dict)

    # Turn NaN to empty lists
    df_table = df_table.fillna('')

    # Convert lists to comma separated strings
    df_table['NCBI'] = df_table['NCBI'].apply(lambda x: ", ".join(x))
    df_table['VFDB'] = df_table['VFDB'].apply(lambda x: ", ".join(x))
    df_table['PlasmidFinder'] = df_table['PlasmidFinder'].apply(lambda x: ", ".join(x))

    df_full_table_code = df_table.to_html()

    # Parse the HTML code using Beautiful Soup
    soup = bs(df_full_table_code, 'html.parser')

    # Find the table element
    df_full_table_code = soup.find('table')

    # Add a class and an ID to the table element
    df_full_table_code['class'] = 'table table-striped'
    df_full_table_code['id'] = 'full_vpr_table'

    for th in df_full_table_code.find_all('th'):
        th['style'] = 'text-align: left'

    df_full_table_code = create_html_widget(df_full_table_code, "VFDB, Plasmid finder, Resistance output - table")


    # ========================= SUBTYPES INFO ========================================================

    subtype_html_string = ""
    df_genes_entry_full = None

    for sp in df["SPECIES"].unique().tolist():

        # Group by species and agregate by making lists of the values
        df_species = df[df["SPECIES"] == sp]

        # Keep only the ["ACCESION", "PRODUCT", "SUBTYPE", "RESISTANCE"] columns
        df_species = df_species[["ACCESSION", "PRODUCT", "SUBTYPE", "RESISTANCE"]]

        df_species = df_species.dropna(subset=["RESISTANCE"])

        if not len(df_species):
            continue

        PC_subtype_subtype_subplot = px.pie(df_species, names='RESISTANCE', hole=0.65, color='RESISTANCE', facet_col='SUBTYPE', facet_col_wrap=4)
        PC_subtype_subtype_subplot.update_traces(textinfo='label')
        PC_subtype_subtype_subplot = configure_pie_chart(PC_subtype_subtype_subplot)
        PC_subtype_subtype_subplot.update_layout(
            height=300*int((len(df_species["SUBTYPE"].unique())+3)/4),
            legend=dict(
                orientation='v',
                yanchor='top',
                y=0.99,
                xanchor='right',
                x=0.99
            ),
            showlegend=False
        )

        PC_subtype_subtype_subplot.update_traces(
                hovertemplate='number of factors pertaining to <b>%{label}</b> resistance: %{value}',
                textposition='inside'
        )

        PC_subtype_subtype_subplot.write_html(f"report/figures/PC_subtype_subtype_subplot.html", full_html=False, include_plotlyjs='cdn')
        f = open(f"report/figures/PC_subtype_subtype_subplot.html")
        PC_subtype_subtype_subplot = f.read()

        subtype_html_string += create_html_card(PC_subtype_subtype_subplot, species_dict[sp])


    df_genes_entry_full = df_species.drop('SUBTYPE', axis=1)
    df_genes_entry_full.to_csv(f"report/full_genes_table.csv", index=False)

    # ___________Pangenome______________

    Heatmap_pangenomic_codes = ""

    for sp in [x for x in os.listdir("pangenome") if os.path.isdir(f"pangenome/{x}")]:

        if "output" not in os.listdir(f"pangenome/{sp}"):
            continue

        # Read the gene_presence_absence.csv file from Roary output
        roary_csv_path = find_file("gene_presence_absence.csv", f"pangenome/{sp}")
        df_pangenome = pd.read_csv(roary_csv_path)
        # Set the index of the DataFrame to the gene names
        df_pangenome = df_pangenome.set_index("Gene")
        df_pangenome = df_pangenome.drop(df_pangenome.columns[:13], axis=1)
        # replace NaN values with 0 and any other value with 1
        df_pangenome = df_pangenome.fillna(0).applymap(lambda x: 0 if x==0 else 1)
        # Transpose the DataFrame so that speciess are rows and genes are columns
        df_pangenome = df_pangenome.T
        df_pangenome = df_pangenome[df_pangenome.sum().sort_values(ascending=False).index]
        df_pangenome = df_pangenome.loc[(df_pangenome.sum(axis=1) != 0), :]
        # Create the heatmap
        Heatmap_pangenomic = go.Figure(
            data=[
                go.Heatmap(
                    x=df_pangenome.columns,
                    y=df_pangenome.index,
                    z=df_pangenome.values,
                    colorscale='Blues',
                    hoverongaps = False,
                    hovertemplate=
                    'VALUE: %{z}<br>' +  # Display the value of the cell
                    'GENE: %{x}<br>' +  # Display the label of the x-axis
                    'SAMPLE: %{y}<br>' +  # Display the label of the y-axis
                    '<extra></extra>'  # Hide the secondary box
                )],
            layout=go.Layout(title="Gene Presence and Absence", xaxis_title="GENES", yaxis_title=f"{species_dict[sp]}")
        )

        # Update the layout of the heatmap, update title and axist title font size
        # Also set title
        Heatmap_pangenomic.update_layout(
            height=max(200 + len(df_pangenome.index)*10, 400),
            title_font_size=20,
            xaxis=dict(title_font_size=18),
            yaxis=dict(title_font_size=18)
        )

        Heatmap_pangenomic = fix_colorbar(Heatmap_pangenomic, facet=False)

        Heatmap_pangenomic_code = create_html_element(Heatmap_pangenomic, f"Roary output - pangenome heatmap - {species_dict[sp]}")

        Heatmap_pangenomic_codes += Heatmap_pangenomic_code

    # ========================= QUALITY CONTROL ======================================================

    df_samples = df.groupby("SAMPLE").agg({
        "CONTIG_LENGTH": "first",
        "CONTIG_COVERAGE": "first",
        "SPECIES": "first",
        "CONTIG_NUMBER": "first"
    }).reset_index()

    # ______________________ Contig Plot _________________________
    Scatter_contig_length = px.scatter(df_samples, x='CONTIG_LENGTH', y='CONTIG_COVERAGE', color='SPECIES',
        color_discrete_map=species_color_mapping,
        hover_data=['CONTIG_LENGTH', 'CONTIG_COVERAGE', 'SPECIES', 'SAMPLE', 'CONTIG_NUMBER'],
    )

    Scatter_contig_length.update_traces(marker_size=25)

    # Update the length of the plot
    Scatter_contig_length.update_layout(
        height=1000
    )

    Scatter_contig_length_code = create_html_element(Scatter_contig_length, "Scatter plot for assembly contigs")

    # ______________________ Contig Box _________________________
    Box_contig_length = px.box(df_samples, y='CONTIG_LENGTH', x='SPECIES', color='SPECIES', color_discrete_map=species_color_mapping)
    Box_contig_length_code = create_html_element(Box_contig_length, "Box plot for assembly contigs")

    # Find all folders named fastqc in the results folder and subfolders recursively
    # without using find_file function
    # Extract the HTML files from the fastqc folders
    html_files = []
    for root, dirs, files in os.walk("results"):
        for file in files:
            if file.endswith(".html"):
                html_files.append(os.path.join(root, file))

    # Define the strings we're looking for
    search_strings = ["Total Sequences", "Sequences flagged as poor quality", "Sequence length", "%GC"]

    df_qc = pd.DataFrame(
        index=list(
            set(
                [os.path.splitext(os.path.basename(file_path))[0].replace("_L001_R1_001", "").replace("_L001_R2_001", "").replace("_fastqc", "") for file_path in html_files])),
        columns=[x+" R1" for x in search_strings]+[x+" R2" for x in search_strings]+['SPECIES']
    )

    # Get the unique index values
    unique_indexes = df_qc.index.unique()

    # Set the index of the DataFrame to the unique values
    df_qc = df_qc.set_index(unique_indexes)

    # Loop over each file path
    for file_path in html_files:

        # Load the HTML file into BeautifulSoup
        with open(file_path, "r") as f:
            soup = bs(f, "html.parser")

        file_name = os.path.splitext(os.path.basename(file_path))[0].replace("_L001_R1_001", "").replace("_L001_R2_001", "").replace("_fastqc", "")

        # Loop over each search string
        for search_string in search_strings:
            # Find the td element with the search string as its contents
            td = soup.find("td", string=search_string)

            # If the td element exists, extract the content of the next td element
            if td is not None:
                next_td = td.find_next_sibling("td")
                if next_td is not None:

                    if "R1" in file_path:
                        df_qc.loc[file_name, search_string+" R1"] = next_td.get_text().strip()

                    if "R2" in file_path:
                        df_qc.loc[file_name, search_string+" R2"] = next_td.get_text().strip()

    for i in df_qc.index.to_list():
        df_qc.loc[i, "SPECIES"] = species_dict[df_samples[df_samples["SAMPLE"]==i]["SPECIES"].values[0]]

    df_qc["Sequences flagged as poor quality R1"] = df_qc["Sequences flagged as poor quality R1"] + " (" + (df_qc["Sequences flagged as poor quality R1"].astype(int) / df_qc["Total Sequences R1"].astype(int) * 100).astype(str) + "%)"
    df_qc["Sequences flagged as poor quality R2"] = df_qc["Sequences flagged as poor quality R2"] + " (" + (df_qc["Sequences flagged as poor quality R2"].astype(int) / df_qc["Total Sequences R2"].astype(int) * 100).astype(str) + "%)"

    def compute_GC_intervalsR1(row):
        
        sp = row["SPECIES"]
        if sp not in gc_content_dict.keys():
            return row["%GC R1"] + "% (expected: N/A)"
        
        return row["%GC R1"] + "% (expected: " + gc_content_dict[sp] + "%)"
    
    df_qc["%GC R1"] = df_qc.apply(compute_GC_intervalsR1, axis=1)

    def compute_GC_intervalsR2(row):
        
        sp = row["SPECIES"]
        if sp not in gc_content_dict.keys():
            return row["%GC R2"] + "% (expected: N/A)"
        
        return row["%GC R2"] + "% (expected: " + gc_content_dict[sp] + "%)"
    
    df_qc["%GC R2"] = df_qc.apply(compute_GC_intervalsR2, axis=1)

    # df.to_csv(f"{TABLES_PATH}/df_reads_table.csv", index=False)
    df_reads_table_code = df_qc.to_html()

    # Parse the HTML code using Beautiful Soup
    soup = bs(df_reads_table_code, 'html.parser')

    # Find the table element
    df_reads_table_code = soup.find('table')

    # Add a class and an ID to the table element
    df_reads_table_code['class'] = 'table table-striped'
    df_reads_table_code['id'] = 'full_qa_table'

    # Define a function to apply colors based on the difference
    def apply_color(td):

        if "N/A" in td.get_text():
            td['style'] = 'background-color: #f8d7da'  # Red

        # Extract the actual and expected percentages using regex
        actual_value = float(td.get_text().split('%')[0])  # Get the actual percentage value
        expected_value = float(td.get_text().split('%')[1].split(" ")[-1])   # Get the expected percentage value
        
        # Calculate the difference
        difference = abs(actual_value - expected_value)
        
        # Apply styles based on the difference
        if difference > 6:
            td['style'] = 'background-color: #f8d7da'  # Red
        else:
            td['style'] = 'background-color: #d4edda'  # Green

    # Iterate through the rows and find the relevant cells
    for row in df_reads_table_code.find_all('tr'):
        for td in row.find_all('td'):
            # Check if the cell matches the format "49% (expected: 50.8%)"
            if "expected" in td.get_text():
                apply_color(td)


    # Define a function to apply colors based on the percentage
    def apply_color(td):
        # Extract the percentage from the text using regex
        match = re.search(r'(\d+) \((\d+\.\d+)%\)', td.get_text())
        if match:
            percent_value = float(match.group(2))  # Get the float percentage value
            
            # Apply styles based on the percentage value
            if percent_value < 5:
                td['style'] = 'background-color: #d4edda'  # Green
            elif 5 <= percent_value < 15:
                td['style'] = 'background-color: #fff3cd'  # Yellow
            elif percent_value >= 15:
                td['style'] = 'background-color: #f8d7da'  # Red

    # Iterate through the rows and find the relevant cells
    for row in df_reads_table_code.find_all('tr'):
        for td in row.find_all('td'):
            # Check if the cell matches the format "integer (float%)"
            if re.search(r'^\d+ \(\d+\.\d+%\)$', td.get_text()):
                apply_color(td)

    for th in df_reads_table_code.find_all('th'):
        th['style'] = 'text-align: left'

    df_reads_table_code = create_html_widget(df_reads_table_code, "FastQC output - table")

    pangenome_pie_chart_codes = ""

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

        # Create a Pie chart trace
        trace = go.Pie(labels=df_pangenome.index, values=df_pangenome['Value'])

        # Create a layout with a title and set height
        # Set the font size of the title to 20
        # Place the title in the center of the chart, on the left
        layout = go.Layout(
            title=f'{species_dict[sp]}',
            title_x=0.1,
            title_y=0.5,
            height=400,
            title_font_size=20
        )

        # Create a Figure object and add the trace and layout
        pangenome_pie_chart = go.Figure(data=[trace], layout=layout)
        configure_pie_chart(pangenome_pie_chart)

        pangenome_pie_chart_code = create_html_element(pangenome_pie_chart, f"Pangenome Summary - {species_dict[sp]}")

        pangenome_pie_chart_codes += pangenome_pie_chart_code


    #================================================================================================================================================
    #============================================ BUILD HTML TEMPLATE ===============================================================================
    #================================================================================================================================================

    # Define the HTML template string
    template = Template(html_string)

    # Define the visualizations dictionary with descriptive titles as keys and visualization codes as values
    visualizations = {
        "MLST Sunburst Plot": sunburst_figure_code,
        "Virulence Heatmap": heatmap_virulence_full_figure_coverage_code,
        "Plasmid Gene Heatmap": heatmap_plasmids_full_figure_coverage_code,
        "Resistance Heatmap": heatmap_resistance_full_figure_code,
        "Resistance Profile": subtype_html_string,
        "Pangenomic Heatmap": Heatmap_pangenomic_codes,
        "Contig Length Scatter Plot": Scatter_contig_length_code,
        "Contig Length Box Plot": Box_contig_length_code,
        "Gene Table": df_full_table_code,
        "FastQC Table": df_reads_table_code,
        "Pangenomic Pie Chart": pangenome_pie_chart_codes
    }

    # Iterate over the items in the visualizations dictionary
    for title, visualization in visualizations.items():

        title_html_string = f"""
            <div class="d-sm-flex align-items-center justify-content-between mb-4">
                <h3 class="h3 mb-0 text-gray-800 mid-height">{title}</h3>
            </div>
        """

        # Render the HTML template with the title and visualization as template variables
        rendered_template = template.render(
            Title_html_string=title_html_string,
            Content_html_string=visualization
        )

        # Open a file with the derived file name and write the rendered template to it
        with open(f"report/{title.replace(' ', '_').replace('_Code', '')}.html", 'w') as file:
            file.write(rendered_template)

    os.system("touch flags/.report")