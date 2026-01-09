import json
import os

from jinja2 import Environment, FileSystemLoader

from report_charts import create_subtype_polar_charts, generate_virulence_factor_bar_chart
from report_contigs import create_contig_plot
from report_data import generate_general_info_table, load_and_preprocess_data
from report_mlst import create_sunburst_chart
from report_pangenome import create_pangenome_heatmaps, create_pangenome_pie_charts
from report_tables import (
    create_annotation_tables,
    create_fastqc_table,
    generate_table_html_for_species,
    generate_virulence_factor_table,
)
from utils import (
    assembly_qc_help_text_html,
    assign_colors,
    create_figure_json,
    create_figure_with_table_json,
    fastqc_help_text_html,
    resistome_help_text_html,
    species_dict,
    sunburst_help_text_html,
)


def build_report(html_template_path, input_dir, output_dir, gc_path):
    """
    Builds the report by generating visualizations and integrating them into the HTML template.

    Parameters:
    - html_template_path (str): Path to the directory containing the HTML template.
    - input_dir (str): Path to the input directory.
    - output_dir (str): Path to the output directory.
    - gc_path (str): Path to the GC content JSON file.

    Returns:
    - report_html (str): The final report as an HTML string.
    """
    # Load and preprocess data.
    df = load_and_preprocess_data()

    # Check if the figures directory exists.
    if not os.path.exists("report/figures"):
        os.makedirs("report/figures")

    species_color_mapping = assign_colors(species_dict)

    # Create Plotly figures and collect their JSON data.
    figures = []

    info_table = generate_general_info_table(df, input_dir, output_dir)

    sunburst_chart_fig = create_sunburst_chart(df, species_dict, species_color_mapping, output_dir)
    sunburst_table = create_figure_with_table_json(
        sunburst_chart_fig,
        info_table,
        title="Sequence typing: sunburst chart",
        fig_id="sunburst_chart"
    )
    figures.append(sunburst_table)

    # Subtype pie charts.
    subtype_pie_charts = create_subtype_polar_charts(df, species_dict)
    figure_table_pairs = []

    for sp, fig in subtype_pie_charts.items():
        # Generate the HTML for the table (you'll need to implement this).
        table_html = generate_table_html_for_species(sp, df)  # Implement this function as needed

        # Create the figure and table JSON.
        figure_table = create_figure_with_table_json(
            fig,
            table_html,
            title=f"{species_dict[sp]}",
            fig_id=f"{sp}"
        )

        figure_table_pairs.append(figure_table)
        figures.append(figure_table)  # Add to the figures list.

    # Pangenome heatmaps.
    pangenome_heatmaps_json = []
    pangenome_heatmaps = create_pangenome_heatmaps(species_dict)
    for sp, fig in pangenome_heatmaps.items():
        pangenome_heatmap = create_figure_json(fig, "", f"pangenome_{sp}")
        pangenome_heatmaps_json.append(pangenome_heatmap)
        figures.append(pangenome_heatmap)

    # Pangenome pie charts.
    pangenome_pie_charts_json = []
    pangenome_pie_charts = create_pangenome_pie_charts(species_dict)
    for sp, fig in pangenome_pie_charts.items():
        pangenome_pie_chart = create_figure_json(fig, "", f"pangenome_pie_{sp}")
        pangenome_pie_charts_json.append(pangenome_pie_chart)
        figures.append(pangenome_pie_chart)

    # Combine pangenome elements by interleaving heatmaps and pie charts.
    pangenome_elements = [elem for pair in zip(pangenome_heatmaps_json, pangenome_pie_charts_json) for elem in pair]

    # Assuming you have your DataFrame 'df' and 'species_dict'.
    tables = create_annotation_tables(df, species_dict)
    resistance_table_html = tables['Resistance']
    virulence_table_html = tables['Virulence']
    plasmid_table_html = tables['Plasmid']

    contig_line_plots, n50_contig_length_box_plot, contig_coverage_box_plot = create_contig_plot(df, species_dict)

    qc_elements = []

    # Create JSON data for the N50 box plot.
    n50_contig_length_box_plot_json = create_figure_json(n50_contig_length_box_plot, "N50 Contig Length Distribution", "n50_contig_length_box_plot")
    figures.append(n50_contig_length_box_plot_json)
    qc_elements.append(n50_contig_length_box_plot_json)

    # Create JSON data for the contig coverage box plot.
    contig_coverage_box_plot_json = create_figure_json(contig_coverage_box_plot, "Average Contig Coverage Distribution", "contig_coverage_box_plot")
    figures.append(contig_coverage_box_plot_json)
    qc_elements.append(contig_coverage_box_plot_json)

    # Create JSON data for the contig plots.
    for sample, fig in contig_line_plots.items():
        contig_line_plot = create_figure_json(fig, f"Contig Lengths for {sample}", f"contig_line_plot_{sample}")
        qc_elements.append(contig_line_plot)
        figures.append(contig_line_plot)

    # Read JSON file for GC content.
    gc_content_dict = {}
    with open(gc_path, 'r') as f:
        gc_content_dict = json.load(f)

    fastqc_table_html = create_fastqc_table(df, species_dict, gc_content_dict)

    virulence_figures = generate_virulence_factor_bar_chart(df, species_dict)
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

    # Prepare the list of tabs.
    tabs = [
        {'id': 'sunburst', 'title': 'MLST', 'figure_tables': [sunburst_table], 'title_text': 'MLST', 'help': sunburst_help_text_html},
        {'id': 'fastqc-table', 'title': 'FastQC', 'content': fastqc_table_html, 'title_text': 'FastQC', 'help': fastqc_help_text_html},
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

    # Load the HTML template using Jinja2.
    env = Environment(loader=FileSystemLoader(html_template_path))
    env.globals.update(zip=zip)
    template = env.get_template('layout.html')

    # Render the template with the variables.
    report_html = template.render(
        tabs=tabs,
        figures=figures
    )

    # Write the report to a file.
    with open(f"report/report.html", 'w') as file:
        file.write(report_html)
    return report_html
