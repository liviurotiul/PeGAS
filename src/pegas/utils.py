# Utility functions for handling files and database operations.
import os
import random

def generate_color():
    # Generate a random hue value between 0 and 360 (in degrees).
    h = random.randint(0, 359)
    # Set saturation and lightness values for a pastel effect.
    s = 55  # Lower saturation for pastel effect.
    l = 70  # Higher lightness for pastel effect.
    # Convert HSL to RGB and format as a hexadecimal color code.
    r, g, b = hsl_to_rgb(h, s, l)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def hsl_to_rgb(h, s, l):
    # Convert HSL values to floating-point numbers between 0 and 1.
    h = h / 360
    s = s / 100
    l = l / 100

    # Calculate intermediate values.
    if s == 0:
        r = g = b = l
    else:
        q = l * (1 + s) if l < 0.5 else l + s - l * s
        p = 2 * l - q
        r = hue_to_rgb(p, q, h + 1/3)
        g = hue_to_rgb(p, q, h)
        b = hue_to_rgb(p, q, h - 1/3)

    # Convert RGB values to integers between 0 and 255.
    r = round(r * 255)
    g = round(g * 255)
    b = round(b * 255)

    return r, g, b


def hue_to_rgb(p, q, t):
    if t < 0:
        t += 1
    elif t > 1:
        t -= 1

    if t < 1/6:
        return p + (q - p) * 6 * t
    elif t < 1/2:
        return q
    elif t < 2/3:
        return p + (q - p) * (2/3 - t) * 6
    else:
        return p


def assign_colors(old_dict):
    colors = set()
    new_dict = {}
    for value in old_dict.values():
        while True:
            color = generate_color()
            if color not in colors:
                colors.add(color)
                new_dict[value] = color
                break
    return new_dict


species_dict = {
    "liberibacter": "Liberibacter",
    "ecoli_achtman_4": "Escherichia coli",
    "schromogenes": "Streptomyces chromogenes",
    "bpseudomallei": "Burkholderia pseudomallei",
    "campylobacter_nonjejuni_3": "Campylobacter (non-Jejuni)",
    "mhaemolytica": "Mannheimia haemolytica",
    "szooepidemicus": "Streptococcus zooepidemicus",
    "orhinotracheale": "Ornithorhynchus anatinus",
    "salmonella": "Salmonella enterica",
    "senterica_achtman_2": "Salmonella enterica",
    "llactis_phage": "Lactococcus lactis phage",
    "magalactiae": "Mycoplasma gallisepticum",
    "sepidermidis": "Staphylococcus epidermidis",
    "vcholerae": "Vibrio cholerae",
    "mhominis_3": "Mycoplasma hominis",
    "psalmonis": "Pasteurella salmonis",
    "brachyspira_4": "Brachyspira hyodysenteriae",
    "hcinaedi": "Helicobacter cinaedi",
    "neisseria": "Neisseria meningitidis",
    "efaecalis": "Enterococcus faecalis",
    "bcc": "Burkholderia cepacia",
    "pgingivalis": "Porphyromonas gingivalis",
    "sbsec": "Streptococcus bovis",
    "shewanella": "Shewanella oneidensis",
    "leptospira": "Leptospira interrogans",
    "mhyopneumoniae": "Mycoplasma hyopneumoniae",
    "scanis": "Streptococcus canis",
    "cfreundii": "Citrobacter freundii",
    "yruckeri": "Yersinia ruckeri",
    "mgallisepticum": "Mycoplasma gallisepticum",
    "kaerogenes": "Klebsiella aerogenes",
    "campylobacter_nonjejuni_8": "Campylobacter lari",
    "msciuri": "Staphylococcus sciuri",
    "fpsychrophilum": "Flavobacterium psychrophilum",
    "mplutonius": "Mycoplasma plutonius",
    "ecloacae": "Enterobacter cloacae",
    "lsalivarius": "Lactobacillus salivarius",
    "saureus": "Staphylococcus aureus",
    "campylobacter_nonjejuni_7": "Campylobacter upsaliensis",
    "campylobacter": "Campylobacter jejuni",
    "leptospira_2": "Leptospira interrogans",
    "dnodosus": "Desulfovibrio nodosus",
    "ranatipestifer": "Riemerella anatipestifer",
    "brachyspira_3": "Brachyspira",
    "pacnes_3": "Propionibacterium acnes",
    "cmaltaromaticum": "Clostridium maltaromaticum",
    "campylobacter_nonjejuni_5": "Campylobacter",
    "oralstrep": "Streptococcus oralis",
    "campylobacter_nonjejuni_9": "Campylobacter",
    "abaumannii": "Acinetobacter baumannii",
    "abaumannii_2": "Acinetobacter baumannii",
    "sagalactiae": "Streptococcus agalactiae",
    "campylobacter_nonjejuni_4": "Campylobacter",
    "leptospira_3": "Leptospira",
    "aphagocytophilum": "Anaplasma phagocytophilum",
    "cperfringens": "Clostridium perfringens",
    "helicobacter": "Helicobacter pylori",
    "spyogenes": "Streptococcus pyogenes",
    "vparahaemolyticus": "Vibrio parahaemolyticus",
    "cdifficile": "Clostridioides difficile",
    "streptomyces": "Streptomyces",
    "ypseudotuberculosis_achtman_3": "Yersinia pseudotuberculosis",
    "mpneumoniae": "Mycoplasma pneumoniae",
    "achromobacter": "Achromobacter",
    "bhenselae": "Bartonella henselae",
    "aeromonas": "Aeromonas",
    "wolbachia": "Wolbachia",
    "hsuis": "Haemophilus suis",
    "miowae": "Mycoplasma iowae",
    "suberis": "Streptococcus suberis",
    "tpallidum": "Treponema pallidum",
    "mcaseolyticus": "Mycoplasma caseolyticus",
    "taylorella": "Taylorella",
    "mabscessus": "Mycobacterium abscessus",
    "shaemolyticus": "Staphylococcus haemolyticus",
    "mbovis_2": "Mycobacterium bovis",
    "edwardsiella": "Edwardsiella",
    "pmultocida_2": "Pasteurella multocida",
    "kingella": "Kingella",
    "sdysgalactiae": "Streptococcus dysgalactiae",
    "pputida": "Pseudomonas putida",
    "xfastidiosa": "Xylella fastidiosa",
    "pdamselae": "Photobacterium damselae",
    "campylobacter_nonjejuni": "Campylobacter species other than C. jejuni",
    "ppentosaceus": "Pediococcus pentosaceus",
    "shominis": "Sphaerochaeta multiformis",
    "ecoli": "Escherichia coli",
    "geotrichum": "Geotrichum species",
    "mcatarrhalis_achtman_6": "Moraxella catarrhalis (Achtman 6)",
    "klebsiella": "Klebsiella pneumoniae",
    "vcholerae_2": "Vibrio cholerae (El Tor biotype)",
    "vibrio": "Vibrio species",
    "hparasuis": "Haemophilus parasuis",
    "ssuis": "Streptococcus suis",
    "pmultocida": "Pasteurella multocida",
    "spneumoniae": "Streptococcus pneumoniae",
    "diphtheria_3": "Corynebacterium diphtheriae (biovar mitis)",
    "cronobacter": "Cronobacter sakazakii",
    "mgallisepticum_2": "Mycoplasma gallisepticum (strain 2)",
    "vtapetis": "Vibrio tapetis",
    "chlamydiales": "Chlamydiales order",
    "bwashoensis": "Bacillus washoeensis",
    "streptothermophilus": "Streptococcus thermophilus",
    "hinfluenzae": "Haemophilus influenzae",
    "brachyspira": "Brachyspira species",
    "spseudintermedius": "Staphylococcus pseudintermedius",
    "bbacilliformis": "Brevibacterium bacilliformis",
    "manserisalpingitidis": "Mycoplasma ansellae (formerly Mycoplasma eris) serovar 2",
    "campylobacter_nonjejuni_6": "Campylobacter species other than C. jejuni (strain 6)",
    "otsutsugamushi": "Orientia tsutsugamushi",
    "mflocculare": "Mycobacterium flocculare",
    "ureaplasma": "Ureaplasma species",
    "gallibacterium": "Gallibacterium species",
    "mycobacteria_2": "Mycobacterium species (strain 2)",
    "staphlugdunensis": "Staphylococcus lugdunensis",
    "mhyorhinis": "Mycoplasma hyorhinis",
    "bcereus": "Bacillus cereus",
    "csepticum": "Clostridium septicum",
    "sthermophilus": "Streptococcus thermophilus",
    "paeruginosa": "Pseudomonas aeruginosa",
    "listeria_2": "Listeria",
    "bfragilis": "Bacteroides fragilis",
    "borrelia": "Borrelia",
    "vvulnificus": "Vibrio vulnificus",
    "arcobacter": "Arcobacter",
    "pfluorescens": "Pseudomonas fluorescens",
    "sinorhizobium": "Sinorhizobium",
    "smaltophilia": "Stenotrophomonas maltophilia",
    "koxytoca": "Klebsiella oxytoca",
    "mcanis": "Mannheimia canis",
    "brucella": "Brucella",
    "campylobacter_nonjejuni_2": "Campylobacter (non-jejuni)",
    "bordetella_3": "Bordetella",
    "rhodococcus": "Rhodococcus",
    "brachyspira_5": "Brachyspira",
    "brachyspira_2": "Brachyspira",
    "cbotulinum": "Clostridium botulinum",
    "efaecium": "Enterococcus faecium",
    "aactinomycetemcomitans": "Aggregatibacter actinomycetemcomitans",
    "msynoviae": "Mycoplasma synoviae",
    "bsubtilis": "Bacillus subtilis",
    "plarvae": "Paenibacillus larvae",
    "sgallolyticus": "Streptococcus gallolyticus",
    "blicheniformis_14": "Bacillus licheniformis",
    "tenacibaculum": "Tenacibaculum",
    "unknown": "Unknown",
    "-": "Unknown",
}

def dilute_hex_color(hex_color, factor):
    # Convert hex color string to RGB tuple.
    hex_color = hex_color.lstrip('#')
    rgb_color = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

    # Apply lightening factor to each RGB component.
    r = min(int(rgb_color[0] + 255 * factor), 255)
    g = min(int(rgb_color[1] + 255 * factor), 255)
    b = min(int(rgb_color[2] + 255 * factor), 255)

    # Convert RGB tuple to hex color string.
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)


def create_html_widget(widget, title, width="8", min_height="300px"):

    collapse_id = "collapse_" + str(random.randint(0, 10**10))

    card_template = f"""
    <div class="row">
        <div class="col-xl-12">
            {widget}                    
        </div>
    </div>
    """
    return card_template


def configure_pie_chart(pie_chart_figure):

    pie_chart_figure.update_layout(
        height=340
    )

    pie_chart_figure.update_layout(
        font=dict(
            family="Courier New, monospace",
            size=20
        ),
        margin=dict(t=20, b=20, l=30, r=30),
    )

    return pie_chart_figure


def create_figure_with_table_json(widget, table_html, title, fig_id):
    """
    Converts a Plotly figure and an HTML table to JSON data for deferred rendering.

    Parameters:
    - widget (go.Figure): The Plotly figure.
    - table_html (str): The HTML string of the data table.
    - title (str): The title for the figure.
    - fig_id (str): A unique identifier for the figure.

    Returns:
    - dict: A dictionary containing the figure data, table HTML, and metadata.
    """
    # Ensure the figure has a unique ID.
    if not fig_id:
        raise ValueError("fig_id is required to uniquely identify the figure.")

    # Update the figure layout if needed.
    widget.update_layout(
        autosize=False,
        margin=dict(l=50, r=50, t=50, b=50)
    )

    # Serialize the figure to JSON.
    fig_json = widget.to_json()

    # Return a dictionary with the figure data and table HTML.
    return {
        'fig_id': fig_id,
        'fig_json': fig_json,
        'title': title,
        'table_html': table_html  # Add the table HTML to the dictionary.
    }


def fix_colorbar(figure, facet, dtick=None):
    if facet is True:
        figure.data[0].colorbar.thicknessmode="pixels"
        figure.data[0].colorbar.lenmode="pixels"
        figure.data[0].colorbar.len=150
        figure.data[0].colorbar.yanchor="top"
        figure.data[0].colorbar.y=1
        if dtick:
            figure.data[0].colorbar.dtick=dtick
        figure.data[1].colorbar.thicknessmode="pixels"
        figure.data[1].colorbar.lenmode="pixels"
        figure.data[1].colorbar.len=150
        figure.data[1].colorbar.yanchor="top"
        figure.data[1].colorbar.y=1
        figure.data[1].colorbar.dtick=dtick
    
    else:
        figure.data[0].colorbar.thicknessmode="pixels"
        figure.data[0].colorbar.lenmode="pixels"
        figure.data[0].colorbar.len=150
        figure.data[0].colorbar.yanchor="top"
        figure.data[0].colorbar.y=1
        if dtick:
            figure.data[0].colorbar.dtick=dtick
    
    return figure

def create_figure_json(widget, title, fig_id):
    """
    Converts a Plotly figure to JSON data for deferred rendering.

    Parameters:
    - widget (go.Figure): The Plotly figure.
    - title (str): The title for the figure.
    - fig_id (str): A unique identifier for the figure.

    Returns:
    - dict: A dictionary containing the figure data and metadata.
    """
    # Ensure the figure has a unique ID.
    if not fig_id:
        raise ValueError("fig_id is required to uniquely identify the figure.")

    # Update the figure layout if needed.
    widget.update_layout(
        autosize=True,
        margin=dict(l=50, r=50, t=50, b=50)
    )

    # Serialize the figure to JSON.
    fig_json = widget.to_json()

    # Return a dictionary with the figure data.
    return {
        'fig_id': fig_id,
        'fig_json': fig_json,
        'title': title
    }

# Write a text HTML description of the sunburst chart that contains the genus and subtype of each sample. Describe that and then the functionalities of the sunburst chart.
sunburst_help_text_html = """
<h3>Sample Composition</h3>
<p>The sunburst chart shows the composition of the collection. Starting from the center, the chart is divided into segments representing the species present in the samples. Each species segment is further divided into subsegments representing the subtypes of the species.
Hover over each segment to view the species and subtype of the sample.</p>
<h3>Interactivity</h3>
<ul>
    <li>Click on a segment to zoom in and view the subtypes of the selected species.</li>
    <li>Click on the center of the chart to zoom out and return to the previous view.</li>
</ul>
"""

fastqc_help_text_html = """
<div>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        The <strong>FastQC summary table</strong> provides an overview of key quality metrics for each of your sequencing samples. 
        This table is designed to help you quickly assess the quality of your raw sequencing data before proceeding to downstream analyses.
    </p>
    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Table Columns:</h4>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>Total Sequences (R1 and R2):</strong> The total number of sequences (reads) for each read direction—forward (R1) and reverse (R2).</li>
        <li><strong>Sequences Flagged as Poor Quality (R1 and R2):</strong> The number and percentage of sequences flagged as poor quality by FastQC, displayed as <code>Number (Percentage%)</code>.</li>
        <li><strong>Sequence Length (R1 and R2):</strong> The length of the sequences for each read direction. Consistent lengths indicate proper sequencing.</li>
        <li><strong>%GC Content (R1 and R2):</strong> The percentage of guanine (G) and cytosine (C) nucleotides in the sequences, displayed as <code>Actual% (expected: Expected%)</code>.</li>
        <li><strong>Species:</strong> The species name associated with each sample.</li>
    </ul>
    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Color Coding in %GC Columns:</h4>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        To assist in data interpretation, the <strong>%GC</strong> content cells are color-coded:
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><span style="background-color: #d4edda; padding: 3px;">Light Green:</span> The actual %GC is within <strong>6%</strong> of the expected value, suggesting consistency with the species' genomic characteristics.</li>
        <li><span style="background-color: #f8d7da; padding: 3px;">Light Red:</span> The actual %GC differs from the expected value by more than <strong>6%</strong>. This does not necessarily indicate a problem but suggests further investigation may be needed.</li>
    </ul>
    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Interpreting the Table:</h4>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        Use the table to evaluate sequencing quality:
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>Low Poor Quality Sequences:</strong> Indicates good sequencing quality.</li>
        <li><strong>High Poor Quality Sequences:</strong> May require investigation into sequencing conditions or sample integrity.</li>
        <li><strong>%GC Deviations (Red):</strong> Could suggest contamination, sequencing bias, or natural genomic variation.</li>
    </ul>
    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Action Steps for Unusual Results:</h4>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li>Review sample labeling and preparation procedures.</li>
        <li>Examine FastQC detailed reports for additional metrics.</li>
        <li>Verify sample purity and check for contamination sources.</li>
        <li>Consult with the sequencing facility if necessary.</li>
    </ul>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        <strong>Note:</strong> A red %GC cell is a prompt to investigate further and does not necessarily indicate an issue. Natural variation or technical factors might account for the deviation.
    </p>
    <p>
    The FastQC summary table is a valuable tool for quickly assessing the quality and integrity of your sequencing data. While the color-coded indicators provide visual cues for potential issues, they are meant to guide further investigation rather than serve as definitive judgments.
Remember: A red-colored cell in the %GC column is a prompt to look closer—it does not necessarily mean there is a problem. It is possible that the observed %GC variation is due to acceptable natural variation or technical factors that can be accounted for in your analysis.
    </p>
</div>

"""

assembly_qc_help_text_html = """
<div>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        The contig assembly quality plots provide insights into the quality and characteristics of the assembled genomes for each species in your dataset. These plots help you assess the assembly process, identify potential issues, and compare assembly metrics across samples and species.
    </p>
    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">1. Cumulative Contig Length Plots</h4>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        For each species, a cumulative contig length plot is generated, displaying the cumulative sum of contig lengths sorted in descending order for each sample.
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>X-Axis (Contig Index):</strong> Represents the index of contigs sorted from largest to smallest.</li>
        <li><strong>Y-Axis (Cumulative Contig Length):</strong> Shows the cumulative length of contigs up to that index.</li>
    </ul>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        The plot highlights the <strong>N50</strong> value for each sample:
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><span style="color: red;"><strong>Red Segment:</strong></span> Indicates the portion of the cumulative length up to the N50 contig. This is the contig length at which half of the total assembly length is reached.</li>
        <li><span style="color: #6c757d;"><strong>Remaining Line:</strong></span> Represents the rest of the contigs beyond the N50.</li>
    </ul>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        <strong>Interpreting the Plot:</strong>
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li>A steeper curve suggests fewer contigs are needed to reach the total assembly length, indicating a more contiguous assembly.</li>
        <li>A longer red segment (up to N50) means larger contigs contribute significantly to the assembly.</li>
        <li>Comparing samples can reveal variations in assembly quality across different samples of the same species.</li>
    </ul>
    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">2. N50 Contig Length Distribution Box Plot</h4>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        This plot displays the distribution of N50 contig lengths across all samples, grouped by species.
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>X-Axis (Species):</strong> Different species in your dataset.</li>
        <li><strong>Y-Axis (N50 Contig Length):</strong> The N50 values for the samples.</li>
    </ul>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        <strong>Interpreting the Plot:</strong>
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li>Higher N50 values indicate better assembly continuity.</li>
        <li>The spread of data points shows variability between samples.</li>
        <li>Outliers can identify samples with unusually high or low N50 values, which may need further investigation.</li>
    </ul>
    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">3. N50 Contig Coverage Box Plot</h4>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        This plot illustrates the distribution of coverage at the N50 contig across all samples, grouped by species.
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>X-Axis (Species):</strong> Different species in your dataset.</li>
        <li><strong>Y-Axis (N50 Contig Coverage):</strong> Coverage values at the N50 contig for the samples.</li>
    </ul>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        <strong>Interpreting the Plot:</strong>
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li>Consistent coverage values suggest uniform sequencing depth.</li>
        <li>Variations in coverage may indicate sequencing biases or issues with library preparation.</li>
        <li>Outliers may highlight samples requiring further quality assessment.</li>
    </ul>
    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Important Notes</h4>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        <strong>N50 Value:</strong> The N50 is a commonly used metric in genomics to assess the quality of genome assemblies. It is the contig length such that 50% of the total assembly length is contained in contigs equal to or larger than this length.
    </p>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        <strong>Interpreting High N50 and Coverage:</strong> While higher N50 and consistent coverage are generally positive indicators, they should be interpreted in the context of the organism's genome and the sequencing technology used.
    </p>
    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Action Steps if Issues are Observed</h4>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li>Review sequencing data quality and consider re-sequencing if necessary.</li>
        <li>Check for contamination or mixed samples if unexpected results are observed.</li>
        <li>Adjust assembly parameters or use alternative assembly tools for improvement.</li>
    </ul>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        If you have questions or need assistance in interpreting these plots, please consult with a bioinformatics specialist.
    </p>
</div>
"""

resistome_help_text_html = """
<div>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        This section provides insights into the distribution of resistance genes across different subtypes within a species. The interactive polar chart and the accompanying table are designed to help you explore and interpret the data effectively.
    </p>

    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">1. Understanding the Subtype Polar Chart</h4>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        The Subtype Polar Chart is a circular barplot that visualizes the average number of resistance genes per sample for each subtype within a species.
    </p>

    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>Bars Represent Resistance Genes:</strong> Each bar corresponds to a subtype and represents the average number of resistance genes per sample for that subtype. The height of the bar indicates the average gene count.</li>
        <li><strong>Ordering of Bars:</strong> The bars are ordered based on the average number of resistance genes per sample, with subtypes having higher averages appearing first.</li>
        <li><strong>Color-Coded Antibiotics Resistance Classes:</strong> Different resistance types are color-coded for easy identification. This allows you to see which resistances are most prevalent in each subtype.</li>
        <li><strong>Center Line (Circular Line Plot):</strong> The center line that loops around the chart represents the number of samples belonging to each subtype. The position and size of the markers on this line indicate the sample count.</li>
        <li><strong>Subtype Labels:</strong> Each subtype is labeled around the chart, making it easy to identify which bar corresponds to which subtype.</li>
    </ul>

    <h5 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Interactivity:</h5>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        You can interact with the chart to get more detailed information:
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>Hover Over Bars:</strong> Hovering over any bar segment will display a tooltip with the following information:
            <ul>
                <li><strong>Subtype:</strong> The subtype associated with the bar.</li>
                <li><strong>Resistance:</strong> The resistance type represented by that segment.</li>
                <li><strong>Genes:</strong> The total number of genes contributing to that segment.</li>
            </ul>
        </li>
        <li><strong>Hover Over Center Line Points:</strong> Hovering over any point on the center line will display:
            <ul>
                <li><strong>Subtype:</strong> The subtype at that position.</li>
                <li><strong>Samples:</strong> The number of samples belonging to that subtype.</li>
            </ul>
        </li>
    </ul>

    <h5 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Interpreting the Chart:</h5>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        The chart provides a visual summary of resistance gene distribution across subtypes:
    </p>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>High Bars:</strong> Subtypes with higher bars have a greater average number of resistance genes per sample.</li>
        <li><strong>Dominant Resistances:</strong> The color segments within the bars indicate which resistance types are most common in each subtype.</li>
        <li><strong>Sample Distribution:</strong> The center line helps you understand how many samples belong to each subtype, which is important for assessing the representativeness of the data.</li>
    </ul>

    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">2. Exploring the Resistance Genes Table</h4>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        Next to the polar chart, there is a table that provides detailed information about the resistance genes identified in each sample.
    </p>

    <h5 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Table Structure:</h5>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        The table is organized hierarchically to facilitate easy navigation and interpretation:
    </p>

    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>SUBTYPE:</strong> The highest level of grouping in the table. Subtypes are displayed with bold text and larger font size for emphasis.</li>
        <li><strong>SAMPLE:</strong> Under each subtype, samples are listed. The sample names are bolded and slightly smaller in font size than the subtypes.</li>
        <li><strong>RESISTANCE:</strong> For each sample, the resistance types identified are listed in bold text.</li>
        <li><strong>GENES:</strong> Under each resistance type, the specific genes detected are displayed. Hovering over a gene will show the percentage identity (%IDENTITY) to reference sequences.</li>
    </ul>

    <h5 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Interactivity:</h5>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        The table is interactive:
    </p>

    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>Hover Over Genes:</strong> Hovering over any gene will display a tooltip with the %IDENTITY, indicating how closely the gene matches known resistance genes.</li>
    </ul>

    <h5 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">Using the Table:</h5>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        The table allows you to:
    </p>

    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>Identify Resistance Profiles:</strong> See which resistance genes are present in each sample and subtype.</li>
        <li><strong>Assess Gene Diversity:</strong> Observe the variety of genes contributing to resistance within samples and subtypes.</li>
        <li><strong>Compare Subtypes and Samples:</strong> Quickly compare resistance profiles across different subtypes and samples.</li>
    </ul>

    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">3. Tips for Interpretation</h4>
    <ul style="font-family: Arial, sans-serif; font-size: 16px;">
        <li><strong>Correlate Chart and Table:</strong> Use the polar chart to identify subtypes of interest and refer to the table for detailed gene information.</li>
        <li><strong>Look for Patterns:</strong> Observe if certain resistance types are prevalent in specific subtypes or if there are unique resistance profiles.</li>
        <li><strong>Sample Representativeness:</strong> Keep in mind the number of samples per subtype when interpreting the average gene counts.</li>
    </ul>

    <h4 style="font-family: 'Arial Narrow', Arial, sans-serif; font-weight: bold;">4. Additional Information</h4>
    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        The polar chart and table are dynamically generated based on your data. The visualization adjusts to highlight the most relevant information for your dataset.
    </p>

    <p style="font-family: Arial, sans-serif; font-size: 16px;">
        If you have questions or need further assistance in interpreting these visualizations, please consult the documentation or reach out to a bioinformatics specialist.
    </p>
</div>
"""
