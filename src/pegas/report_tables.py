import os
import re

import pandas as pd
from bs4 import BeautifulSoup as bs
from tqdm import tqdm

from utils import create_html_widget


def create_annotation_tables(df, species_dict):
    """
    Creates three HTML tables summarizing annotations from NCBI, VFDB, and PlasmidFinder,
    with custom CSS to set the table width to 100%.

    Parameters:
    - df (DataFrame): The DataFrame containing the data.
    - species_dict (dict): Mapping from species codes to common names.

    Returns:
    - tables (dict): A dictionary with keys 'Resistance', 'Virulence', 'Plasmid' and values as HTML code for the tables.
    """
    # Map species codes to common names in the DataFrame.
    df['SPECIES_NAME'] = df['SPECIES'].map(species_dict)

    # Define predictors and their corresponding table names.
    predictors = {
        'NCBI': 'Resistance',
        'VFDB': 'Virulence',
        'PlasmidFinder': 'Plasmid'
    }

    # Initialize a dictionary to hold the HTML tables.
    tables = {}

    for predictor, table_name in predictors.items():
        # Filter the DataFrame for the current predictor.
        df_predictor = df[df['PREDICTION_SOURCE'] == predictor]

        # Check if the filtered DataFrame is not empty.
        if not df_predictor.empty:
            # Select relevant columns (exclude unnecessary ones).
            df_table = df_predictor.drop(columns=["CONTIG_LENGTH", "CONTIG_COVERAGE", "CONTIG_NUMBER"], errors='ignore')

            # Remove duplicates if any.
            df_table = df_table.drop_duplicates()

            # Sort the table by SAMPLE for better readability.
            df_table = df_table.sort_values(by=['SAMPLE'])

            # Generate HTML code for the table.
            table_html = df_table.to_html(
                classes='table table-striped',  # Might need to comment this; it shows search and pagination twice.
                index=False,
                border=0,
                table_id=f'{table_name}_table'
            )

            # Parse the HTML to modify the table tag.
            soup = bs(table_html, 'html.parser')

            # Set the style attribute to make the table width 100%.
            table_tag = soup.find('table')
            if table_tag:
                table_tag['style'] = 'width: 100%;'

            # Convert back to string.
            table_html = str(soup)

            # Store the HTML table in the dictionary.
            tables[table_name] = create_html_widget(table_html, f"{table_name}")
        else:
            # If there is no data for this predictor, return a message.
            tables[table_name] = f'<p>No data available for {table_name}</p>'

    return tables


def create_fastqc_table(df, species_dict, gc_content_dict):
    """
    Creates an HTML table summarizing FastQC outputs for all samples.

    Parameters:
    - df (DataFrame): DataFrame containing sample information.
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

    # Find all FastQC HTML files.
    html_files = []
    for root, dirs, files in os.walk("fastqc"):
        for file in files:
            if file.endswith(".html"):
                html_files.append(os.path.join(root, file))

    # Define the strings we're looking for in the FastQC reports.
    search_strings = ["Total Sequences", "Sequences flagged as poor quality", "Sequence length", "%GC"]

    # Function to extract sample name.
    def get_core_sample_name(filename):
        """Extracts the core sample name by removing _R1 or _R2 and other suffixes."""
        return os.path.basename(filename).replace("_R1", "").replace("_R2", "").replace(".fastq.gz", "").replace("_fastqc.html", "")

    # Function to extract read direction from file path.
    def get_read_direction_from_file(file_path):
        base_name = os.path.basename(file_path)
        match = re.search(r'_R([12])', base_name)
        if match:
            return 'R' + match.group(1)
        else:
            # Try to infer from the directory structure if not in filename.
            if 'R1' in file_path or 'read1' in file_path.lower():
                return 'R1'
            elif 'R2' in file_path or 'read2' in file_path.lower():
                return 'R2'
            else:
                return None

    # Build a set of unique sample names.
    sample_names = set()
    for file_path in html_files:
        sample_name = get_core_sample_name(file_path)
        sample_names.add(sample_name)

    # Create the DataFrame with the corrected index.
    df_qc = pd.DataFrame(
        index=sorted(sample_names),
        columns=[f"{x} R1" for x in search_strings] + [f"{x} R2" for x in search_strings] + ['SPECIES']
    )

    # Loop over each file path and extract data.
    for file_path in html_files:
        # Load the HTML file into BeautifulSoup.
        with open(file_path, "r", encoding='utf-8') as f:
            soup = bs(f, "html.parser")

        sample_name = get_core_sample_name(file_path)
        read_direction = get_read_direction_from_file(file_path)

        if read_direction is None:
            tqdm.write(f"Warning: Could not determine read direction for file '{file_path}'. Skipping.")
            continue  # Skip files without read direction.

        # Loop over each search string.
        for search_string in search_strings:
            # Find the <td> element with the search string as its content.
            td = soup.find("td", string=search_string)

            # If the <td> element exists, extract the content of the next <td> element.
            if td is not None:
                next_td = td.find_next_sibling("td")
                if next_td is not None:
                    column_name = f"{search_string} {read_direction}"
                    df_qc.loc[sample_name, column_name] = next_td.get_text().strip()
                else:
                    tqdm.write(f"Warning: Next <td> element not found for '{search_string}' in file '{file_path}'.")
            else:
                tqdm.write(f"Warning: '{search_string}' not found in file '{file_path}'.")

    # Assign 'SPECIES' to each sample in df_qc based on df_samples and species_dict.
    for sample_name in df_qc.index:
        sample_row = df_samples[df_samples["SAMPLE"] == sample_name]

        if not sample_row.empty:
            species_code = sample_row["SPECIES"].values[0]
            species_common_name = species_dict.get(species_code, 'Unknown')
            df_qc.loc[sample_name, "SPECIES"] = species_common_name
        else:
            tqdm.write(f"Warning: Sample '{sample_name}' not found in df_samples. Assigning 'Unknown' to SPECIES.")
            df_qc.loc[sample_name, "SPECIES"] = 'Unknown'

    # Convert relevant columns to numeric, handling non-numeric values.
    numeric_columns = [col for col in df_qc.columns if any(x in col for x in ["Total Sequences", "Sequences flagged as poor quality"])]
    for col in numeric_columns:
        df_qc[col] = pd.to_numeric(df_qc[col], errors='coerce')

    # Calculate percentages and format the 'Sequences flagged as poor quality' columns.
    for read_direction in ['R1', 'R2']:
        total_sequences_col = f"Total Sequences {read_direction}"
        poor_quality_col = f"Sequences flagged as poor quality {read_direction}"
        if total_sequences_col in df_qc.columns and poor_quality_col in df_qc.columns:
            # Avoid division by zero and handle missing data.
            total_sequences = df_qc[total_sequences_col].replace(0, pd.NA)
            poor_quality_sequences = df_qc[poor_quality_col]
            with pd.option_context('mode.use_inf_as_na', True):
                percentage = (poor_quality_sequences / total_sequences * 100).round(2)
            # Combine counts and percentages.
            df_qc[poor_quality_col] = poor_quality_sequences.fillna(0).astype(int).astype(str) + ' (' + percentage.fillna(0).astype(str) + '%)'

    # Compute %GC intervals and annotate with expected values.
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

    samples = df_qc.index.tolist()

    aseembly_qc = []
    for sample in samples:
        contigs_file = "results/" + sample + "/shovill/contigs.fa"
        with open(contigs_file, "r") as f:
            contigs_content = f.read()
        if ">" in contigs_content:
            aseembly_qc.append("Success")
        else:
            aseembly_qc.append("Fail")
    df_qc["Aseembly"] = aseembly_qc
    # Generate HTML table code.
    df_reads_table_code = df_qc.to_html(classes='table table-striped', table_id='full_qa_table', index=True)

    # Parse the HTML code using BeautifulSoup.
    soup = bs(df_reads_table_code, 'html.parser')

    # Style the table headers.
    for th in soup.find_all('th'):
        th['style'] = 'text-align: left'

    # Apply color coding to cells based on conditions.
    # Define a function to apply colors based on the difference in %GC.
    def apply_gc_color(td):
        text = td.get_text().strip()
        if "N/A" in text:
            td['style'] = 'background-color: #f8d7da;'  # Light red.
            return

        match = re.search(r'([\d\.]+)% \(expected: ([\d\.]+)%\)', text)
        if match:
            actual_gc = float(match.group(1))
            expected_gc = float(match.group(2))
            difference = abs(actual_gc - expected_gc)
            if difference > 6:
                color = '#f8d7da'  # Light red.
            else:
                color = '#d4edda'  # Light green.
            td['style'] = f'background-color: {color};'
        else:
            td['style'] = 'background-color: #ffffff;'

    for row in soup.find('table').find_all('tr'):
        for td in row.find_all('td'):
            # Assuming the %GC columns are at specific positions.
            for gc_col in [f"%GC R1", f"%GC R2"]:
                if gc_col in df_qc.columns:
                    apply_gc_color(td)

    # Convert the modified table back to HTML.
    df_reads_table_code = str(soup)

    # Wrap the table in a widget.
    fastqc_table_html = create_html_widget(df_reads_table_code, "FastQC output - table")

    return fastqc_table_html


def generate_table_html_for_species(sp, df):
    """
    Generates a grouped HTML table for the given species, organized by SUBTYPE, SAMPLE, and RESISTANCE.

    Parameters:
    - sp (str): The species identifier.
    - df (DataFrame): The DataFrame containing the data.

    Returns:
    - str: The HTML string of the grouped table.
    """
    # Filter for the specified species and drop rows with missing RESISTANCE.
    df_species = df[df['SPECIES'] == sp].dropna(subset=['RESISTANCE'])

    # Define the columns you want to display.
    columns_to_display = ['SUBTYPE', 'SAMPLE', 'RESISTANCE', 'GENE', '%IDENTITY']
    df_species = df_species[columns_to_display]

    # Sort the data for consistent grouping.
    df_species.sort_values(['SUBTYPE', 'SAMPLE', 'RESISTANCE'], inplace=True)

    # Initialize the HTML table structure with Bootstrap classes.
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

    # Group by SUBTYPE, SAMPLE, and RESISTANCE.
    grouped = df_species.groupby(['SUBTYPE', 'SAMPLE', 'RESISTANCE'])

    # Initialize variables to keep track of rowspan and previous values.
    prev_subtype = None
    prev_sample = None
    subtype_rowspan = {}
    sample_rowspan = {}

    # Compute rowspan for SUBTYPE and SAMPLE.
    for (subtype, sample), group in df_species.groupby(['SUBTYPE', 'SAMPLE']):
        subtype_rowspan[subtype] = subtype_rowspan.get(subtype, 0) + len(group['RESISTANCE'].unique())
        sample_rowspan[(subtype, sample)] = sample_rowspan.get((subtype, sample), 0) + len(group['RESISTANCE'].unique())

    # Iterate over groups and build table rows.
    for (subtype, sample, resistance), group in grouped:
        genes_html = ""
        for _, row in group.iterrows():
            genes_html += f"""
            <span style='display: inline-block; margin-right: 5px;' title='%IDENTITY: {row['%IDENTITY']}%'>
                {row['GENE']}
            </span>
            """

        # Start table row.
        table_html += "<tr>"

        # SUBTYPE cell.
        if subtype != prev_subtype:
            rowspan = subtype_rowspan[subtype]
            table_html += f"<td rowspan='{rowspan}' class='subtype-cell'>{subtype}</td>"
            prev_subtype = subtype
            prev_sample = None  # Reset prev_sample when subtype changes.

        else:
            # If not first occurrence, skip the subtype cell.
            pass

        # SAMPLE cell.
        if sample != prev_sample:
            rowspan = sample_rowspan[(subtype, sample)]
            table_html += f"<td rowspan='{rowspan}' class='sample-cell'>{sample}</td>"
            prev_sample = sample
        else:
            # If not first occurrence, skip the sample cell.
            pass

        # RESISTANCE cell.
        table_html += f"<td class='resistance-cell'>{resistance}</td>"

        # GENES cell.
        table_html += f"<td class='genes-cell'>{genes_html}</td>"

        table_html += "</tr>"

    # Close the table.
    table_html += """
        </tbody>
    </table>
    """

    # Add custom CSS for styling.
    custom_css = """
    <style>
    /* Table styling */
    .table {
        width: 100%;
        text-align: left;
        border-collapse: collapse;
    }
    .table td, .table th {
        padding: 8px;
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


def generate_virulence_factor_table(df, species_dict):
    """
    Generates HTML tables of virulence factors for each species.

    Parameters:
    - df (DataFrame): The DataFrame containing the data.
    - species_dict (dict): Dictionary mapping species codes to species names.

    Returns:
    - tables (dict): Dictionary with species names as keys and HTML code for tables as values.
    """
    # Filter the DataFrame for virulence factors.
    df_virulence = df[df['PREDICTION_SOURCE'] == 'VFDB']

    species = df_virulence['SPECIES'].unique()

    tables = {}

    for sp in species:
        species_name = species_dict.get(sp, sp)

        df_species = df_virulence[df_virulence['SPECIES'] == sp]
        df_species = df_species.drop(columns=['SPECIES', 'CONTIG_LENGTH', 'CONTIG_COVERAGE', 'PREDICTION_SOURCE', 'RESISTANCE', 'CONTIG_NUMBER'])
        df_species.reset_index(drop=True, inplace=True)

        # Generate HTML code for the table.
        html_table = df_species.to_html(
            classes='table table-striped table-bordered',
            table_id=f'virulence_table_{sp}',
            index=False,
            border=0,
            justify='left'
        )

        # Store in the dictionary.
        tables[species_name] = html_table

    return tables
