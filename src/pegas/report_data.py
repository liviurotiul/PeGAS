import warnings
from datetime import datetime

import matplotlib
import pandas as pd
def load_and_preprocess_data():
    """
    Loads the results.csv file and preprocesses the DataFrame.

    Returns:
    - df (DataFrame): Preprocessed DataFrame.
    """
    warnings.filterwarnings("ignore")
    matplotlib.use('Agg')

    # Read results.csv.
    df = pd.read_csv("dataframe/results.csv", dtype={'SAMPLE': str})

    # Set the columns to the correct data types.
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

    df["SEQUENCE"] = df["SEQUENCE"].str.replace("contig", "")

    # Convert to int.
    df["CONTIG"] = df["SEQUENCE"].astype(int)

    return df.drop(columns=["SEQUENCE"])


def generate_general_info_table(df, input_dir, output_dir):
    """
    Generates an HTML table containing general information about the run.

    Parameters:
    - df (DataFrame): The DataFrame containing the data.
    - input_dir (str): The path of the input directory.
    - output_dir (str): The path of the output directory.

    Returns:
    - str: The HTML string of the table.
    """
    # Calculate the required values.
    report_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    num_samples = len(df['SAMPLE'].unique())
    num_genes = len(df['GENE'])

    # Build the table data.
    table_data = [
        ['Date', report_date],
        ['Number of Samples', num_samples],
        ['Number of Unique Genes', num_genes],
        ['Input Directory', input_dir],
        ['Output Directory', output_dir]
    ]

    # Initialize the HTML table.
    html = '''
    <table style="margin: 0 auto; font-family: Arial, sans-serif; border-collapse: collapse;">
    '''

    # Build the table rows.
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
