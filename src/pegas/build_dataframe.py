import pandas as pd
import os
import re

from tqdm import tqdm

# This script iterates over all folders in results representing each sample and builds a DataFrame containing information about those samples and genes.
# The DataFrame is then saved as a CSV file in results.
def build_dataframe():
    try:
        # Iterate over all the folders in results.
        for sample in tqdm(os.listdir("results")):
            # Get the abricate TSVs, MLST TSVs, and contigs FASTA in each folder without using find_extension because we need to keep track of the paths.
            abricate_ncbi, abricate_vfdb, abricate_plasmidfinder = None, None, None
            mlst = None
            n50_value = 0
            n50_coverage = 0.0
            contig_number = 0

            # Read contigs.fa once per sample and extract lengths/coverage.
            contigs_file = f"results/{sample}/shovill/contigs.fa"
            if os.path.isfile(contigs_file):
                with open(contigs_file, "r") as f:
                    contigs_content = f.read()

                contig_entries = contigs_content.strip().split(">")[1:]

                contig_lengths = []
                contig_coverages = []
                total_assembly_length = 0

                for contig in contig_entries:
                    lines = contig.strip().split("\n")
                    header = lines[0]
                    sequence = ''.join(lines[1:])

                    len_match = re.search(r'len=(\d+)', header)
                    cov_match = re.search(r'cov=([\d\.]+)', header)

                    if len_match:
                        length = int(len_match.group(1))
                        contig_lengths.append(length)
                        total_assembly_length += length
                    else:
                        continue

                    if cov_match:
                        coverage = float(cov_match.group(1))
                        contig_coverages.append(coverage)
                    else:
                        contig_coverages.append(0.0)

                if not contig_lengths:
                    tqdm.write(f"No contigs found in sample {sample}.")
                else:
                    sorted_contig_lengths = sorted(contig_lengths, reverse=True)
                    cumulative_length = 0
                    half_total_length = total_assembly_length / 2
                    n50_value = None

                    for length in sorted_contig_lengths:
                        cumulative_length += length
                        if cumulative_length >= half_total_length:
                            n50_value = length
                            break

                    if n50_value is not None:
                        n50_index = contig_lengths.index(n50_value)
                        n50_coverage = contig_coverages[n50_index]
                        contig_number = len(contig_lengths)

            # Iterate over all files in each folder, read the abricate TSVs, MLST, and contigs FASTA, and add the information to the DataFrame.
            for file in os.listdir(f"results/{sample}"):

                # Read abricate_ncbi.tsv and extract the GENE, ACCESSION, and PRODUCT columns.
                if file == "abricate_ncbi.tsv":
                    abricate_ncbi = pd.read_csv(f"results/{sample}/{file}", sep="\t")
                    abricate_ncbi = abricate_ncbi[["GENE", "ACCESSION", "PRODUCT", "%IDENTITY", "%COVERAGE", "RESISTANCE", "START", "END", "SEQUENCE", "STRAND"]]

                    # Add the prediction source column.
                    abricate_ncbi["PREDICTION_SOURCE"] = "NCBI"
                
                # Read abricate_vfdb.tsv and extract the GENE, ACCESSION, and PRODUCT columns.
                if file == "abricate_vfdb.tsv":
                    abricate_vfdb = pd.read_csv(f"results/{sample}/{file}", sep="\t")
                    abricate_vfdb = abricate_vfdb[["GENE", "ACCESSION", "PRODUCT", "%IDENTITY", "%COVERAGE", "RESISTANCE", "START", "END", "SEQUENCE", "STRAND"]]

                    # Add the prediction source column.
                    abricate_vfdb["PREDICTION_SOURCE"] = "VFDB"
                
                # Read abricate_plasmidfinder.tsv and extract the GENE, ACCESSION, and PRODUCT columns.
                if file == "abricate_plasmidfinder.tsv":
                    abricate_plasmidfinder = pd.read_csv(f"results/{sample}/{file}", sep="\t")
                    abricate_plasmidfinder = abricate_plasmidfinder[["GENE", "ACCESSION", "PRODUCT", "%IDENTITY", "%COVERAGE", "RESISTANCE", "START", "END", "SEQUENCE", "STRAND"]]

                    # Add the prediction source column.
                    abricate_plasmidfinder["PREDICTION_SOURCE"] = "PlasmidFinder"
                
                # Read mlst.tsv and extract the second and third columns.
                # mlst.tsv has no header.
                if file == "mlst.tsv":
                    mlst = pd.read_csv(f"results/{sample}/{file}", sep="\t", header=None)

                    # Keep only the second and third columns by number.
                    mlst = mlst.iloc[:, 1:3]

                    # Name the columns.
                    mlst.columns = ["SPECIES", "SUBTYPE"]
            
            # Concatenate the abricate DataFrames and add a sample name column.
            df = pd.concat([abricate_ncbi, abricate_vfdb, abricate_plasmidfinder])
            df["SAMPLE"] = sample

            # Add the contig length and coverage columns.
            df["CONTIG_LENGTH"] = n50_value
            df["CONTIG_COVERAGE"] = n50_coverage
            df["CONTIG_NUMBER"] = contig_number

            # Add the SPECIES and SUBTYPE columns.
            df["SPECIES"] = mlst["SPECIES"][0]
            df["SUBTYPE"] = mlst["SUBTYPE"][0]

            # Add df to the results DataFrame without using a try/except clause.
            if "results" not in locals():
                results = df
            else:
                results = pd.concat([results, df])
            

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

        # Save the results DataFrame as a CSV file.
        results.to_csv("dataframe/results.csv", index=False)

    except Exception as e:
        raise e
