import pandas as pd
import os
import re

from tqdm import tqdm

# This script iterates over all the folders in results representing each sample and builds a 
# dataframe containing all the information regarding those sample and genes found there
# The dataframe is then saved as a csv file in results
def build_dataframe():
    try:
        # Iterate over all the folders in results
        for sample in tqdm(os.listdir("results")):
            # Get the abricate tsvs, mlst tsvs, and the contigs fasta in each folder without
            # using find_extension because we need to keep track of the paths
            abricate_ncbi, abricate_vfdb, abricate_plasmidfinder = None, None, None
            mlst = None

            # Iterate over all the files in each folder, read the abricate tsvs, mlst and contigs fasta
            # and add the information to the dataframe
            for file in os.listdir(f"results/{sample}"):

                # Read abricate_ncbi.tsv ànd extract the GENE, ACCESSSION and PRODUCT columns
                if file == "abricate_ncbi.tsv":
                    abricate_ncbi = pd.read_csv(f"results/{sample}/{file}", sep="\t")
                    abricate_ncbi = abricate_ncbi[["GENE", "ACCESSION", "PRODUCT", "%IDENTITY", "%COVERAGE", "RESISTANCE"]]

                    # Add the prediction source column
                    abricate_ncbi["PREDICTION_SOURCE"] = "NCBI"
                
                # Read abricate_vfdb.tsv and extract the GENE, ACCESSSION and PRODUCT columns
                if file == "abricate_vfdb.tsv":
                    abricate_vfdb = pd.read_csv(f"results/{sample}/{file}", sep="\t")
                    abricate_vfdb = abricate_vfdb[["GENE", "ACCESSION", "PRODUCT", "%IDENTITY", "%COVERAGE", "RESISTANCE"]]

                    # Add the prediction source column
                    abricate_vfdb["PREDICTION_SOURCE"] = "VFDB"
                
                # Read abricate_plasmidfinder.tsv and extract the GENE, ACCESSSION and PRODUCT columns
                if file == "abricate_plasmidfinder.tsv":
                    abricate_plasmidfinder = pd.read_csv(f"results/{sample}/{file}", sep="\t")
                    abricate_plasmidfinder = abricate_plasmidfinder[["GENE", "ACCESSION", "PRODUCT", "%IDENTITY", "%COVERAGE", "RESISTANCE"]]

                    # Add the prediction source column 
                    abricate_plasmidfinder["PREDICTION_SOURCE"] = "PlasmidFinder"
                
                # Read contigs.fa from the shovill folder and extract the the length of the
                # longest contig and the coverage that it has; Check for the shovill folder
                # The header looks like this:
                # >contig00001 len=129428 cov=64.2 corr=0 spades=NODE_1_length_129428_cov_64.218227_pilon
                if os.path.isdir(f"results/{sample}/shovill"):
                    
                    # Check for the contigs.fa file
                    if os.path.isfile(f"results/{sample}/shovill/contigs.fa"):
                        contigs = open(f"results/{sample}/shovill/contigs.fa", "r")
                        contigs = contigs.read()

                        # Get the first line of the file
                        contig_data = contigs.split("\n")[0]

                        # Get the length of the longest contig len=...
                        contig_length = int(re.findall(r"len=(\d+)", contig_data)[0])

                        # Get the coverage of the longest contig cov=...
                        contig_coverage = float(re.findall(r"cov=(\d+.\d+)", contig_data)[0])

                        # Number of contigs
                        contig_number = len(re.findall(r">", contigs))
                
                # Read mlst.tsv and extract thesecond and third columns
                # mlst.tsv has no header        
                if file == "mlst.tsv":
                    mlst = pd.read_csv(f"results/{sample}/{file}", sep="\t", header=None)

                    # Keep only the second and third columns by number
                    mlst = mlst.iloc[:, 1:3]

                    # Name the columns
                    mlst.columns = ["SPECIES", "SUBTYPE"]
            
            # Concatenate the abricate dataframes and add a sample name column
            df = pd.concat([abricate_ncbi, abricate_vfdb, abricate_plasmidfinder])
            df["SAMPLE"] = sample

            # Add the contig length and coverage columns
            df["CONTIG_LENGTH"] = contig_length
            df["CONTIG_COVERAGE"] = contig_coverage
            df["CONTIG_NUMBER"] = contig_number

            # Add the SPECIES and SUBTYPE columns
            df["SPECIES"] = mlst["SPECIES"][0]
            df["SUBTYPE"] = mlst["SUBTYPE"][0]

            # Add df to the results dataframe without using try except clause
            if "results" not in locals():
                results = df
            else:
                results = pd.concat([results, df])
            

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

        # Save the results dataframe as a csv file
        results.to_csv("dataframe/results.csv", index=False)

    except Exception as e:
        raise e