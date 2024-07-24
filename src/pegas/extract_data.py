import pandas as pd
import os
from app.paths import *
import re
from tqdm import tqdm
from app.utils import *
# length_(.*)_cov


def generate_tables():

    tsv_files, paths = find_extension(DATABASE_PATH, '.tsv')

    dfs = []

    for file, path in zip(tsv_files, paths):
        if not file.startswith("abricate_"):
            continue

        df = pd.read_csv(os.path.join(path, file), sep='\t')
        df = df.drop(["#FILE", "SEQUENCE", "START", "END", "STRAND", "COVERAGE", "COVERAGE_MAP", "GAPS", "DATABASE"], axis=1)
        df["FILE"] = path.split('/')[-1]
        dfs.append(df)

    df_genes = pd.concat(dfs)
    df_genes = df_genes.reset_index()

    accession_product = {}
    accession_resistance = {}

    def extract_abricate_data(source):

        def filter_list_an(tsv):
            for (acc, pro, res) in zip(tsv['ACCESSION'].to_list(), tsv['PRODUCT'].to_list(), tsv['RESISTANCE'].to_list()):
                accession_product[acc] = pro

                try:
                    res = res.split(';')
                    for med in res:
                        if accession_resistance.get(med, False):
                            accession_resistance[acc] += med
                        else:
                            accession_resistance[acc] = [med]
                except:
                    continue

            an_list = [*set(tsv['ACCESSION'].to_list()), ]
            an_list = [x for x in an_list if str(x) != 'nan']

            return an_list

        accesions_ncbi, accesions_plasmids, accesions_vfdb = [], [], []
        products_ncbi, products_plasmids, products_vfdb = [], [], []

        accesions_ncbi += filter_list_an(pd.read_csv(f"{source}/abricate_ncbi.tsv", sep='\t'))

        accesions_plasmids += filter_list_an(pd.read_csv(f"{source}/abricate_plasmidfinder.tsv", sep='\t'))

        accesions_vfdb += filter_list_an(pd.read_csv(f"{source}/abricate_vfdb.tsv", sep='\t'))

        return accesions_ncbi, accesions_plasmids, accesions_vfdb

    abricate = None
    species = os.listdir(DATABASE_PATH)
    try:
        species.remove('.gitignore')
    except:
        pass

    df_species = {}
    df_hierarchical_full = {}
    df_hierarchical_short = {}
    df_hierarchical_subtype = {}

    samples_samples = []
    samples_subtypes = []
    samples_ncbi = []
    samples_plasmids = []
    samples_vfdb = []
    samples_species = []
    samples_max_contig_len = []
    samples_max_contig_cov = []

    for sp in species:
        subtypes = os.listdir(f"{DATABASE_PATH}/{sp}")

        species_max_contig_len = []
        species_max_contig_cov = []
        species_samples = []
        species_abricate = []

        for subtype in subtypes:
            samples = os.listdir(f"{DATABASE_PATH}/{sp}/{subtype}")
            try:
                species.remove('.gitignore')
            except:
                pass

            abricate_hits = 0

            for sample in samples:

                if not species_mapping.get(sp, False):
                    tqdm.write(f"{sp} is not in the internal species dictionary, either add it or request an update to the dictionary on the github page or author's personal email")
                    continue
                    
                lengths = read_contig(f"{DATABASE_PATH}/{sp}/{subtype}/{sample}/shovill/contigs.fasta")
                species_max_contig_len.append(lengths[0][0])
                species_max_contig_cov.append(lengths[0][1])
                ncbi, plasmids, vfdb = extract_abricate_data(f"{DATABASE_PATH}/{sp}/{subtype}/{sample}")
                abricate = ncbi + plasmids + vfdb
                
                df_hierarchical_full[(species_mapping[sp], subtype, sample)] = {'abricate': ', '.join(abricate)}
                species_abricate += abricate
                no_ncbi, no_plasmids, no_vfdb = len(ncbi), len(plasmids), len(vfdb)
                abricate_hits = abricate_hits + no_ncbi + no_plasmids + no_vfdb

                df_hierarchical_short[(species_mapping[sp], subtype, sample)] = {
                    'abricate-ncbi': no_ncbi,
                    'abricate-plasmid-finder': no_plasmids,
                    'arbicate-vfdb': no_vfdb
                }

                # samples_abricate.append(abricate)
                samples_ncbi.append(ncbi)
                samples_plasmids.append(plasmids)
                samples_vfdb.append(vfdb)
                samples_species.append(species_mapping[sp])
                samples_subtypes.append(subtype)
                samples_samples.append(sample)
                samples_max_contig_len.append(lengths[0][0])
                samples_max_contig_cov.append(lengths[0][1])

            df_hierarchical_subtype[(species_mapping[sp], subtype)] = {
                'abricate hits': abricate_hits
            }

        df_species[species_mapping[sp]] = {
            'subtypes': subtypes,
            'contigs': species_max_contig_len,
            'samples': species_samples,
            'abricate': species_abricate
        }

    df_samples = pd.DataFrame.from_dict({
        'name': samples_samples,
        'subtype': samples_subtypes,
        'species': samples_species,
        'maximum_length_contig': samples_max_contig_len,
        'maximum_length_contig_coverage': samples_max_contig_cov,
        'ncbi_resistance': samples_ncbi,
        'plasmids': samples_plasmids,
        'virulence_factors': samples_vfdb
    })

    df_species = pd.DataFrame.from_dict(df_species)
    df_hierarchical_short = pd.DataFrame.from_dict(df_hierarchical_short)
    df_hierarchical_full = pd.DataFrame.from_dict(df_hierarchical_full)
    df_hierarchical_subtype = pd.DataFrame.from_dict(df_hierarchical_subtype)

    return (df_samples, df_species, df_hierarchical_short, df_hierarchical_full, df_hierarchical_subtype, accession_product, accession_resistance, df_genes)


def read_contig(contig_path):

    file_contents = None

    with open(contig_path) as f:
        file_contents = f.read()

    lengths = [x.split('_') for x in re.findall('length_(.*)_pilon', file_contents)]
    lengths = [(int(x), float(z)) for (x,y,z) in lengths]

    return lengths

