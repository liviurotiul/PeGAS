# Here you can find different useful functions for handling files and database operations
import os
import fnmatch
import shutil
import subprocess
import random

from tqdm import tqdm
from subprocess import call
from typing import List, Dict, Tuple


def generate_color():
    # Generate a random hue value between 0 and 360 (in degrees)
    h = random.randint(0, 359)
    # Set saturation and lightness values to fixed values
    s = 90
    l = 60
    # Convert HSL to RGB and format as a hexadecimal color code
    r, g, b = hsl_to_rgb(h, s, l)
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def hsl_to_rgb(h, s, l):
    # Convert HSL values to floating-point numbers between 0 and 1
    h = h / 360
    s = s / 100
    l = l / 100

    # Calculate intermediate values
    if s == 0:
        r = g = b = l
    else:
        q = l * (1 + s) if l < 0.5 else l + s - l * s
        p = 2 * l - q
        r = hue_to_rgb(p, q, h + 1/3)
        g = hue_to_rgb(p, q, h)
        b = hue_to_rgb(p, q, h - 1/3)

    # Convert RGB values to integers between 0 and 255
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
    "klebsiella": "Klebsiella species",
    "vcholerae_2": "Vibrio cholerae (El Tor biotype)",
    "vibrio": "Vibrio species",
    "hparasuis": "Haemophilus parasuis",
    "ssuis": "Streptococcus suis",
    "pmultocida": "Pasteurella multocida",
    "spneumoniae": "Streptococcus pneumoniae",
    "diphtheria_3": "Corynebacterium diphtheriae (biovar mitis)",
    "cronobacter": "Cronobacter species",
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

gc_content_dict = {
    "Escherichia coli": "50.8",
    "Streptomyces chromogenes": "72.1",
    "Burkholderia pseudomallei": "68.8",
    "Campylobacter (non-Jejuni)": "31.5",
    "Mannheimia haemolytica": "41.5",
    "Streptococcus zooepidemicus": "41.6",
    "Salmonella enterica": "52.2",
    "Lactococcus lactis phage": "39.2",
    "Klebsiella pneumoniae": "57.2",
    "Enterococcus faecalis": "37.9",
    "Enterobacter cloacae": "55.1",
    "Acinetobacter baumannii": "39.0",
    "Streptococcus agalactiae": "35.8",
    "Enterococcus faecium": "37.9",
    "Streptococcus suis": "41.8",
    "Streptococcus pneumoniae": "39.7",
    "Clostridium perfringens": "28.6"
}

def dilute_hex_color(hex_color, factor):
    # Convert hex color string to RGB tuple
    hex_color = hex_color.lstrip('#')
    rgb_color = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

    # Apply lightening factor to each RGB component
    r = min(int(rgb_color[0] + 255 * factor), 255)
    g = min(int(rgb_color[1] + 255 * factor), 255)
    b = min(int(rgb_color[2] + 255 * factor), 255)

    # Convert RGB tuple to hex color string
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)


def find_file(filename, search_path):
    """Recursively searches for a file in a given directory and its subdirectories."""

    # Iterate through all items in the search path
    for root, dirnames, filenames in os.walk(search_path):

        # Check if the file exists in the current directory
        if filename in filenames:
            return os.path.join(root, filename)

    # If the file was not found in any of the directories, return None
    return None


def create_html_card(widget, title, width="8", min_height="300px"):

    collapse_id = "collapse_" + str(random.randint(0, 10**10))

    card_template = f"""
    <div class="row">
        <div class="col-xl-12">
            <div class="accordion" id="accordionExample">
                <div class="accordion-item">
                    <h2 class="accordion-header" id="headingOne">
                        <button class="accordion-button" type="button" data-bs-toggle="collapse" data-bs-target="#{collapse_id}" aria-expanded="true">
                            {title}
                        </button>
                    </h2>
                    <div id="{collapse_id}" class="accordion-collapse collapse show" data-bs-parent="#accordionExample">
                        <div class="accordion-body">
                            <div class="card shadow mb-4 h-90">
                                <!-- Card Body -->
                                <div class="card-body">
                                {widget}
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>                        
        </div>
    </div>
    """
    return card_template

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
            size=15
        ),
        margin=dict(t=20, b=20, l=30, r=30),
    )

    return pie_chart_figure


def connection_graph_string_format(hover_str, df, accession_product):

    ncbi, plasmids, vfdb = [], [], []
    ret_str = ""
    for item in hover_str:
        if item in df['ncbi_resistance'].sum():
            ncbi.append(accession_product[item])
        if item in df['plasmids'].sum():
            plasmids.append(accession_product[item])
        if item in df['virulence_factors'].sum():
            vfdb.append(accession_product[item])

    if len(ncbi):
        ret_str = ret_str + "<b>resistance factors</b>: <br> " + str(ncbi)
    if len(plasmids):
        ret_str = ret_str + "<br><b>plasmids</b>:<br> " + str(plasmids)
    if len(ncbi):
        ret_str = ret_str + "<br><b>virulence_factors</b>: <br> " + str(vfdb)

    ret_str = ret_str.translate(str.maketrans({'[': '', ']': '', '\'': '', ',': '<br>'}))

    return ret_str

def find_html_files(path):
    html_files = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".html"):
                html_files.append(os.path.join(root, file))
    return html_files


def create_html_element(widget, title, width="8", min_height="300px"):

    # If /report doesn't exist, create it
    if not os.path.exists("report"):
        os.system("mkdir report")

    # If /report/figures/ doesn't exist, create it
    if not os.path.exists("report/figures"):
        os.system("mkdir report/figures")

    widget.write_html(f"report/figures/{title}.html", full_html=False, include_plotlyjs='cdn')
    f = open(f"report/figures/{title}.html")
    code = create_html_widget(f.read(), title, width, min_height)
    return code


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


def find_file(filename, search_path):
    """Recursively searches for a file in a given directory and its subdirectories."""

    # Iterate through all items in the search path
    for root, dirnames, filenames in os.walk(search_path):

        # Check if the file exists in the current directory
        if filename in filenames:
            return os.path.join(root, filename)

    # If the file was not found in any of the directories, return None
    return None

