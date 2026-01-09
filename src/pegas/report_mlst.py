import os

import pandas as pd
import plotly.graph_objects as go

from utils import dilute_hex_color


def create_sunburst_chart(df, species_dict, species_color_mapping, output_dir):
    """
    Creates a sunburst chart based on sample, subtype, and species data.

    Parameters:
    - df (DataFrame): The DataFrame containing the data.
    - species_dict (dict): Mapping from species codes to common names.
    - species_color_mapping (dict): Mapping from species common names to colors.

    Returns:
    - sunburst_figure (go.Figure): Plotly sunburst figure.
    """
    label = []
    parent = []
    value = []
    node_info = []
    colors = []

    # Find all mlst.tsv in the directory.
    mlst_df = pd.DataFrame()
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file == "mlst.tsv":
                mlst = pd.read_csv(os.path.join(root, file), sep="\t", header=None)

                # Keep only the second and third columns by number.
                mlst = mlst.iloc[:, 1:3]

                # Name the columns.
                mlst.columns = ["SPECIES", "SUBTYPE"]

                mlst["SAMPLE"] = os.path.basename(root)

                mlst_df = pd.concat([mlst_df, mlst], ignore_index=True)

    mlst_df["SPECIES_COMMON_NAME"] = mlst_df["SPECIES"].map(species_dict)

    seen_species, seen_subtype = set(), set()
    for _, row in mlst_df.iterrows():
        sp = row["SPECIES_COMMON_NAME"]
        st = row["SUBTYPE"]
        smp = row["SAMPLE"]
        st_label = f"ST{st}({sp})"

        if sp not in seen_species:
            label.append(sp)
            parent.append("")
            value.append(0)
            seen_species.add(sp)
        if st_label not in seen_subtype:
            label.append(st_label)
            parent.append(sp)
            value.append(0)
            seen_subtype.add(st_label)
        label.append(smp)
        parent.append(st_label)
        value.append(1)

    node_info, colors = [], []
    sample_list = mlst_df["SAMPLE"].nunique()

    mlst_df["SUBTYPE"] = (
        mlst_df["SUBTYPE"]
        .astype(str)
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)  # Handle 3.0 -> "3".
    )

    for par, itm in zip(parent, label):
        if par == "":  # Species node.
            n = (mlst_df["SPECIES_COMMON_NAME"] == itm).sum()
            node_info.append(
                f"Number of isolates: {n} <br> {round(n / sample_list * 100, 2)}% of total isolates"
            )
            colors.append(species_color_mapping.get(itm, "#999999"))
        elif "(" in par:  # Sample node.
            sp = par.split("(", 1)[1].rstrip(")")
            node_info.append("")
            colors.append(dilute_hex_color(species_color_mapping.get(sp, "#999999"), 0.1))
        else:  # Subtype node.
            sp = par
            st = str(itm.split("(", 1)[0].replace("ST", "")).strip()
            mask = (mlst_df["SPECIES_COMMON_NAME"] == sp) & (mlst_df["SUBTYPE"] == st)
            n = mask.sum()
            node_info.append(
                f"Number of isolates: {n} <br> {round(n / sample_list * 100, 2)}% of total isolates"
            )
            colors.append(dilute_hex_color(species_color_mapping.get(sp, "#999999"), 0.05))

    sunburst_data = {
        "label": label,
        "parent": parent,
        "value": value,
    }

    # Define the trace for the sunburst chart.
    sunburst_trace = go.Sunburst(
        labels=sunburst_data["label"],
        parents=sunburst_data["parent"],
        values=sunburst_data["value"],
        customdata=node_info,
        hovertemplate="<b>%{label}</b><br>%{customdata}<extra></extra>",
        marker=dict(colors=colors),
    )

    sunburst_figure = go.Figure(data=sunburst_trace)
    sunburst_figure.update_layout(
        font=dict(
            family="Arial, sans-serif",  # Use a professional sans-serif font like Arial.
            size=14,  # Adjust font size as needed.
            color="black",  # Set font color.
        )
    )
    sunburst_figure.update_layout(width=1500)
    sunburst_figure.update_layout(
        hoverlabel_namelength=0,
        title="Sequence typing: sunburst chart",
    )

    return sunburst_figure
