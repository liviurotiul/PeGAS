import os

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.colors import qualitative


def create_contig_plot(df, species_dict):
    species_contig_plots = {}
    n50_rows, cov_rows = [], []

    df = df[['SAMPLE', 'SPECIES']].dropna().astype(str)
    df = df.drop_duplicates(subset=['SAMPLE'], keep='first')
    df['SpeciesName'] = df['SPECIES'].map(species_dict).fillna(df['SPECIES'])

    pastel = (qualitative.Pastel1 + qualitative.Pastel2)
    set2 = qualitative.Set2
    species_order = sorted(df['SpeciesName'].unique())
    per_species_idx = {sp: 0 for sp in species_order}

    for _, row in df.iterrows():
        sample = row['SAMPLE']
        sp_name = row['SpeciesName']
        contigs_file = f"results/{sample}/shovill/contigs.fa"
        if not os.path.exists(contigs_file):
            continue

        lens, covs = [], []
        with open(contigs_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    parts = line.strip().split()
                    try:
                        L = int([x for x in parts if x.startswith('len=')][0].split('=')[1])
                        C = float([x for x in parts if x.startswith('cov=')][0].split('=')[1])
                    except Exception:
                        continue
                    lens.append(L)
                    covs.append(C)
        if not lens:
            continue

        pairs = sorted(zip(lens, covs), key=lambda t: t[0], reverse=True)
        sorted_lens = np.array([p[0] for p in pairs], dtype=int)
        sorted_covs = np.array([p[1] for p in pairs], dtype=float)

        # Assumes you have this helper defined elsewhere.
        n50_value, n50_idx = calculate_n50_and_index(sorted_lens)
        n50_cov = float(sorted_covs[n50_idx])

        n50_rows.append({'Sample': sample, 'N50': int(n50_value), 'Species': sp_name})
        cov_rows.append({'Sample': sample, 'N50_Coverage': n50_cov, 'Species': sp_name})

        cum = np.insert(np.cumsum(sorted_lens), 0, 0)
        x = np.arange(cum.size)
        x_n50 = x[:n50_idx + 2]
        y_n50 = cum[:n50_idx + 2]
        x_rest = x[n50_idx + 1:]
        y_rest = cum[n50_idx + 1:]

        fig = species_contig_plots.setdefault(
            sp_name,
            go.Figure(layout=dict(
                title=f'Cumulative Contig Lengths for Species {sp_name}',
                xaxis_title='Contig Index',
                yaxis_title='Cumulative Contig Length (bp)',
                legend_title='Samples',
                template='plotly_white',
                height=600,
                xaxis=dict(range=[0, None]),
                yaxis=dict(range=[0, None])
            ))
        )
        ci = per_species_idx[sp_name] % len(pastel)
        color = pastel[ci]
        per_species_idx[sp_name] += 1

        fig.add_trace(go.Scatter(
            x=x_n50, y=y_n50, mode='lines',
            name=f'{sample} (to N50)',
            line=dict(color='red', width=3),
            hovertemplate='Sample: ' + sample + '<br>Contig Index: %{x}<br>Cumulative Length: %{y} bp<extra></extra>'
        ))
        fig.add_trace(go.Scatter(
            x=x_rest, y=y_rest, mode='lines',
            name=sample,
            line=dict(color=color, width=3),
            hovertemplate='Sample: ' + sample + '<br>Contig Index: %{x}<br>Cumulative Length: %{y} bp<extra></extra>'
        ))

    n50_df = pd.DataFrame(n50_rows)
    coverage_df = pd.DataFrame(cov_rows)

    # Guard against empty data.
    if n50_df.empty:
        n50_box_fig = go.Figure()
        coverage_box_fig = go.Figure()
        return species_contig_plots, n50_box_fig, coverage_box_fig

    n50_df['Species'] = pd.Categorical(n50_df['Species'], categories=species_order, ordered=True)
    coverage_df['Species'] = pd.Categorical(coverage_df['Species'], categories=species_order, ordered=True)

    n50_box_fig = px.box(
        n50_df, x='Species', y='N50', color='Species',
        points='all', title='N50 Contig Length Distribution Across Species',
        template='plotly_white', color_discrete_sequence=set2,
        custom_data=['Sample']
    )
    n50_box_fig.update_traces(
        jitter=0.3, pointpos=-1.8, marker_size=6,
        hovertemplate='Sample: %{customdata[0]}<br>Species: %{x}<br>N50: %{y} bp<extra></extra>'
    )
    n50_box_fig.update_layout(yaxis_title='N50 Contig Length (bp)', xaxis_title='Species', height=400)

    coverage_box_fig = px.box(
        coverage_df, x='Species', y='N50_Coverage', color='Species',
        points='all', title='N50 Contig Coverage Across Species',
        template='plotly_white', color_discrete_sequence=set2,
        custom_data=['Sample']
    )
    coverage_box_fig.update_traces(
        jitter=0.3, pointpos=-1.8, marker_size=6,
        hovertemplate='Sample: %{customdata[0]}<br>Species: %{x}<br>Coverage: %{y:.2f}X<extra></extra>'
    )
    coverage_box_fig.update_layout(yaxis_title='N50 Contig Coverage (X)', xaxis_title='Species', height=400)

    return species_contig_plots, n50_box_fig, coverage_box_fig


def calculate_n50_and_index(sorted_lengths):
    """
    Calculates the N50 value and the index where N50 is reached from a list of contig lengths sorted in descending order.

    Parameters:
    - sorted_lengths (list): List of contig lengths sorted in descending order.

    Returns:
    - n50 (int): The N50 contig length.
    - index (int): The index at which N50 is reached.
    """
    total_length = sum(sorted_lengths)
    half_total = total_length / 2
    running_total = 0
    for idx, length in enumerate(sorted_lengths):
        running_total += length
        if running_total >= half_total:
            return length, idx
    return 0, 0
