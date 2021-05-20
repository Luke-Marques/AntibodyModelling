import pandas as pd
import numpy as np
import itertools
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import seaborn as sns
import get_mutations as mut

%matplotlib inline
sns.set_theme()
sns.set_context('paper')
sns.set_palette('crest')

colours = sns.color_palette(palette = 'crest', n_colors = 6)

patient_id_list = [
    'YF189',
    'YF191',
    'YF192',
    'YF199',
    'YF200',
    'YF203'
]


def draw_MutPosition_hist(df, axes):
    max_MutPosition = df.MutPosition.max()
    return sns.histplot(
        ax = axes,
        data = df,
        weights = np.ones(len(df)) / len(df),
        x = 'MutPosition',
        kde = True,
        bins = max_MutPosition,
        color = colours[4])


def create_MutPosition_hists(patient_id_list):
    '''Creates histograms of MutPosition occurances for each patient_id and
    across all patient_ids given the patient_id_list arg'''

    # creates dir for plots if it doesn't already exists
    plot_dir_path = os.getcwd() + os.sep + 'results' + os.sep + \
        'mutation_data' + os.sep + 'MutPosition_plots'
    if not os.path.exists(plot_dir_path):
        os.makedirs(plot_dir_path)

    # creates list to store individual patient_id dfs in for aggregate histogram
    df_light_list = []
    df_heavy_list = []

    # loops through patient_ids, appends dataframes to df_lists and creates
    # histograms of MutPosition distributions for both heavy and light chains
    for patient_id in patient_id_list:
        for chain_type in ['Light', 'Heavy']:
            if chain_type == 'Light':
                df = mut.create_mut_dataframe_light(patient_id)
                df_light_list.append(df)
            elif chain_type == 'Heavy':
                df = mut.create_mut_dataframe_heavy(patient_id)
                df_heavy_list.append(df)
            fig, ax = plt.subplots(figsize=(7,5))
            draw_MutPosition_hist(df, ax)
            ax.set_xlabel('Mutation Position')
            ax.set_ylabel('Counts')
            ax.set_xlim(left=0, right=120)
            ax.set_ylim(bottom=0, top=0.06)
            ax.yaxis.set_major_formatter(PercentFormatter(1))
            ax.set_title(patient_id + \
                ' - Distribution of Mutation Posistions within' + chain_type \
                + ' Chain')
            plot_path = plot_dir_path + os.sep + patient_id + chain_type + \
                'MutPosition_hist.pdf'
            fig.savefig(plot_path)

    # combines all patient's light chain mutation dataframes to single dataframe
    df_light_all = pd.concat(df_light_list)

    # combines all patient's heavy chain mutation dataframes to single dataframe
    df_heavy_all = pd.concat(df_heavy_list)

    # draws histogram of mutation positions of heavy and light chains of
    # all patients
    for chain_type in ['Light', 'Heavy']:
        if chain_type == 'Light':
            df = df_light_all
        elif chain_type == 'Heavy':
            df = df_heavy_all
        fig, ax = plt.subplots(figsize=(7,5))
        draw_MutPosition_hist(df, ax)
        ax.set_xlabel('Mutation Position')
        ax.set_ylabel('Counts')
        ax.set_xlim(left=0, right=120)
        ax.set_ylim(bottom=0, top=0.06)
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.set_title('All Patients - Distribution of Mutation Posistions ' + \
            'within ' + chain_type + ' Chain')
        plot_path = plot_dir_path + os.sep + 'AllPatients' + chain_type + \
            'MutPosition_hist.pdf'
        fig.savefig(plot_path)


def draw_Region_barplot(df, axes):

    df = df.groupby('Region').count()[['SeqID']] \
        .rename({'SeqID': 'Counts'}, axis=1).reset_index()

    # stores lengths of each region as variables
    FR1_len = 26-0
    FR2_len = 55-38
    FR3_len = 104-65
    CDR1_len = 38-26
    CDR2_len = 65-55
    CDR3_len = 117-104

    df.loc[df.Region == 'FR1', 'Prop_Counts'] = df.Counts / FR1_len
    df.loc[df.Region == 'FR2', 'Prop_Counts'] = df.Counts / FR2_len
    df.loc[df.Region == 'FR3', 'Prop_Counts'] = df.Counts / FR3_len
    df.loc[df.Region == 'CDR1', 'Prop_Counts'] = df.Counts / CDR1_len
    df.loc[df.Region == 'CDR2', 'Prop_Counts'] = df.Counts / CDR2_len
    df.loc[df.Region == 'CDR3', 'Prop_Counts'] = df.Counts / CDR3_len

    Prop_Counts_total = df.Prop_Counts.sum()

    df.loc[df.Region == 'FR1', 'Prop_pct'] = \
        df.Prop_Counts / Prop_Counts_total
    df.loc[df.Region == 'FR2', 'Prop_pct'] = \
        df.Prop_Counts / Prop_Counts_total
    df.loc[df.Region == 'FR3', 'Prop_pct'] = \
        df.Prop_Counts / Prop_Counts_total
    df.loc[df.Region == 'CDR1', 'Prop_pct'] = \
        df.Prop_Counts / Prop_Counts_total
    df.loc[df.Region == 'CDR2', 'Prop_pct'] = \
        df.Prop_Counts / Prop_Counts_total
    df.loc[df.Region == 'CDR3', 'Prop_pct'] = \
        df.Prop_Counts / Prop_Counts_total

    df = df.reindex([3, 0, 4, 1, 5, 2]).reset_index()

    return sns.barplot(
        ax = axes,
        x = 'Region',
        y = 'Prop_pct',
        data = df,
        palette = 'crest')


def create_Region_barplots(patient_id_list):
    '''Creates bar plots of the mutation Region occurances for each patient_id
    and across all patient_ids, proportional to the number of amino acids in the
    region'''

    # creates dir for plots if it doesn't already exists
    plot_dir_path = os.getcwd() + os.sep + 'results' + os.sep + \
        'mutation_data' + os.sep + 'Region_plots'
    if not os.path.exists(plot_dir_path):
        os.makedirs(plot_dir_path)

    # creates list to store individual patient_id dfs in for aggregate stats
    df_light_list = []
    df_heavy_list = []

    # loops through patient_ids, appends dataframes to df_list and creates
    # barplots of Region mutation distributions
    for patient_id in patient_id_list:
        for chain_type in ['Light', 'Heavy']:
            if chain_type == 'Light':
                df = mut.create_mut_dataframe_light(patient_id)
                df = df[df.Region != 'CDR3'] # can remove later on
                df_light_list.append(df)
            elif chain_type == 'Heavy':
                df = mut.create_mut_dataframe_heavy(patient_id)
                df = df[df.Region != 'CDR3'] # can remove later on
                df_heavy_list.append(df)

            fig, ax = plt.subplots(figsize=(7,5))
            draw_Region_barplot(df, ax)
            ax.set_xlabel('Mutation Region')
            ax.set_ylabel('Proportional Percentage of Mutations')
            ax.set_ylim(bottom=0, top=0.4)
            ax.yaxis.set_major_formatter(PercentFormatter(1))
            ax.set_title(patient_id + \
                ' - Proportional Distribution of Mutations within Regions')
            plot_path = plot_dir_path + os.sep + patient_id + \
                        '_muts_by_Region_within_ ' + chain_type + '_Chain.pdf'
            fig.savefig(plot_path)

    # combines all patient's light chain mutation dataframes to single dataframe
    df_light_all = pd.concat(df_light_list)

    # combines all patient's heavy chain mutation dataframes to single dataframe
    df_heavy_all = pd.concat(df_heavy_list)

    for chain_type in ['Light', 'Heavy']:

        if chain_type == 'Light':
            df = df_light_all
        elif chain_type == 'Heavy':
            df = df_heavy_all

        fig, ax = plt.subplots(figsize=(7,5))
        draw_Region_barplot(df, ax)
        ax.set_xlabel('Mutation Region')
        ax.set_ylabel('Proportional Percentage of Mutations')
        ax.set_ylim(bottom=0, top=0.4)
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.set_title('All Patients ' + \
            ' - Proportional Distribution of Mutations within Regions')
        plot_path = plot_dir_path + os.sep + \
            'AllPatients_muts_by_Region_within_' + chain_type + '_Chain.pdf'
        fig.savefig(plot_path)


def draw_aa_heatmap(df, axes):
    '''Takes the mutation dataframe from get_mutations.py create_mut_dataframe
    (light or heavy) and axes for plot to be assigned to as args and returns a
    heatmap of amino-acid mutations'''

    # sets up lists for ordering of columns and index of heatmap df
    charged_aminoacids = ['R', 'H', 'K', 'D', 'E']
    poscharged_aminoacids = ['R', 'H', 'K']
    negcharged_aminoacids = ['D', 'E']
    uncharged_aminoacids = ['Q', 'N', 'S', 'T', 'Y', 'C']
    nonpolar_aminoacids = ['A', 'I', 'L', 'M', 'F', 'V', 'P', 'G', 'W']
    ordered_aminoacids = nonpolar_aminoacids + \
        uncharged_aminoacids + \
        negcharged_aminoacids + \
        poscharged_aminoacids
    reversed_aminoacids = list(reversed(ordered_aminoacids))

    # creates pivot table from mutation df and converts count to percentages
    df = df[df.sequence_aa != '.']
    df = pd.crosstab(df.sequence_aa, df.germline_aa)
    df = df.div(df.to_numpy().sum())
    df = df * 100
    df = df.reindex(reversed_aminoacids)
    df = df.reindex(columns = ordered_aminoacids)

    return sns.heatmap(data=df, cmap='crest')


def draw_grouped_aa_heatmap(df, axes):
    '''Takes the mutation dataframe from get_mutations.py create_mut_dataframe
    (light or heavy) and axes for plot to be assigned to as args and returns a
    heatmap of amino-acid type mutations'''

    # sets up lists for ordering of columns and index of heatmap df
    ordering = ['Non-Polar',
                'Uncharged',
                'Positively Charged',
                'Negatively Charged']
    reversed_ordering = list(reversed(ordering))

    # creates pivot table from mutation df and converts count to percentages
    df = df[df.sequence_aa != '.']
    df = pd.crosstab(df.germline_aa_type, df.sequence_aa_type)
    df = df.div(df.to_numpy().sum())
    df = df * 100
    df = df.reindex(reversed_ordering)
    df = df.reindex(columns = ordering)

    return sns.heatmap(data=df, annot=True, cmap='crest')


def create_mutation_heatmaps_per_patient_per_chain(patient_id_list):
    '''Creates and saves mutation heatmaps (grouped and ungrouped) of each chain
    type (heavy and light) for each patient in the patient_id_list arg'''

    # loops through patient_ids, appends dataframes to df_list and creates
    # 'histograms' of Region distributions for light and heavy chains
    for patient_id in patient_id_list:

        for chain_type in ['Light', 'Heavy']:

            # sets up list to store region dataframes and overall dataframes of
            # each patient and chain type
            df_list = []

            if chain_type == 'Light':

                # generates mutation dataframe
                df = mut.create_mut_dataframe_light(patient_id)
                df = df[df.Region != 'CDR3'] # can remove to include CDR3
                df_list.append(df)

            elif chain_type == 'Heavy':

                # generates mutation dataframe
                df = mut.create_mut_dataframe_heavy(patient_id)
                df = df[df.Region != 'CDR3'] # can remove to include CDR3
                df_list.append(df)

            # splits dataframe by 'Region' values
            df_CDR1 = df[df.Region == 'CDR1'].copy()
            df_list.append(df_CDR1)
            df_CDR2 = df[df.Region == 'CDR2'].copy()
            df_list.append(df_CDR2)
            df_FR1 = df[df.Region == 'FR1'].copy()
            df_list.append(df_FR1)
            df_FR2 = df[df.Region == 'FR2'].copy()
            df_list.append(df_FR2)
            df_FR3 = df[df.Region == 'FR3'].copy()
            df_list.append(df_FR3)

            # create dir for patient_id plots if it doesn't already exists
            patient_plot_dir_path = os.getcwd() + os.sep + 'results' + \
                os.sep + 'mutation_data' + os.sep + 'Mut_Heatmaps' + \
                os.sep + patient_id
            chain_plot_dir_path = patient_plot_dir_path + os.sep + \
                chain_type
            if not os.path.exists(chain_plot_dir_path):
                os.makedirs(chain_plot_dir_path)

            # stores names of regions in list
            df_names = ['All_Regions', 'CDR1', 'CDR2', 'FR1', 'FR2', 'FR3']

            counter = 0
            for df in df_list:

                # gets name of region from df_names list
                name = df_names[counter]

                # draws ungrouped amino acid mutation heatmap
                fig_ungrouped, ax = plt.subplots(figsize=(7,5))
                draw_aa_heatmap(df, ax)
                ax.set_xlabel('Germline Amino Acid')
                ax.set_ylabel('Mutated Amino Acid')
                title = patient_id + ' -  ' + chain_type + ' Chain' + ' - ' + \
                    name + ' - Individual Amino-Acid Mutation Heatmap'
                ax.set_title(title)
                plot_path = chain_plot_dir_path + os.sep + patient_id + \
                    '_' + chain_type + '_' + name + '_Mut_Heatmap.pdf'
                fig_ungrouped.savefig(plot_path, bbox_inches = "tight")

                # draws grouped amino acid mutation heatmap
                fig_grouped, ax = plt.subplots(figsize=(7,5))
                draw_grouped_aa_heatmap(df, ax)
                ax.set_xlabel('Germline Amino Acid')
                ax.set_ylabel('Mutated Amino Acid')
                title = patient_id + ' - ' + chain_type + ' Chain' + ' - ' + \
                    name + ' - Grouped Amino-Acid Mutation Heatmap'
                ax.set_title(title)
                plot_path = chain_plot_dir_path + os.sep + patient_id + \
                    '_' + chain_type + '_' + name + '_Grouped_Mut_Heatmap.pdf'
                fig_grouped.savefig(plot_path, bbox_inches = "tight")

                counter += 1

            # resets df_list to empty
            df_list = []


def create_mutation_heatmaps_all_patients_per_chain(patient_id_list):
    '''Creates and saves mutation heatmaps (grouped and ungrouped) of each chain
    type (heavy and light) from mutation data of all patients the
    patient_id_list arg, concatenated'''


    # creates lists to store individual patient_id dfs in for aggregate stats
    # light chain dfs:
    df_CDR1_light_list = []
    df_CDR2_light_list = []
    df_FR1_light_list = []
    df_FR2_light_list = []
    df_FR3_light_list = []
    df_all_light_list = []
    # heavy chain dfs:
    df_CDR1_heavy_list = []
    df_CDR2_heavy_list = []
    df_FR1_heavy_list = []
    df_FR2_heavy_list = []
    df_FR3_heavy_list = []
    # all regions dfs:
    df_all_heavy_list = []
    df_all_light_list = []

    # loops through patient_ids, appends dataframes to df_list and creates
    # 'histograms' of Region distributions for light and heavy chains
    for patient_id in patient_id_list:
        for chain_type in ['Light', 'Heavy']:
            if chain_type == 'Light':

                # generates mutation dataframe
                df = mut.create_mut_dataframe_light(patient_id)
                df = df[df.Region != 'CDR3'] # can remove to include CDR3

                # appends dataframe of all regions to list
                df_all_light_list.append(df)

                # splits dataframe by 'Region' values and appends each
                # regions dataframe to list
                df_CDR1 = df[df.Region == 'CDR1'].copy()
                df_CDR1_light_list.append(df_CDR1)
                df_CDR2 = df[df.Region == 'CDR2'].copy()
                df_CDR2_light_list.append(df_CDR2)
                df_FR1 = df[df.Region == 'FR1'].copy()
                df_FR1_light_list.append(df_FR1)
                df_FR2 = df[df.Region == 'FR2'].copy()
                df_FR2_light_list.append(df_FR2)
                df_FR3 = df[df.Region == 'FR3'].copy()
                df_FR3_light_list.append(df_FR3)

            elif chain_type == 'Heavy':

                # generates mutation dataframe
                df = mut.create_mut_dataframe_heavy(patient_id)
                df = df[df.Region != 'CDR3'] # can remove to include CDR3

                # appends dataframe of all regions to list
                df_all_heavy_list.append(df)

                # splits dataframe by 'Region' values and appends each
                # regions dataframe to list
                df_CDR1 = df[df.Region == 'CDR1'].copy()
                df_CDR1_heavy_list.append(df_CDR1)
                df_CDR2 = df[df.Region == 'CDR2'].copy()
                df_CDR2_heavy_list.append(df_CDR2)
                df_FR1 = df[df.Region == 'FR1'].copy()
                df_FR1_heavy_list.append(df_FR1)
                df_FR2 = df[df.Region == 'FR2'].copy()
                df_FR2_heavy_list.append(df_FR2)
                df_FR3 = df[df.Region == 'FR3'].copy()
                df_FR3_heavy_list.append(df_FR3)

    # combines all dataframes in lists
    df_CDR1_light = pd.concat(df_CDR1_light_list)
    df_CDR1_heavy = pd.concat(df_CDR1_heavy_list)
    df_CDR2_light = pd.concat(df_CDR2_light_list)
    df_CDR2_heavy = pd.concat(df_CDR2_heavy_list)
    df_FR1_light = pd.concat(df_FR1_light_list)
    df_FR1_heavy = pd.concat(df_FR1_heavy_list)
    df_FR2_light = pd.concat(df_FR2_light_list)
    df_FR2_heavy = pd.concat(df_FR2_heavy_list)
    df_FR3_light = pd.concat(df_FR3_light_list)
    df_FR3_heavy = pd.concat(df_FR3_heavy_list)
    df_all_regions_heavy = pd.concat(df_all_heavy_list)
    df_all_regions_light = pd.concat(df_all_light_list)

    # stores each region df of same chain type in list to loop through
    df_heavy_list = [
        df_CDR1_heavy,
        df_CDR2_heavy,
        df_FR1_heavy,
        df_FR2_heavy,
        df_FR3_heavy,
        df_all_regions_heavy]
    df_light_list = [
        df_CDR1_light,
        df_CDR2_light,
        df_FR1_light,
        df_FR2_light,
        df_FR3_light,
        df_all_regions_light]

    # stores names of regions in list
    df_names = ['CDR1', 'CDR2', 'FR1', 'FR2', 'FR3', 'All_Regions']

    for chain_type in ['Light', 'Heavy']:
        if chain_type == 'Light':
            df_list = df_light_list
        elif chain_type == 'Heavy':
            df_list = df_heavy_list

        # create dir for plots if it doesn't already exists
        all_patients_plot_dir_path = os.getcwd() + os.sep + 'results' + \
            os.sep + 'mutation_data' + os.sep + 'Mut_Heatmaps' + \
            os.sep + 'All_Patients'
        chain_plot_dir_path = all_patients_plot_dir_path + os.sep + \
            chain_type
        if not os.path.exists(chain_plot_dir_path):
            os.makedirs(chain_plot_dir_path)

        counter = 0
        for df in df_list:

            # gets name of region from df_names list
            name = df_names[counter]

            # draws ungrouped amino acid mutation heatmap
            fig_ungrouped, ax = plt.subplots(figsize=(7,5))
            draw_aa_heatmap(df, ax)
            ax.set_xlabel('Germline Amino Acid')
            ax.set_ylabel('Mutated Amino Acid')
            title = 'All Patients' + ' - ' + chain_type + ' Chain' + ' - ' + \
                name + ' - Individual Amino-Acid Mutation Heatmap'
            ax.set_title(title)
            plot_path = chain_plot_dir_path + os.sep + 'All-Patients_' + \
                chain_type + '_' + name + '_Mut_Heatmap.pdf'
            fig_ungrouped.savefig(plot_path, bbox_inches = "tight")

            # draws grouped amino acid mutation heatmap
            fig_grouped, ax = plt.subplots(figsize=(7,5))
            draw_grouped_aa_heatmap(df, ax)
            ax.set_xlabel('Germline Amino Acid')
            ax.set_ylabel('Mutated Amino Acid')
            title = 'All Patients' + ' - ' + chain_type + ' Chain' + ' - ' + \
                name + ' - Grouped Amino-Acid Mutation Heatmap'
            ax.set_title(title)
            plot_path = chain_plot_dir_path + os.sep + 'All-Patients_' + \
                chain_type + '_' + name + '_Grouped_Mut_Heatmap.pdf'
            fig_grouped.savefig(plot_path, bbox_inches = "tight")

            counter += 1
