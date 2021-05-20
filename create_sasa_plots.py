import pandas as pd
import numpy as np
import itertools
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import seaborn as sns
import  clean_repertoire_data as crd
import read_pops as rp

pd.options.display.html.table_schema = True
pd.options.display.max_rows = None
pd.options.display.max_columns = None

patient_id_list = ['YF189',
                   'YF191',
                   'YF192',
                   'YF199',
                   'YF200',
                   'YF203']


def get_SASA_trend_dataframe(patient_id):

    # read per_residue_SASA.csv file and convert to dataframe
    SASA_dir = os.getcwd() + os.sep + 'results' + os.sep + 'SASA_results'
    SASA_filepath = SASA_dir + os.sep + patient_id + '_per_residue_SASA.csv'
    df_SASA = pd.read_csv(SASA_filepath)

    df_SASA = df_SASA[['jobname', 'Chain', 'Q(SASA)']] \
                         .groupby(['jobname', 'Chain'], as_index=False) \
                         .agg({'Q(SASA)': 'sum'}) \
                         .reset_index(drop=True)

    df_SASA = df_SASA[['jobname', 'Q(SASA)']]
    df_SASA['g'] = (df_SASA.groupby('jobname').cumcount() + 1).astype(str)
    df_SASA = df_SASA.set_index(['jobname','g']).unstack()
    df_SASA.sort_index(axis=1, level=1, inplace=True)
    df_SASA.columns = [''.join(col) for col in df_SASA.columns]
    df_SASA.reset_index(inplace=True)
    df_SASA = df_SASA.rename(columns = {'Q(SASA)1': 'H_SASA',
                                        'Q(SASA)2': 'L_SASA'})

    # gets sequence dataframe to get v.identity for each jobname chain
    df_sequence = crd.prepare_sequence_dataframe(patient_id) \
        [['sequence_id', 'v_identity']]

    # gets paired cell dataframe to get CloneIDs for each jobname
    df_paired = crd.reshape_paired_cell_dataframe(patient_id)

    # adds SASA columns to df_paired
    df = df_paired.merge(df_SASA,
                         left_on = 'jobname',
                         right_on = 'jobname')

    # adds H_v_identity and L_v_identity columns to df
    df = df.merge(df_sequence,
                  left_on = 'H_SeqID',
                  right_on = 'sequence_id')
    df = df.rename(columns = {'v_identity': 'H_v_identity'})
    df = df.merge(df_sequence,
                  left_on = 'L_SeqID',
                  right_on = 'sequence_id')
    df = df.rename(columns = {'v_identity': 'L_v_identity'})
    df = df.drop(columns = ['sequence_id_x', 'sequence_id_y'])

    # sorts by H_CloneID and H_v_identity to see trends
    df = df.sort_values(['L_CloneID', 'L_v_identity'], ascending=[True, False])
    df = df.reset_index(drop=True)

    return df # redundant


def get_SASA_identity_correlation(patient_id):

    # gets trend DataFrame
    df = get_SASA_trend_dataframe(patient_id)[['jobname',
                                                'H_CloneID',
                                                'H_SASA',
                                                'H_v_identity',
                                                'L_CloneID',
                                                'L_SASA',
                                                'L_v_identity']]

    # creates lists of unique CloneID values
    H_CloneID_list = list(df.H_CloneID.unique())
    L_CloneID_list = list(df.L_CloneID.unique())

    # creates new empty columns for % SASA increase which will show the increase
    # in SASA values for a given antibody compared to the antibody of the same
    # clone with the highest v_identity (closest to germline)
    df['H_SASA_pc_increase'] = ''
    df['L_SASA_pc_increase'] = ''

    # creates empty lists to store correlations for each CloneID
    H_corr_list = []
    L_corr_list = []

    # sorts by H_CloneID and H_v_identity to see trends
    df = df.sort_values(['H_CloneID', 'H_v_identity'], ascending=[True, False])
    df = df.reset_index(drop=True)

    # loops through unique CloneIDs in list and slices dataframe
    for CloneID in H_CloneID_list:
        df_slice = df[df.H_CloneID == CloneID]
        df_slice = df_slice[['H_SASA', 'H_v_identity']]
        df_slice = df_slice.groupby('H_v_identity', as_index=False) \
                           .agg({'H_SASA': 'mean'}) \
                           .sort_values('H_v_identity', ascending=False) \
                           .dropna() \
                           .reset_index(drop=True)
        if len(df_slice) > 1:
            correlation = df_slice.H_v_identity.corr(df_slice.H_SASA)
            data_points = len(df_slice) - 1
            for num in range(data_points):
                H_corr_list.append(correlation)

    # sorts by H_CloneID and H_v_identity to see trends
    df = df.sort_values(['L_CloneID', 'L_v_identity'], ascending=[True, False])
    df = df.reset_index(drop=True)

    for CloneID in L_CloneID_list:
        df_slice = df[df.L_CloneID == CloneID]
        df_slice = df_slice[['L_SASA', 'L_v_identity']]
        df_slice = df_slice.groupby('L_v_identity', as_index=False) \
                           .agg({'L_SASA': 'mean'}) \
                           .sort_values('L_v_identity', ascending=False) \
                           .dropna() \
                           .reset_index(drop=True)
        if len(df_slice) > 1:
            correlation = df_slice.L_v_identity.corr(df_slice.L_SASA)
            data_points = len(df_slice) - 1
            for num in range(data_points):
                L_corr_list.append(correlation)

    # calculates weighted mean correlation from all clones
    H_corr = np.mean(H_corr_list)
    L_corr = np.mean(L_corr_list)

    return H_corr, L_corr # redundant


def show_SASA_identity_correlation(patient_id_list):

    H_corr_list_all = []
    L_corr_list_all = []
    H_corr_list_young = []
    L_corr_list_young = []
    H_corr_list_old = []
    L_corr_list_old = []

    for patient_id in patient_id_list:

        old_list = ['YF189',
                    'YF191',
                    'YF192']

        H_corr = get_SASA_identity_correlation(patient_id)[0]
        if H_corr != 1.0 or -1.0:
            H_corr_list_all.append(H_corr)
        L_corr = get_SASA_identity_correlation(patient_id)[1]
        if L_corr != 1.0 or -1.0:
            L_corr_list_all.append(L_corr)

        if patient_id in old_list:
            if H_corr != 1.0 or -1.0:
                H_corr_list_old.append(H_corr)
            if L_corr != 1.0 or -1.0:
                L_corr_list_old.append(L_corr)
        else:
            if H_corr != 1.0 or -1.0:
                H_corr_list_young.append(H_corr)
            if L_corr != 1.0 or -1.0:
                L_corr_list_young.append(L_corr)

        print(patient_id + ':')
        print('    H_correlation: ' + str(H_corr))
        print('    L_correlation: ' + str(L_corr))

    avg_H_corr = np.mean(H_corr_list_all)
    avg_L_corr = np.mean(H_corr_list_all)
    avg_H_corr_young = np.mean(H_corr_list_young)
    avg_L_corr_young = np.mean(L_corr_list_young)
    avg_H_corr_old = np.mean(H_corr_list_old)
    avg_L_corr_old = np.mean(L_corr_list_old)

    print('all_patients:')
    print('    H_correlation: ' + str(avg_H_corr))
    print('    L_correlation: ' + str(avg_L_corr))

    print('young_patients:')
    print('    H_correlation: ' + str(avg_H_corr_young))
    print('    L_correlation: ' + str(avg_L_corr_young))

    print('old_patients:')
    print('    H_correlation: ' + str(avg_H_corr_old))
    print('    L_correlation: ' + str(avg_L_corr_old)) # redundant


def sort_alphanumeric_list(list):

    # convert list to string, with each object joined with comma
    list_string = ', '.join(list)

    # replace letters in string with decimals
    list_string = list_string.replace('A', '.1') \
                             .replace('B', '.2') \
                             .replace('C', '.3') \
                             .replace('D', '.4') \
                             .replace('E', '.5') \
                             .replace('F', '.6') \
                             .replace('G', '.7') \
                             .replace('H', '.8') \
                             .replace('I', '.9')

    # convert string back to list
    listOfStrings = list_string.split(', ')

    # convert objects in list to floats
    listOfIntergers = []
    for value in listOfStrings:
        value_float = float(value)
        value_int = int(value_float *  10)
        listOfIntergers.append(value_int)

    # sort interger list numerically
    listOfIntergers_sorted = sorted(listOfIntergers)

    # convert objects in sorted list back to strings
    listOfStrings_sorted = []
    for value in listOfIntergers_sorted:
        value_float = value / 10
        value_string = str(value_float)
        listOfStrings_sorted.append(value_string)

    # convert list to string, with each object joined with comma
    list_string = ', '.join(listOfStrings_sorted)

    # replaces decimals with original letters
    list_string = list_string.replace('.1', 'A') \
                             .replace('.2', 'B') \
                             .replace('.3', 'C') \
                             .replace('.4', 'D') \
                             .replace('.5', 'E') \
                             .replace('.6', 'F') \
                             .replace('.7', 'G') \
                             .replace('.8', 'H') \
                             .replace('.9', 'I') \
                             .replace('.0', '')

    # convert string back to list
    list_sorted = list_string.split(', ')

    return list_sorted


def draw_perResSASA_heatmap(patient_id_list):

    # df_H = rp.create_SASA_meta_dataframe_all_patients(patient_id_list)[0]
    # df_L = rp.create_SASA_meta_dataframe_all_patients(patient_id_list)[1]

    SASA_dir = os.getcwd() + os.sep + 'results' + os.sep + 'SASA_results'
    filepath_H = SASA_dir + os.sep + 'all_patients_SASA_meta_H.csv'
    filepath_L = SASA_dir + os.sep + 'all_patients_SASA_meta_L.csv'
    df_H = pd.read_csv(filepath_H)
    df_L = pd.read_csv(filepath_L)

    for chain_type in 'Heavy', 'Light':

        if chain_type == 'Heavy':
            df = df_H
        elif chain_type == 'Light':
            df = df_L

        reorder_cols = 'n'
        if df.Position_Insertion.dtypes != 'int64':
            insertion_list = ['111A', '111B', '111C', '111D', '111E', '111F',
                              '112F', '112E', '112D', '112C', '112B', '112A']
            cols_sorted = list(range(1, 112, 1)) + insertion_list + list(range(112, 129, 1))
            cols_sorted = [x if (type(x) is str) else str(x) for x in cols_sorted]
            reorder_cols = 'y'

        df = df.pivot(index = 'jobname',
                      columns = 'Position_Insertion',
                      values = 'Q(SASA)')

        if reorder_cols == 'y':
            df = df.reindex(columns=cols_sorted)

        # draws seaborn heatmap from dataframe
        fig, ax = plt.subplots(figsize=(12,4))
        ax = sns.heatmap(data=df,
                         cmap='mako_r',
                         yticklabels=False,
                         vmin = 0,
                         vmax = 1.8,
                         cbar_kws = {'label': 'Residue SASA'})
        ax.set_xticks(['1','27','39','56','66','105','118','128'], ' ')
        ax.set_xlabel('Residue Position', size=14)
        ax.set_ylabel('Ab V-Region', size=14)
        title = 'Q(SASA) Across All Residues in All Ab V-Region Models - ' + \
                chain_type + ' Chain'
        ax.set_title(title, size=18)

        # creates dir for plot
        plot_dir = os.getcwd() + os.sep + 'results' + os.sep + 'SASA_results' \
                    + os.sep + 'plots'
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        plot_path = plot_dir + os.sep + 'SASA_heatmap_' + chain_type + '.png'

        # saves plot to dir
        fig.savefig(plot_path, bbox_inches = "tight")


def get_avgNumExposedPerRegion(patient_id, chain):

    # get SASA dataframe
    filepath = os.getcwd() + os.sep + 'results' + os.sep + 'SASA_results' + \
               os.sep + 'all_patients_SASA_meta_' + chain + '.csv'
    df = pd.read_csv(filepath)

    # slice dataframe on patient_id
    df = df[df.patient_id == patient_id]

    # adds column 'Position' which does not use insertion codes
    print('trying')
    Positions = df.Position_Insertion.values.tolist()
    Positions = [str(x) for x in Positions]
    Positions = [x[:3] for x in Positions]
    df['Position'] = Positions
    df.Position = df.Position.astype(int)
    print('done')

    # adds 'Region' columns using conditions based on IMGT numbering scheme
    conditions = [
        (df.Position > 0) & (df.Position <= 26),
        (df.Position > 26) & (df.Position <= 38),
        (df.Position > 38) & (df.Position <= 55),
        (df.Position > 55) & (df.Position <= 65),
        (df.Position > 65) & (df.Position <= 104),
        (df.Position > 104) & (df.Position <= 117),
        (df.Position > 117),]
    Regions = [
        'FR1',
        'CDR1',
        'FR2',
        'CDR2',
        'FR3',
        'CDR3',
        'FR4']
    df['Region'] = np.select(conditions, Regions)

    # loops through CloneIDs and finds number of exposed/buried residues in each
    # region and stores data in dictionary
    d = {}
    for Region in Regions:
        d[Region] = {'Exposed': [], 'Buried': [], 'Unknown': []}
    if chain == 'H':
        df['CloneID'] = df.H_CloneID
    elif chain == 'L':
        df['CloneID'] = df.L_CloneID
    CloneID_list = df.CloneID.unique()
    for CloneID in CloneID_list:
        df_CloneID = df[df.CloneID == CloneID]
        for Region in Regions:
            df_Region = df_CloneID[df_CloneID.Region == Region]
            numNaN = df_Region.Exposed.isna().sum()
            numTrue = df_Region.Exposed.sum()
            numFalse = len(df_Region) - numNaN - numTrue
            d[Region]['Unknown'].append(numNaN)
            d[Region]['Exposed'].append(numTrue)
            d[Region]['Buried'].append(numFalse)

    for Region in Regions:
        avgExposed = np.mean(d[Region]['Exposed'])
        avgBuried = np.mean(d[Region]['Buried'])
        avgUnknown = np.mean(d[Region]['Unknown'])
        d[Region]['Exposed'] = avgExposed
        d[Region]['Buried'] = avgBuried
        d[Region]['Unknown'] = avgUnknown

    print(d)

    return d


draw_perResSASA_heatmap(patient_id_list)
