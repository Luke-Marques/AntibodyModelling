import os
import re
import glob
import numpy as np
import pandas as pd
import clean_repertoire_data as crd

pd.options.display.html.table_schema = True
pd.options.display.max_rows = None

patient_id_list = ['YF189',
                   'YF191',
                   'YF192',
                   'YF199',
                   'YF200',
                   'YF203']


def extractSASA(pops_file_path):
    '''Parses POPS output file to extract SASA values for each residue and
    returns dictionary of residue, chain type, position and SASA'''

    pat = re.compile('===[\s{1}]RESIDUE[\s{1}]SASAs[\s{1}]===(.*?)' + \
                     '\n===[\s{1}]MOLECULE[\s{1}]SASAs', re.M|re.S)

    # gets jobname from pops_file_path
    filename = pops_file_path.split(os.sep)[-1]
    jobname = filename.split('.')[0]

    # opens and reads POPS out file
    if os.path.isfile(pops_file_path):
        with open(pops_file_path) as f:
                text = f.read()
                for match in pat.finditer(text):
                        lines = match.group(0).split('\n')

        # get list of ResidNRs (positions) from POPS output file
        positions = list(set([l.split()[2] for l in lines[1:] \
                            if len(l.split()) > 3 \
                            and l.split()[2] != "ResidNr"]))


        # gets list of relevant lines from POPS output file, with each column of
        # each line stored as an item in a list
        subset = [line.split() for line in lines]
        subset = [line for line in subset \
                  if len(line) == 10]

        # removes first line from subset which contains column names and stores
        # removed line as a list, col_names
        col_names = subset.pop(0)

        # converts list of lines to df
        df = pd.DataFrame(subset, columns=col_names)

        # creates column position which combines columns ResidNr and
        # iCode when iCode != '-'
        df['ResidNr'] = df['ResidNr'].astype(str)
        df['Position'] = df['ResidNr']
        df['Position_Insertion'] = df['ResidNr']
        df.loc[df['iCode'] != '-', 'Position_Insertion'] = \
            df['ResidNr'] + df['iCode']

        # selects relevant columns from def
        df['jobname'] = jobname
        df['Q(SASA)'] = df['Q(SASA)'].fillna(0).astype(float)
        df.loc[df['Q(SASA)'] > 0.15, 'Exposed'] = True
        df.loc[df['Q(SASA)'] <= 0.15, 'Exposed'] = False
        df = df[['jobname',
                 'Position',
                 'Position_Insertion',
                 'Chain',
                 'Q(SASA)',
                 'Exposed']]

    return df


def create_SASA_dataframe(patient_id, save=False):
    '''Takes a string representing the patient_id as arg and generates
    dataframes for all POPS out file in the pops directory, combines all
    dataframes and (optional) saves the dataframe as a .csv file in the
    pops directory'''

    # finds directory of POPS out files
    cwd = os.getcwd()
    pops_dir = cwd + os.sep + 'results' + os.sep + 'ABodyBuilder_results' + \
               os.sep + patient_id + os.sep + 'pops'
    pops_file_paths = glob.glob(pops_dir + os.sep + '*.pops')

    # creates empty list to store POPS out file dataframes
    pops_df_list = []

    # loops through each POPS out file, creates dataframe and stores in list
    for pops_file_path in pops_file_paths:
        df_SASA = extractSASA(pops_file_path)
        pops_df_list.append(df_SASA)

    # combines all dfs into one df of all jobs for the patient
    df_SASA_all = pd.concat(pops_df_list)

    # if save arg is set to True then the following code will save the DataFrame
    # as a .csv file in SASA_results
    if save == True:
        SASA_dir = cwd + os.sep + 'results' + os.sep + 'SASA_results'
        filepath = SASA_dir + os.sep + patient_id + '_per_residue_SASA.csv'
        df_SASA_all.to_csv(filepath, index=False)

    return df_SASA_all


def create_SASA_meta_dataframe(patient_id, save=False):

    # gets SASA dataframe
    df_SASA = create_SASA_dataframe(patient_id)

    # slices SASA dataframe on Chain value
    df_SASA_H = df_SASA[df_SASA.Chain == 'H']
    df_SASA_L = df_SASA[df_SASA.Chain == 'L']

    # gets metadata of each jobname from clean_repertoire_data.py dataframe
    df_meta = crd.reshape_paired_cell_dataframe(patient_id)

    # gets sequence dataframe to get v.identity for each jobname chain
    df_identity = crd.prepare_sequence_dataframe(patient_id) \
                  [['sequence_id', 'v_identity']]

    # adds H_v_identity and L_v_identity columns to df_meta
    df_meta = df_meta.merge(df_identity,
                  left_on = 'H_SeqID',
                  right_on = 'sequence_id')
    df_meta = df_meta.rename(columns = {'v_identity': 'H_v_identity'})
    df_meta = df_meta.merge(df_identity,
                  left_on = 'L_SeqID',
                  right_on = 'sequence_id')
    df_meta = df_meta.rename(columns = {'v_identity': 'L_v_identity'})
    df_meta = df_meta.drop(columns = ['sequence_id_x', 'sequence_id_y'])

    # merges dataframes on jobname column
    df_SASA_H = df_SASA_H.merge(df_meta, left_on='jobname', right_on='jobname')
    df_SASA_L = df_SASA_L.merge(df_meta, left_on='jobname', right_on='jobname')

    # drops unnecessary columns
    df_SASA_H = df_SASA_H[['jobname',
                           'Chain',
                           'H_CloneID',
                           'H_SeqID',
                           'H_v_identity',
                           'Position_Insertion',
                           'Q(SASA)',
                           'Exposed']]
    df_SASA_L = df_SASA_L[['jobname',
                           'Chain',
                           'L_CloneID',
                           'L_SeqID',
                           'L_v_identity',
                           'Position_Insertion',
                           'Q(SASA)',
                           'Exposed']]

    # if save arg is set to True then the following code will save the DataFrame
    # as a .csv file in SASA_results
    if save == True:
        SASA_dir = os.getcwd() + os.sep + 'results' + os.sep + 'SASA_results'
        filepath_H = SASA_dir + os.sep + patient_id + '_SASA_meta_H.csv'
        df_SASA_H.to_csv(filepath_H, index=False)
        filepath_L = SASA_dir + os.sep + patient_id + '_SASA_meta_L.csv'
        df_SASA_L.to_csv(filepath_L, index=False)

    return df_SASA_H, df_SASA_L


def create_SASA_meta_dataframe_all_patients(patient_id_list, save=False):

    # creates empty lists to store dataframes
    df_H_list = []
    df_L_list = []

    # loops through patient IDs and creates per residue SASA dataframe for each
    # patient and appends dataframe to df_list
    for patient_id in patient_id_list:
        df_H = create_SASA_meta_dataframe(patient_id)[0]
        df_H['patient_id'] = patient_id
        df_H_list.append(df_H)
        df_L = create_SASA_meta_dataframe(patient_id)[1]
        df_L['patient_id'] = patient_id
        df_L_list.append(df_L)

    # concatenates dataframes in df_list to single dataframe
    df_SASA_H_all_patients = pd.concat(df_H_list)
    df_SASA_L_all_patients = pd.concat(df_L_list)

    # if save arg is set to True then the following code will save the DataFrame
    # as a .csv file in SASA_results
    if save == True:
        SASA_dir = os.getcwd() + os.sep + 'results' + os.sep + 'SASA_results'
        filepath_H = SASA_dir + os.sep + 'all_patients_SASA_meta_H.csv'
        df_SASA_H_all_patients.to_csv(filepath_H, index=False)
        filepath_L = SASA_dir + os.sep + 'all_patients_SASA_meta_L.csv'
        df_SASA_L_all_patients.to_csv(filepath_L, index=False)

    return df_SASA_H_all_patients, df_SASA_L_all_patients
