import pandas as pd
import os
import numpy as np
import clean_repertoire_data as crd


pd.options.display.html.table_schema = True
pd.options.display.max_rows = None
pd.options.display.max_columns = None


def create_RegionAASeq_DataFrame(patient_id):
    '''Takes the patient_id string (e.g. 'YF189') as an arg and converts that
    patient's sequence data file to a dataframe and returns the dataframe
    (could be modified to take full name of data file as an input instead)'''


    # can change to name of your sequence alignment data file
    file_path = os.getcwd() + os.sep + 'inputs' + os.sep + 'sequence_data' \
                + os.sep + patient_id + '_sequence_data.tsv'

    # reads file and converts to dataframe
    df = pd.read_csv(file_path,
                     sep='\t',
                     low_memory=False,
                     usecols=['sequence_id', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa',
                              'fwr1_aa', 'fwr2_aa', 'fwr3_aa', 'fwr4_aa',
                              'vj_in_frame', 'productive',
                              'sequence_alignment_aa'])

    # filters DataFrame to only include sequences where the vj gene segments are
    # in frame and thus we have the complete sequence for them
    df = df[df.vj_in_frame == 'T']
    df = df[df.productive == 'T']
    df = df[~df.sequence_alignment_aa.str.contains('*', regex=False)]
    df = df.drop(columns=['vj_in_frame', 'productive', 'sequence_alignment_aa'])

    return df


def find_RegionAALengths(patient_id):
    '''Takes a patient_id string as an arg and generates a DataFrame of amino-
    acid sequences of all regions, using the create_RegionAASeq_DataFrame
    function, and returns a the DataFrame with additional Region length
    columns appended.'''

    # creates the sequence DataFrame by reading .tsv file
    df = create_RegionAASeq_DataFrame(patient_id)

    # gets list of column names (Region names) in DataFrame
    Regions = list(df)

    # loops through each Region in Regions list
    for Region in Regions:

        # creates new column for storing Region length
        col_name = Region + '_Length'
        df[col_name] = df[Region].str.len()

    df = df.drop(columns=['sequence_id_Length'])

    return df


def merge_RegionAA_Repertoire(patient_id):

    # reads repertoire data .csv file to DataFrame
    file_path = os.getcwd() + os.sep + 'inputs' + os.sep + 'repertoire_data' \
                + os.sep + 'Pacbio_Oct19_lib_' + patient_id \
                + '_withClones_collapsedByCell_annotated_withCellType.txt'
    dfR = pd.read_csv(file_path,
                      sep='\t',
                      low_memory=False,
                      usecols=['Seq_ID', 'Class', 'UseAsRef'])
    dfR = dfR[dfR.UseAsRef == True]
    dfR = dfR[~dfR.Seq_ID.str.contains('UN', regex=False)]
    dfR = dfR[~dfR.Seq_ID.str.contains('NA', regex=False)]
    dfR = dfR.drop(columns=['UseAsRef'])

    # gets RegionAALengths DataFrame
    dfL = find_RegionAALengths(patient_id)

    # merges dfR and dfL on Seq_ID column
    df = dfR.merge(dfL, left_on='Seq_ID', right_on='sequence_id') \
            .drop(columns=['sequence_id'])

    # creates new column 'Chain' which is 'H' for heavy chain sequences and 'L'
    # for light chain sequences
    df['Chain'] = df.Class
    df.Chain = df.Chain.replace({'A': 'H',
                                 'G': 'H',
                                 'M': 'H',
                                 'K': 'L'})

    return df


def create_RegionAALength_DataFrame(patient_id):

    df = merge_RegionAA_Repertoire(patient_id)

    Regions = ['cdr1', 'cdr2', 'cdr3', 'fwr1', 'fwr2', 'fwr3', 'fwr4']

    for Region in Regions:
        for Chain in 'H', 'L':
            col_name = Region + '_' + Chain + '_aa_Length'
            df[col_name] = np.nan
            df.loc[df.Chain == Chain, col_name] = df[Region + '_aa_Length']

    return df


def find_avgCDR3LengthPerClone(patient_id):

    df_agg = crd.aggregate_dataframe(patient_id)
    df_agg = df_agg.sort_values('no_seqs', ascending=False)

    df_Clones = crd.prepare_repertoire_dataframe(patient_id)
    df_Clones = df_Clones[['CloneID', 'SeqID', 'Class']]
    df_Clones['Chain'] = df_Clones.Class
    df_Clones.Chain = df_Clones.Chain.replace({'A': 'H',
                                               'G': 'H',
                                               'M': 'H',
                                               'L': 'L',
                                               'K': 'L'})

    df_RegLen = find_RegionAALengths(patient_id)

    df = df_RegLen.merge(df_Clones, left_on='sequence_id', right_on='SeqID') \
         .drop(columns='sequence_id')

    df = df[['CloneID', 'Chain', 'cdr3_aa_Length']]
    df = df[df.Chain == 'H']
    df = df.drop(columns='Chain').reset_index(drop=True)
    df = df.groupby('CloneID', as_index=False).agg({'cdr3_aa_Length': 'mean'})

    df = df.merge(df_agg, left_on='CloneID', right_on='CloneID')
    df = df.sort_values('no_seqs', ascending=False)
    df = df[['no_seqs', 'cdr3_aa_Length']]

    return df

old_list = ['YF189', 'YF191', 'YF192']
young_list = ['YF199', 'YF200', 'YF203']

dfOs = []
dfYs = []
for patient_id in old_list:
    df = create_RegionAALength_DataFrame(patient_id)
    df = df[df.Class.isin(['A', 'G', 'M'])]
    dfOs.append(df)
for patient_id in young_list:
    df = create_RegionAALength_DataFrame(patient_id)
    df = df[df.Class.isin(['A', 'G', 'M'])]
    dfYs.append(df)
dfO = pd.concat(dfOs)
dfY = pd.concat(dfYs)
num_old_seqs = len(dfO)
num_young_seqs = len(dfY)
avg_CDRH3Length_old = dfO.cdr3_aa_Length.mean()
avg_CDRH3Length_young = dfY.cdr3_aa_Length.mean()
std_CDRH3Length_old = dfO.cdr3_aa_Length.std()
std_CDRH3Length_young = dfY.cdr3_aa_Length.std()

print('OLD PATIENTS')
print('Number of Sequences: ' + str(num_old_seqs))
print('Mean Length: ' + str(avg_CDRH3Length_old))
print('STD Length: ' + str(std_CDRH3Length_old))
print('YOUNG PATIENTS')
print('Number of Sequences: ' + str(num_young_seqs))
print('Mean Length: ' + str(avg_CDRH3Length_young))
print('STD Length: ' + str(std_CDRH3Length_young))

old_std_div_N = (std_CDRH3Length_old**2)/num_old_seqs
young_std_div_N = (std_CDRH3Length_young**2)/num_young_seqs
SE = np.sqrt(old_std_div_N + young_std_div_N)

print('SE: ' + str(SE))

Z = (abs(avg_CDRH3Length_old - avg_CDRH3Length_young))/SE

print('Z: ' + str(Z))
