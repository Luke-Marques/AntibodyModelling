import pandas as pd
import glob
import os
from Levenshtein import distance as lev
from Bio.SeqIO.FastaIO import SimpleFastaParser
import clean_repertoire_data as crd


pd.options.display.html.table_schema = True
pd.options.display.max_rows = None
pd.options.display.max_columns = None


patient_id_list = ['YF189',
                   'YF191',
                   'YF192',
                   'YF199',
                   'YF200',
                   'YF203']


def read_fastaFile_to_DataFrames(patient_id):

    # defines PDBencode directory
    PDBencodeDir = os.getcwd() + os.sep + 'results' + os.sep + \
                   'PDBencode' + os.sep + patient_id + os.sep + \
                   'PDBencodeOutFiles'

    # creates list to store dataframes
    dfs = []

    # loops through CDR length bins
    for bin in 'shortCDR', 'mediumCDR', 'longCDR':

        # defines bin directory
        binDir = PDBencodeDir + os.sep + bin

        # get list of .fasta files
        fastaFiles = glob.glob(binDir + os.sep + '*.fasta')

        jobnames = []
        chains = []
        PDBencodes = []

        # loops through fasta files in
        for fastaFile in fastaFiles:

            # opens fastaFile
            with open(fastaFile) as file:

                # parses fasta file to obtain jobname, Chains and PDBencode strings
                for title, sequence in SimpleFastaParser(file):
                    jobnames.append(title.split("|")[0])
                    chains.append(title.split("|")[1])
                    PDBencodes.append(sequence)

        # creates dataframe
        data = {'jobname': jobnames,
                'chain': chains,
                'code': PDBencodes,
                'CDR_bin': bin}
        df = pd.DataFrame(data, columns=['jobname', 'chain', 'code', 'CDR_bin'])
        df = df.replace({'shortCDR': 'S',
                         'mediumCDR': 'M',
                         'longCDR': 'L'})

        # appends dataframe to dataframe list
        dfs.append(df)

    df_allBins = pd.concat(dfs)

    return df_allBins


def write_DataFrame_to_csv(df, dir_path, file_name):
    '''Takes pandas dataframe object as argument and writes to .csv file in
    directory'''

    file_path = dir_path + os.sep + file_name
    df.to_csv(file_path, index=False)


def write_fastaFiles_to_csv(patient_id_list):

    PDBencodeDir = os.getcwd() + os.sep + 'results' + os.sep + 'PDBencode'
    dfs = []
    for patient_id in patient_id_list:
        df = read_fastaFile_to_DataFrames(patient_id)
        dfs.append(df)
    df_allPatients = pd.concat(dfs)
    filename = 'AllPatients_BinnedCDRH3_PDBencodeStrings.csv'
    write_DataFrame_to_csv(df_allPatients, PDBencodeDir, filename)


def get_jobnamePairs(patient_id, CDRbin):

    # defines PDBencode directory
    PDBencodeDir = os.getcwd() + os.sep + 'results' + os.sep + \
                   'PDBencode' + os.sep + 'All_Patients' + os.sep + \
                   'PDBencodeOutFiles' + os.sep + CDRbin

    # gets list of all .fasta files in directory
    filepaths = glob.glob(PDBencodeDir + os.sep + '*.fasta')

    pairs = []

    for file1 in filepaths:
        filename1 = file1.split(os.sep)[-1].split('.')[0]
        jobname1 = filename1.split('_')[0] + '_' + filename1.split('_')[1]
        for file2 in filepaths:
            filename2 = file2.split(os.sep)[-1].split('.')[0]
            jobname2 = filename2.split('_')[0] + '_' + filename2.split('_')[1]
            pair = [jobname1, jobname2]
            pairs.append(pair)

    return pairs


def create_levDistance_dataframe(patient_id_list, chain, save=False):

    # gets PDBencode dataframe
    PDBencodeDir = os.getcwd() + os.sep + 'results' + os.sep + 'PDBencode'
    filepath = PDBencodeDir + os.sep + \
               'AllPatients_BinnedCDRH3_PDBencodeStrings.csv'
    df_PDBencode = pd.read_csv(filepath)

    df_allBins_list = []

    for patient_id in patient_id_list:

        dfs = []

        for CDRbin in 'shortCDR', 'mediumCDR', 'longCDR':

            jobname1_list = []
            jobname2_list = []
            levDistance_list = []

            pairs = get_jobnamePairs(patient_id, CDRbin)

            for pair in pairs:

                jobname1_list.append(pair[0])
                jobname2_list.append(pair[1])

                df_cond1 = df_PDBencode[(df_PDBencode.jobname==pair[0]) & \
                                        (df_PDBencode.chain==chain)]
                code1 = df_cond1.code.values[0]

                df_cond2 = df_PDBencode[(df_PDBencode.jobname==pair[1]) & \
                                        (df_PDBencode.chain==chain)]
                code2 = df_cond2.code.values[0]

                levDistance = lev(code1, code2)
                levDistance_list.append(levDistance)

            # creates dataframe
            data = {'jobname1': jobname1_list,
                    'jobname2': jobname2_list,
                    'CDR_bin': CDRbin,
                    'conf_levD': levDistance_list}
            df = pd.DataFrame(data, columns=['jobname1',
                                             'jobname2',
                                             'CDR_bin',
                                             'conf_levD'])
            df = df.replace({'shortCDR': 'S',
                             'mediumCDR': 'M',
                             'longCDR': 'L'})

            dfs.append(df)

        df_allBins = pd.concat(dfs)
        df_allBins_list.append(df_allBins)

    dff = pd.concat(df_allBins_list)

    dff['patient_id1'] = dff['jobname1'].str.split('D').str[0]
    dff['patient_id2'] = dff['jobname2'].str.split('D').str[0]

    youngList = ['YF199',
                 'YF200',
                 'YF203']
    dff['age1'] = 'Old'
    dff['age2'] = 'Old'
    dff.loc[dff.jobname1.isin(youngList), 'age1'] = 'Young'
    dff.loc[dff.jobname2.isin(youngList), 'age2'] = 'Young'

    jobnames = dff.jobname1.unique()

    dfS_list = []
    dfR_list = []

    for patient_id in patient_id_list:
        dfS = crd.reshape_paired_cell_dataframe(patient_id)
        dfS_list.append(dfS)

    dfS = pd.concat(dfS_list)
    dfS = dfS[dfS.jobname.isin(jobnames)]

    if chain == 'H':
        Seq_IDs = dfS.H_SeqID.unique()
    elif chain == 'L':
        Seq_IDs = dfS.L_SeqID.unique()

    for patient_id in patient_id_list:
        file_path = os.getcwd() + os.sep + 'inputs' + os.sep + \
                    'repertoire_data' + os.sep + 'Pacbio_Oct19_lib_' + \
                    patient_id + \
                    '_withClones_collapsedByCell_annotated_withCellType.txt'
        dfR = pd.read_csv(file_path,
                          sep='\t',
                          low_memory=False,
                          usecols=['Seq_ID',
                                   'Vfamily',
                                   'Vgene',
                                   'Class'])
        dfR.Class = dfR.Class.replace({'A': 'H',
                                       'G': 'H',
                                       'M': 'H',
                                       'L': 'L',
                                       'K': 'L',})
        dfR = dfR[dfR.Class == chain]
        dfR['jobname'] = dfR['Seq_ID'].str.split('_').str[0] + '_' + \
                         dfR['Seq_ID'].str.split('_').str[3]
        dfR = dfR[dfR.Seq_ID.isin(Seq_IDs)]
        dfR = dfR[['jobname', 'Seq_ID', 'Vfamily', 'Vgene']]
        dfR_list.append(dfR)

    dfR = pd.concat(dfR_list)

    # merge to add SeqID and gene info to DataFrame
    dff = dff.merge(dfR, left_on='jobname1', right_on='jobname')
    dff = dff.rename(columns={'Seq_ID': 'Seq_ID1',
                              'Vfamily': 'Vfamily1',
                              'Vgene': 'Vgene1'}) \
             .drop(columns=['jobname'])
    dff = dff.merge(dfR, left_on='jobname2', right_on='jobname')
    dff = dff.rename(columns={'Seq_ID': 'Seq_ID2',
                              'Vfamily': 'Vfamily2',
                               'Vgene': 'Vgene2'}) \
             .drop(columns=['jobname'])

    # merge to add amino-acid sequence to DataFrame
    aa_seq_col = chain + '_Seq'
    dfS = dfS[['jobname', aa_seq_col]]
    dff = dff.merge(dfS, left_on='jobname1', right_on='jobname')
    dff = dff.rename(columns={aa_seq_col : 'aa_sequence1'}) \
             .drop(columns=['jobname'])
    dff = dff.merge(dfS, left_on='jobname2', right_on='jobname')
    dff = dff.rename(columns={aa_seq_col : 'aa_sequence2'}) \
             .drop(columns=['jobname'])

    # merge to add conformational code to DataFrame
    df_PDBencode = df_PDBencode[df_PDBencode.chain == chain]
    df_PDBencode = df_PDBencode[['jobname', 'code']]
    dff = dff.merge(df_PDBencode, left_on='jobname1', right_on='jobname')
    dff = dff.rename(columns={'code': 'conf_code1'}) \
             .drop(columns=['jobname'])
    dff = dff.merge(df_PDBencode, left_on='jobname2', right_on='jobname')
    dff = dff.rename(columns={'code': 'conf_code2'}) \
             .drop(columns=['jobname'])

    # calculate Levenshtein distance between aa sequences for each row
    levDistance_list = []
    for row in dff.itertuples():
        levDistance = lev(row.aa_sequence1, row.aa_sequence2)
        levDistance_list.append(levDistance)
    dff['aa_levD'] = levDistance_list

    dff = dff.drop_duplicates()

    # saves DataFrame if save arg is True
    if save == True:

        # defines name of plot file
        filename = 'CDR' + chain + '3_levDistance.csv'

        # saves figure as pdf
        dff.to_csv(PDBencodeDir + os.sep + filename)


print(lev('ARDAVYDESSEGGYVRLDP', 'AIEGYCYSTSCLTHDAFDI'))
print(lev('GEIEIGALRKWMCXSRCKMNBI', 'GEIEIGALNKVMCXSRCKMNBI'))
