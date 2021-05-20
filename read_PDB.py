import os
import glob
import pandas as pd
from biopandas.pdb import PandasPdb
import seaborn as sns
import matplotlib.pyplot as plt
import shutil
from scipy.stats import shapiro


pd.options.display.html.table_schema = True
pd.options.display.max_rows = None
pd.options.display.max_columns = None


sns.set_theme()
sns.set_context('paper')
sns.set_palette('crest')


patient_id_list = ['YF189',
                   'YF191',
                   'YF192',
                   'YF199',
                   'YF200',
                   'YF203']


def get_singleChainPDBs(patient_id, save=False):

    # creates empty list to store jobnames of .pdb files which do not contain
    # info on both H and L chains
    singleChainPDBs = []

    # gets absolute filepaths of .pdb files
    patient_dir = os.getcwd() + os.sep + 'results' + os.sep + \
                  'ABodyBuilder_results' + os.sep + patient_id
    pdb_dir = patient_dir + os.sep + 'pdb_structures'
    pdb_filepaths = glob.glob(pdb_dir + os.sep + '*.pdb')

    # loops through .pdb files to find those with single chain
    for pdb_filepath in pdb_filepaths:
        pdb_filename = pdb_filepath.split(os.sep)[-1]
        jobname = pdb_filename.split('.')[0]
        ppdb = PandasPdb().read_pdb(pdb_filepath)
        df = ppdb.df['ATOM']
        if len(df.chain_id.unique()) != 2:
            singleChainPDBs.append(jobname)

    # saves singleChainPDBs list to .txt file in pdb_dir
    if save == True:
        with open(pdb_dir + os.sep + 'singleChainPDBs.txt', 'w') as output:
            output.write(str(singleChainPDBs))

    return singleChainPDBs


def create_modelsCDR3Length_DataFrame(patient_id):

    # creates dictionary to store CDR lengths for each model
    d = {}

    # gets absolute filepaths of .pdb files
    patient_dir = os.getcwd() + os.sep + 'results' + os.sep + \
                  'ABodyBuilder_results' + os.sep + patient_id
    pdb_dir = patient_dir + os.sep + 'pdb_structures'
    pdb_filepaths = glob.glob(pdb_dir + os.sep + '*.pdb')

    # loops through .pdb files and converts info to dataframe
    for pdb_filepath in pdb_filepaths:
        pdb_filename = pdb_filepath.split(os.sep)[-1]
        jobname = pdb_filename.split('.')[0]
        d[jobname] = {'CDRH3_len': 0, 'CDRL3_len': 0}
        ppdb = PandasPdb().read_pdb(pdb_filepath)
        df = ppdb.df['ATOM']
        df['Position_Insertion'] = \
            df['residue_number'].astype(str) + df['insertion']
        for chain_id in ['H', 'L']:
            df_chain = df[df.chain_id == chain_id]
            df_chain = df_chain[['residue_number',
                     'insertion',
                     'residue_name',
                     'chain_id',
                     'Position_Insertion']]
            df_chain = df_chain.groupby('Position_Insertion', as_index=False) \
                               .agg('first')
            df_chain = df_chain.sort_values('residue_number')
            df_chain = df_chain[(df_chain.residue_number > 103) & \
                                (df_chain.residue_number < 121)]
            CDR3_len = len(df_chain)
            if chain_id == 'H':
                d[jobname]['CDRH3_len'] = CDR3_len
            elif chain_id == 'L':
                d[jobname]['CDRL3_len'] = CDR3_len

    df = pd.DataFrame.from_dict(d, orient='index')
    df['CDRH3_len'] = df['CDRH3_len'].astype(int)
    df['CDRL3_len'] = df['CDRL3_len'].astype(int)
    df.sort_values(['CDRH3_len', 'CDRL3_len'])
    df = df[df.CDRH3_len != 0]

    df = df.reset_index(drop=False).rename(columns={'index': 'jobname'})

    return df


def plot_CDRH3Length_kde(patient_id_list, save=False):

    # creates empty list to store DataFrames
    dfs = []

    # loops through patients and generates CDR3 length DataFrame for each id
    for patient_id in patient_id_list:

        # gets DataFrame of CDR3 lengths
        df = create_modelsCDR3Length_DataFrame(patient_id)

        # adds DataFrame to list
        dfs.append(df)

    # concatenates DataFrames to single DataFrame of all patients
    df = pd.concat(dfs)

    # creates figure for plot
    fig = plt.figure(figsize=(10,8))

    # draws kernel density estimation plot on figure
    sns.displot(data=df,
                x='CDRH3_len',
                kind='kde',
                bw_adjust=1,
                palette='crest')

    # sets figure attributes
    title = '''Kernel Density Estimation of CDRH3 Lengths of All Patients'''
    plt.title(title, size=18)
    xlabel = 'Length of CDRH3 Region'
    plt.xlabel(xlabel, size=14)
    ylabel = 'Density'
    plt.ylabel(ylabel, size=14)
    plt.xlim([0, 35])
    plt.ylim([0, 0.15])

    # saves figure to plots directory if save arg is True
    if save == True:

        # creates directory for plot if it doesn't already exist
        plot_dir = os.getcwd() + os.sep + 'results' + os.sep + \
                   'plots' + os.sep + 'CDR_plots'
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        # defines name of plot file
        plot_filename = 'CDRH3Length_KDE.pdf'

        # saves figure as pdf
        fig.savefig(plot_dir + os.sep + plot_filename, bbox_inches='tight')


plot_CDRH3Length_kde(patient_id_list)


def find_ShapiroResult(df, column):
    '''Use Shapiro-Wilk test from scipy.stats to test the null hypothesis that
    data was drawn from a normal distribution, p-value < 0.05 means not enough
    confidence.'''

    results = shapiro(df[column])

    return results


def find_SeriesQuartiles(patient_id_list):

    dfs = []

    for patient_id in patient_id_list:

        df = create_modelsCDR3Length_DataFrame(patient_id)
        dfs.append(df)

    df = pd.concat(dfs)

    lower_quartile = df.CDRH3_len.quantile(0.25)
    middle_quartile = df.CDRH3_len.quantile(0.5)
    upper_quartile = df.CDRH3_len.quantile(0.75)

    quartiles = [lower_quartile, middle_quartile, upper_quartile]

    return quartiles


def create_modelsBinnedOnCDRH3Length_DataFrame(patient_id, patient_id_list):
    '''Takes patient_id string as arg, generates a DataFrame of CDR3 lengths for
    each model, determines quartiles for CDRH3 lengths and bins each model into
    3 catagories; short, medium and long CDRH3 regions.'''

    df = create_modelsCDR3Length_DataFrame(patient_id)

    quartiles = find_SeriesQuartiles(patient_id_list)

    lower_quartile = quartiles[0]
    upper_quartile = quartiles[2]

    # adds CDR_bin column to dataframe
    df['CDR_bin'] = 'M'
    df.loc[df.CDRH3_len < lower_quartile, 'CDR_bin'] = 'S'
    df.loc[df.CDRH3_len > upper_quartile, 'CDR_bin'] = 'L'

    return df


def copy_modelsOnCDRH3Length(patient_id, patient_id_list):

    # gets binned CDR dataframe
    df = create_modelsBinnedOnCDRH3Length_DataFrame(patient_id, patient_id_list)

    # defines PDB directory
    pdbDir = os.getcwd() + os.sep + 'results' + os.sep + \
             'ABodyBuilder_results' + os.sep + patient_id + os.sep + \
             'pdb_structures'

    # creates directories for binned PDBs
    PDBencodeDir = os.getcwd() + os.sep + 'results' + os.sep + \
                   'PDBencode' + os.sep + patient_id
    shortDir = PDBencodeDir + os.sep + 'shortCDR'
    mediumDir = PDBencodeDir + os.sep + 'mediumCDR'
    longDir = PDBencodeDir + os.sep + 'longCDR'
    if not os.path.exists(shortDir):
        os.makedirs(shortDir)
    if not os.path.exists(mediumDir):
        os.makedirs(mediumDir)
    if not os.path.exists(longDir):
        os.makedirs(longDir)

    # loops through rows in dataframe and copies PDB to correct dir
    for row in df.itertuples():
        pdb_filepath = pdbDir + os.sep + row.jobname + '.pdb'
        if row.CDR_bin == 'S':
            shutil.copy(pdb_filepath, shortDir)
        elif row.CDR_bin == 'M':
            shutil.copy(pdb_filepath, mediumDir)
        elif row.CDR_bin == 'L':
            shutil.copy(pdb_filepath, longDir)


def plot_CDR3Length_hist(patient_id, patient_id_list, save=False):

    # gets dataframe of CDR3 lengths
    df = create_modelsCDR3Length_DataFrame(patient_id)

    # gets bin limits
    quartiles = find_SeriesQuartiles(patient_id_list)
    lower_quartile = quartiles[0]
    upper_quartile = quartiles[2]

    # creates figure for plot
    fig, axs = plt.subplots(1, 2, figsize=(15,5))

    # draws barplots on figure
    sns.histplot(ax = axs[0],
                 x = 'CDRH3_len',
                 binwidth = 3,
                 kde = True,
                 data = df,
                 palette = 'crest')
    sns.histplot(ax = axs[1],
                 x = 'CDRL3_len',
                 binwidth = 3,
                 data = df,
                 palette = 'crest')

    # sets figure attributes
    title = 'Distribution of CDR3 Lengths for Heavy and Light Chains of ' + \
            'Patient ' + patient_id
    fig.suptitle(title, size=18)
    axs[0].set_xlim([5,35])
    axs[1].set_xlim([5,35])
    xlabel = 'CDR3 Length'
    axs[0].set_xlabel(xlabel, size=14)
    axs[1].set_xlabel(xlabel, size=14)
    H_title = 'Heavy Chain'
    L_title = 'Light Chain'
    axs[0].set_title(H_title, size=14)
    axs[1].set_title(L_title, size=14)

    # adds vertical lines on Heavy Chain subplot at bin limits
    axs[0].axvline(lower_quartile, color='red')
    axs[0].axvline(upper_quartile, color='red')

    if save == True:

        # creates directory for plot if it doesn't already exist
        plot_dir = os.getcwd() + os.sep + 'results' + os.sep + \
                    'ABodyBuilder_results' + os.sep + 'CDR_plots'
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        # defines name of plot file
        plot_filename = patient_id + '_CDR3Length_Hist.pdf'

        # saves figure as pdf
        fig.savefig(plot_dir + os.sep + plot_filename)


print(get_singleChainPDBs('YF189'))
