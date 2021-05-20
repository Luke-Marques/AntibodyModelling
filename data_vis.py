import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme(color_codes=True)
from scipy.cluster import hierarchy
from fastcluster import linkage
from matplotlib.colors import rgb2hex, colorConverter
from collections import Counter
import find_CDRlengths as fl

%matplotlib inline

pd.options.display.html.table_schema = True
pd.options.display.max_rows = None
pd.options.display.max_columns = None


patient_id_list = ['YF189', 'YF191', 'YF192', 'YF199', 'YF200', 'YF203']


def draw_confCmap(CDRbin, save=False):

    # reads csv of CDRH3 PDBencode Levenshtein Distances to DataFrame
    df = pd.read_csv(os.getcwd() + os.sep + 'results' + os.sep + 'PDBencode' \
                     + os.sep + 'CDRH3_levDistance.csv')

    # slices DataFrame on CDRbin
    df = df[df.CDR_bin == CDRbin].reset_index(drop=True)

    # pivots DataFrame so that values are Levenshtein Distances
    df_pivot = pd.pivot_table(df,
                              index = ['jobname1'],
                              columns = ['jobname2'],
                              values = 'conf_levD')

    # resets pivoted DataFrame to unstacked index
    df_pivot_original = df_pivot
    df_pivot.columns.name = None
    df_pivot = df_pivot.reset_index().rename(columns={'jobname1': 'jobname'})

    # gets repertoire DataFrames TEMP -> ADD DATA TO levDistance DATAFRAME
    patient_id_list = ['YF189',
                       'YF191',
                       'YF192',
                       'YF199',
                       'YF200',
                       'YF203']
    dfR_list = []
    for patient_id in patient_id_list:
        file_path = os.getcwd() + os.sep + 'inputs' + os.sep + \
                    'repertoire_data' + os.sep + 'Pacbio_Oct19_lib_' + \
                    patient_id + \
                    '_withClones_collapsedByCell_annotated_withCellType.txt'
        dfR = pd.read_csv(file_path,
                          sep='\t',
                          low_memory=False,
                          usecols=['Seq_ID', 'Jfamily', 'Dfamily'])
        dfR_list.append(dfR)
    dfR = pd.concat(dfR_list)

    # creates patient_id column on df_pivot
    df_pivot['patient_id'] = df_pivot['jobname'].str.split('D').str[0]

    # # creates age row_colors column on df_pivot
    # youngList = ['YF199', 'YF200', 'YF203']
    # df_pivot['age'] = 'Old'
    # df_pivot.loc[df_pivot.patient_id.isin(youngList), 'age'] = 'Young'

    # creates Jfamily column on df_pivot
    Vfamily_list = []
    for row in df_pivot.itertuples():
        dff = df[df.jobname1 == row.jobname].head(1)
        Vfamily = dff.iloc[0]['Vfamily1']
        Vfamily_list.append(Vfamily)
    df_pivot['Vfamily'] = Vfamily_list

    # creates Jfamily and Dfamily columns on df_pivot
    Seq_ID_list = []
    for row in df_pivot.itertuples():
        dff = df[df.jobname1 == row.jobname].head(1)
        Seq_ID = dff.iloc[0]['Seq_ID1']
        Seq_ID_list.append(Seq_ID)
    df_pivot['Seq_ID'] = Seq_ID_list
    Jfamily_list = []
    Dfamily_list = []
    for row in df_pivot.itertuples():
        dff = dfR[dfR.Seq_ID == row.Seq_ID].head(1)
        Jfamily = dff.iloc[0]['Jfamily']
        Jfamily_list.append(Jfamily)
        Dfamily = dff.iloc[0]['Dfamily']
        Dfamily_list.append(Dfamily)
    df_pivot['Jfamily'] = Jfamily_list
    df_pivot['Dfamily'] = Dfamily_list

    # creates empty list to store col_colors maps
    col_colors = []

    # creates patient_id row_colors palette
    num_patient_id_uniqVals= len(df_pivot['patient_id'].unique())
    old_pal = list(sns.color_palette('Greens', 3))
    young_pal = list(sns.color_palette('Purples', 3))
    patient_id_pal = old_pal + young_pal
    patient_id_labels = df_pivot.pop('patient_id')
    patient_id_lut = dict(zip(patient_id_labels.unique(), patient_id_pal))
    patient_id_colors = patient_id_labels.map(patient_id_lut)
    col_colors.append(patient_id_colors)

    display(sns.color_palette('Greens', 3))
    display(sns.color_palette('Purples', 3))
    print(patient_id_labels.unique())

    # # creates age row_colors palette
    # num_age_uniqVals= len(df_pivot['age'].unique())
    # age_pal = list(sns.color_palette('Paired_r', num_age_uniqVals))
    # age_labels = df_pivot.pop('age')
    # age_lut = dict(zip(age_labels.unique(), age_pal))
    # age_colors = age_labels.map(age_lut)
    # col_colors.append(age_colors)

    # creates Vfamily row_colors palette
    num_Vfamily_uniqVals= len(df_pivot['Vfamily'].unique())
    Vfamily_pal = list(sns.color_palette('Blues', num_Vfamily_uniqVals))
    Vfamily_labels = df_pivot.pop('Vfamily')
    Vfamily_lut = dict(zip(Vfamily_labels.unique(), Vfamily_pal))
    Vfamily_colors = Vfamily_labels.map(Vfamily_lut)
    col_colors.append(Vfamily_colors)

    display(sns.color_palette('Blues', num_Vfamily_uniqVals))
    print(Vfamily_labels.unique())

    # creates Jfamily row_colors palette
    num_Jfamily_uniqVals= len(df_pivot['Jfamily'].unique())
    Jfamily_pal = list(sns.color_palette('Oranges', num_Jfamily_uniqVals))
    Jfamily_labels = df_pivot.pop('Jfamily')
    Jfamily_lut = dict(zip(Jfamily_labels.unique(), Jfamily_pal))
    Jfamily_colors =Jfamily_labels.map(Jfamily_lut)
    col_colors.append(Jfamily_colors)

    display(sns.color_palette('Oranges', num_Jfamily_uniqVals))
    print(Jfamily_labels.unique())

    # creates Dfamily row_colors palette
    num_Dfamily_uniqVals= len(df_pivot['Dfamily'].unique())
    Dfamily_pal = list(sns.color_palette('Reds', num_Dfamily_uniqVals))
    Dfamily_labels = df_pivot.pop('Dfamily')
    Dfamily_lut = dict(zip(Dfamily_labels.unique(), Dfamily_pal))
    Dfamily_colors =Dfamily_labels.map(Dfamily_lut)
    col_colors.append(Dfamily_colors)

    display(sns.color_palette('Reds', num_Dfamily_uniqVals))
    print(Dfamily_labels.unique())

    # drops columns
    df_pivot = df_pivot.drop(columns=['Seq_ID'])

    # sets index and removes names of columns and rows
    df_pivot = df_pivot.set_index('jobname')
    df_pivot.index.name = None
    df_pivot = df_pivot.reset_index(drop=True)
    num_cols = len(df_pivot.columns)
    colNames = list(range(0, num_cols, 1))
    df_pivot.columns = colNames

    # normalise df_pivot values
    levD_max = df_pivot.max().max()
    df_pivot = df_pivot.div(levD_max)

    # # gets clustering linkage for df_pivot
    # link = hierarchy.linkage(df_pivot, method='ward', metric='euclidean')
    # print(link)

    # draws clusermap from
    g = sns.clustermap(df_pivot,
                       method='ward',
                       metric='euclidean',
                       cmap='mako',
                       xticklabels=False,
                       yticklabels=False,
                       col_colors=col_colors)

    # den = hierarchy.dendrogram(link,
    #                            labels=df_pivot.index,
    #                            abv_threshold_color='#AAAAAA')

    # removes row dendrogram whilst maintaining row clustering
    g.ax_row_dendrogram.set_visible(False)

    if save == True:

        plot_path = os.getcwd() + os.sep + CDRbin + '_CDRH3_conf_clustermap.png'
        g.savefig(plot_path)


draw_confCmap('S', save=True)
draw_confCmap('M', save=True)
draw_confCmap('L', save=True)


def draw_seqCmap(CDRbin, save=False):

    # reads csv of CDRH3 PDBencode Levenshtein Distances to DataFrame
    df = pd.read_csv(os.getcwd() + os.sep + 'results' + os.sep + 'PDBencode' \
                     + os.sep + 'CDRH3_levDistance.csv')

    # slices DataFrame on CDRbin
    df = df[df.CDR_bin == CDRbin].reset_index(drop=True)

    # pivots DataFrame so that values are Levenshtein Distances
    df_pivot = pd.pivot_table(df,
                              index = ['jobname1'],
                              columns = ['jobname2'],
                              values = 'aa_levD')

    # resets pivoted DataFrame to unstacked index
    df_pivot_original = df_pivot
    df_pivot.columns.name = None
    df_pivot = df_pivot.reset_index().rename(columns={'jobname1': 'jobname'})

    # gets repertoire DataFrames TEMP -> ADD DATA TO levDistance DATAFRAME
    patient_id_list = ['YF189',
                       'YF191',
                       'YF192',
                       'YF199',
                       'YF200',
                       'YF203']
    dfR_list = []
    for patient_id in patient_id_list:
        file_path = os.getcwd() + os.sep + 'inputs' + os.sep + \
                    'repertoire_data' + os.sep + 'Pacbio_Oct19_lib_' + \
                    patient_id + \
                    '_withClones_collapsedByCell_annotated_withCellType.txt'
        dfR = pd.read_csv(file_path,
                          sep='\t',
                          low_memory=False,
                          usecols=['Seq_ID', 'Jfamily', 'Dfamily'])
        dfR_list.append(dfR)
    dfR = pd.concat(dfR_list)

    # creates patient_id column on df_pivot
    df_pivot['patient_id'] = df_pivot['jobname'].str.split('D').str[0]

    # # creates age row_colors column on df_pivot
    # youngList = ['YF199', 'YF200', 'YF203']
    # df_pivot['age'] = 'Old'
    # df_pivot.loc[df_pivot.patient_id.isin(youngList), 'age'] = 'Young'

    # creates Jfamily column on df_pivot
    Vfamily_list = []
    for row in df_pivot.itertuples():
        dff = df[df.jobname1 == row.jobname].head(1)
        Vfamily = dff.iloc[0]['Vfamily1']
        Vfamily_list.append(Vfamily)
    df_pivot['Vfamily'] = Vfamily_list

    # creates Jfamily and Dfamily columns on df_pivot
    Seq_ID_list = []
    for row in df_pivot.itertuples():
        dff = df[df.jobname1 == row.jobname].head(1)
        Seq_ID = dff.iloc[0]['Seq_ID1']
        Seq_ID_list.append(Seq_ID)
    df_pivot['Seq_ID'] = Seq_ID_list
    Jfamily_list = []
    Dfamily_list = []
    for row in df_pivot.itertuples():
        dff = dfR[dfR.Seq_ID == row.Seq_ID].head(1)
        Jfamily = dff.iloc[0]['Jfamily']
        Jfamily_list.append(Jfamily)
        Dfamily = dff.iloc[0]['Dfamily']
        Dfamily_list.append(Dfamily)
    df_pivot['Jfamily'] = Jfamily_list
    df_pivot['Dfamily'] = Dfamily_list

    # creates empty list to store col_colors maps
    col_colors = []

    # creates patient_id row_colors palette
    num_patient_id_uniqVals= len(df_pivot['patient_id'].unique())
    old_pal = list(sns.color_palette('Greens', 3))
    young_pal = list(sns.color_palette('Purples', 3))
    patient_id_pal = old_pal + young_pal
    patient_id_labels = df_pivot.pop('patient_id')
    patient_id_lut = dict(zip(patient_id_labels.unique(), patient_id_pal))
    patient_id_colors = patient_id_labels.map(patient_id_lut)
    col_colors.append(patient_id_colors)

    display(sns.color_palette('Greens', 3))
    display(sns.color_palette('Purples', 3))
    print(patient_id_labels.unique())

    # # creates age row_colors palette
    # num_age_uniqVals= len(df_pivot['age'].unique())
    # age_pal = list(sns.color_palette('Paired_r', num_age_uniqVals))
    # age_labels = df_pivot.pop('age')
    # age_lut = dict(zip(age_labels.unique(), age_pal))
    # age_colors = age_labels.map(age_lut)
    # col_colors.append(age_colors)

    # creates Vfamily row_colors palette
    num_Vfamily_uniqVals= len(df_pivot['Vfamily'].unique())
    Vfamily_pal = list(sns.color_palette('Blues', num_Vfamily_uniqVals))
    Vfamily_labels = df_pivot.pop('Vfamily')
    Vfamily_lut = dict(zip(Vfamily_labels.unique(), Vfamily_pal))
    Vfamily_colors = Vfamily_labels.map(Vfamily_lut)
    col_colors.append(Vfamily_colors)

    display(sns.color_palette('Blues', num_Vfamily_uniqVals))
    print(Vfamily_labels.unique())

    # creates Jfamily row_colors palette
    num_Jfamily_uniqVals= len(df_pivot['Jfamily'].unique())
    Jfamily_pal = list(sns.color_palette('Oranges', num_Jfamily_uniqVals))
    Jfamily_labels = df_pivot.pop('Jfamily')
    Jfamily_lut = dict(zip(Jfamily_labels.unique(), Jfamily_pal))
    Jfamily_colors =Jfamily_labels.map(Jfamily_lut)
    col_colors.append(Jfamily_colors)

    display(sns.color_palette('Oranges', num_Jfamily_uniqVals))
    print(Jfamily_labels.unique())

    # creates Dfamily row_colors palette
    num_Dfamily_uniqVals= len(df_pivot['Dfamily'].unique())
    Dfamily_pal = list(sns.color_palette('Reds', num_Dfamily_uniqVals))
    Dfamily_labels = df_pivot.pop('Dfamily')
    Dfamily_lut = dict(zip(Dfamily_labels.unique(), Dfamily_pal))
    Dfamily_colors =Dfamily_labels.map(Dfamily_lut)
    col_colors.append(Dfamily_colors)

    display(sns.color_palette('Reds', num_Dfamily_uniqVals))
    print(Dfamily_labels.unique())

    # drops columns
    df_pivot = df_pivot.drop(columns=['Seq_ID'])

    # sets index and removes names of columns and rows
    df_pivot = df_pivot.set_index('jobname')
    df_pivot.index.name = None
    df_pivot = df_pivot.reset_index(drop=True)
    num_cols = len(df_pivot.columns)
    colNames = list(range(0, num_cols, 1))
    df_pivot.columns = colNames

    # normalise df_pivot values
    levD_max = df_pivot.max().max()
    df_pivot = df_pivot.div(levD_max)

    # # gets clustering linkage for df_pivot
    # link = hierarchy.linkage(df_pivot, method='ward', metric='euclidean')
    # print(link)

    # draws clusermap from
    g = sns.clustermap(df_pivot,
                       method='ward',
                       metric='euclidean',
                       cmap='mako',
                       xticklabels=False,
                       yticklabels=False,
                       col_colors=col_colors)

    # den = hierarchy.dendrogram(link,
    #                            labels=df_pivot.index,
    #                            abv_threshold_color='#AAAAAA')

    # removes row dendrogram whilst maintaining row clustering
    g.ax_row_dendrogram.set_visible(False)

    if save == True:

        plot_path = os.getcwd() + os.sep + CDRbin + '_CDRH3_aaSeq_clustermap.png'
        g.savefig(plot_path)


def draw_confDist(CDR=None, CDRbin=None, save=False):

    if CDRbin != None:
        CDR = 'CDRH3'
    # reads csv of CDRH3 PDBencode Levenshtein Distances to DataFrame
    df = pd.read_csv(os.getcwd() + os.sep + 'results' + os.sep + 'PDBencode' \
                     + os.sep + CDR + '_levDistance.csv')

    # slices DataFrame on CDRbin
    df = df[df.CDR_bin == CDRbin].reset_index(drop=True)

    # creates figure and axis
    fig, ax =  plt.subplots(figsize=(7,5))

    # draws histogram on axis
    ax = sns.histplot(data=df,
                      x=conf_levD,
                      palette='mako')

    # sets figure attributes
    ax.set_xlabel('Pairwise Conformation Code Levenshtein Distance')
    ax.set_ylabel('Count')

    if save == True:
        plotDir = os.getcwd() + os.sep + 'results' + os.sep + 'plots' + os.sep + \
                  'CDRConfDistributions'
        if not os.path.exists(plotDir):
            os.makedirs(plotDir)
        if CDRbin != None:
            plotPath = plotDir + os.sep + CDR + os.sep + CDRbin
        else:
            plotPath = plotDir + os.sep + CDR
        if not os.path.exists(plotPath):
            os.makedirs(plotPath)
        fig.savefig(plotPath + os.sep + CDR + '_' + CDRbin + '_Conf_Hist.pdf')


def count_confCodes(CDRbin):
    '''Takes levD DataFrame (either for heavy or light chain sequences) as arg
    and adds columns for each letter in the conf_code, with proportion of
    conf_code as that letter as values.'''

    # reads csv of CDRH3 PDBencode Levenshtein Distances to DataFrame
    df = pd.read_csv(os.getcwd() + os.sep + 'results' + os.sep + 'PDBencode' \
                     + os.sep + 'CDRH3_levDistance.csv')

    # filters DataFrame for CDRbin
    df = df[df.CDR_bin == CDRbin]

    # selects relevant columns
    df = df[['jobname1', 'conf_code1']] \
           .rename(columns={'jobname1': 'jobname',
                            'conf_code1': 'conf_code'})

    # gets list of jobname values from Levenshtein Dataframe
    jobname_list = df.jobname.to_list()

    # gets list of conf_code values from Levenshtein DataFrame
    conf_code_list = df.conf_code.to_list()

    # creates an empty list to store counter dictionaries
    counter_list = []

    # loops through each row of dataframe and generates Counter dictionary from
    # conformation code string
    for row in df.itertuples():
        conf_code = row.conf_code
        counter = Counter(conf_code)
        for item, count in counter.items():
            counter[item] /= len(conf_code)
        counter_list.append(counter)

    # converts each counter dict in counter_list to single pandas DataFrame\
    dfC_list = []
    for counter in counter_list:
        dfC = pd.DataFrame(counter, index=['count']).reset_index(drop=True)
        dfC_list.append(dfC)
    dfC = pd.concat(dfC_list)

    # adds jobname and conf_code columns
    dfC.index = conf_code_list
    dfC = dfC.reset_index(drop=False).rename(columns={'index': 'conf_code'})
    dfC.index = jobname_list
    dfC = dfC.reset_index(drop=False).rename(columns={'index': 'jobname'})

    dfC = dfC.fillna(0)

    dfC = dfC[['jobname', 'conf_code', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
               'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
               'V', 'W', 'X', 'Y']]

    return dfC


def draw_confCodesBoxplot(CDRbin, save=False):

    # gets conf_code count DataFrame from count_confCodes function
    df = count_confCodes(CDRbin)

    # drops unused columns
    df = df.drop(columns=['jobname', 'conf_code'])

    # converts wide-form DataFrame to long-form DataFrame
    dfM = pd.melt(df)

    # draws boxplot
    fig, ax = plt.subplots(figsize=(10,5))
    ax = sns.boxplot(x='variable',
                     y='value',
                     data=dfM,
                     color=sns.color_palette('mako', 5)[3])


def draw_confCodesBarplot(CDRbin, save=False):

    # gets conf_code count DataFrame from count_confCodes function
    df = count_confCodes(CDRbin)

    # drops unused columns
    df = df.drop(columns=['jobname', 'conf_code'])

    # converts wide-form DataFrame to long-form DataFrame
    dfM = pd.melt(df)

    # draws barplot
    fig, ax = plt.subplots(figsize=(10,5))
    ax = sns.barplot(x='variable',
                     y='value',
                     data=dfM,
                     color=sns.color_palette('mako', 5)[3])


def draw_CDRLengths_Histograms(patient_id_list, save=False):

    # defines directory to store histograms
    plotDir = os.getcwd() + os.sep + 'results' + os.sep + 'plots' + os.sep + \
              'CDR_Lengths'
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)

    # defines patient_id values in age groups
    YoungPatients = ['YF199', 'YF200', 'YF203']
    OldPatients = ['YF189', 'YF191', 'YF192']

    # defines names of Regions and stores in list
    CDRs = ['cdr1', 'cdr2', 'cdr3']

    # defines color palette to use for Region histograms
    L_colors = sns.color_palette('hls')[:3]
    H_colors = sns.color_palette('hls')[3:6]

    dfs = []
    dfOlds = []
    dfYoungs = []

    for patient_id in patient_id_list:

        # generates Region Length DataFrame for patient
        df = fl.merge_RegionAA_Repertoire(patient_id)

        # appends DataFrame to all patients DataFrame list
        dfs.append(df)

        # checks what age group patient is in and appends DataFrame to age list
        if patient_id in YoungPatients:
            dfYoungs.append(df)
        elif patient_id in OldPatients:
            dfOlds.append(df)

        for Chain in 'H', 'L':

            # filters DataFrame on Chain
            dfChain = df[df.Chain == Chain]

            for CDR in CDRs:

                # finds position of CDR in CDRs list
                index = CDRs.index(CDR)

                # gets color from colors list with position of index
                if Chain == 'H':
                    color = H_colors[index]
                elif Chain == 'L':
                    color = L_colors[index]

                # gets name of Region Length column
                col_name = CDR + '_aa_Length'

                # creates figure and axis
                fig, ax =  plt.subplots(figsize=(7,5))

                # draws histogram on axis
                ax = sns.histplot(data=dfChain,
                                   x=col_name,
                                   color=color)

                # sets figure attributes
                ax.set_xlabel('Region Amino-Acid Length')
                ax.set_ylabel('Count')
                # ax.set_ylim()

                if save == True:
                    plotPath = plotDir + os.sep + patient_id + os.sep + Chain
                    if not os.path.exists(plotPath):
                        os.makedirs(plotPath)
                    fig.savefig(plotPath + os.sep + CDR + '_Length_Hist.pdf')

        for Class in 'L', 'K':

            # filters DataFrame on Class
            dfClass = df[df.Class == Class]

            for CDR in CDRs:

                # finds position of CDR in CDRs list
                index = CDRs.index(CDR)

                # gets color from colors list with position of index
                color = L_colors[index]

                # gets name of Region Length column
                col_name = CDR + '_aa_Length'

                # creates figure and axis
                fig, ax =  plt.subplots(figsize=(7,5))

                # draws histogram on axis
                ax = sns.histplot(data=dfChain,
                                   x=col_name,
                                   color=color)

                # sets figure attributes
                ax.set_xlabel('Region Amino-Acid Length')
                ax.set_ylabel('Count')
                # ax.set_ylim()

                if save == True:
                    plotPath = plotDir + os.sep + patient_id + os.sep + 'L' \
                               + os.sep + Class
                    if not os.path.exists(plotPath):
                        os.makedirs(plotPath)
                    fig.savefig(plotPath + os.sep + CDR + '_Length_Hist.pdf')

    # concatenates DataFrames stored in lists to single DataFrames
    dfAll = pd.concat(dfs)
    dfOld = pd.concat(dfOlds)
    dfYoung = pd.concat(dfYoungs)

    # adds DataFrames to list to loop through
    df_list = [dfAll, dfOld, dfYoung]
    groups = ['AllPatients', 'Old', 'Young']

    # loops through df_list
    n = 0
    for df in df_list:

        # get group name for DataFrame
        group = groups[n]
        n += 1

        for Chain in 'H', 'L':

            # filters DataFrame on Chain
            dfChain = df[df.Chain == Chain]

            for CDR in CDRs:

                # finds position of CDR in CDRs list
                index = CDRs.index(CDR)

                # gets color from colors list with position of index
                if Chain == 'H':
                    color = H_colors[index]
                elif Chain == 'L':
                    color = L_colors[index]

                # gets name of Region Length column
                col_name = CDR + '_aa_Length'

                # creates figure and axis
                fig, ax =  plt.subplots(figsize=(7,5))

                # draws histogram on axis
                ax = sns.histplot(data=dfChain,
                                   x=col_name,
                                   color=color)

                # sets figure attributes
                ax.set_xlabel('Region Amino-Acid Length')
                ax.set_ylabel('Count')
                # ax.set_ylim()

                if save == True:
                    plotPath = plotDir + os.sep + group + os.sep + Chain
                    if not os.path.exists(plotPath):
                        os.makedirs(plotPath)
                    fig.savefig(plotPath + os.sep + CDR + '_Length_Hist.pdf')

        for Class in 'L', 'K':

            # filters DataFrame on Class
            dfClass = df[df.Class == Class]

            for CDR in CDRs:

                # finds position of CDR in CDRs list
                index = CDRs.index(CDR)

                # gets color from colors list with position of index
                color = L_colors[index]

                # gets name of Region Length column
                col_name = CDR + '_aa_Length'

                # creates figure and axis
                fig, ax =  plt.subplots(figsize=(7,5))

                # draws histogram on axis
                ax = sns.histplot(data=dfChain,
                                   x=col_name,
                                   color=color)

                # sets figure attributes
                ax.set_xlabel('Region Amino-Acid Length')
                ax.set_ylabel('Count')
                # ax.set_ylim()

                if save == True:
                    plotPath = plotDir + os.sep + group + os.sep + 'L' \
                               + os.sep + Class
                    if not os.path.exists(plotPath):
                        os.makedirs(plotPath)
                    fig.savefig(plotPath + os.sep + CDR + '_Length_Hist.pdf')


def draw_CDRLengths_Boxplots(patient_id_list, save=False):

    # defines directory to store histograms
    plotDir = os.getcwd() + os.sep + 'results' + os.sep + 'plots' + os.sep + \
              'CDR_Lengths'
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)

    # defines names of Regions and stores in list
    CDRs = ['cdr1_L', 'cdr2_L', 'cdr3_L', 'cdr1_H', 'cdr2_H', 'cdr3_H']

    # defines patient_id values in age groups
    YoungPatients = ['YF199', 'YF200', 'YF203']
    OldPatients = ['YF189', 'YF191', 'YF192']

    dfs = []
    dfOlds = []
    dfYoungs = []

    for patient_id in patient_id_list:

        # generates Region Length DataFrame for patient
        df = fl.create_RegionAALength_DataFrame(patient_id)

        # selects relevant columns from DataFrame
        col_names = []
        for CDR in CDRs:
            col_name = CDR + '_aa_Length'
            col_names.append(col_name)
        dfCDRs = df[col_names]

        # converts wide form DataFrame to long form DataFrame
        dfCDRs = pd.melt(dfCDRs).reset_index(drop=True)

        # appends DataFrame to all patients DataFrame list
        dfs.append(dfCDRs)

        # checks what age group patient is in and appends DataFrame to age list
        if patient_id in YoungPatients:
            dfYoungs.append(dfCDRs)
        elif patient_id in OldPatients:
            dfOlds.append(dfCDRs)

        # creates figure and axis
        fig, ax = plt.subplots(figsize=(10,5))

        # draws boxplots on axis
        ax = sns.boxplot(data=dfCDRs, x='variable', y='value', palette='Paired')

        # sets figure attributes
        ticklabels = ['CDRL1', 'CDRL2', 'CDRL3', 'CDRH1', 'CDRH2', 'CDRH3']
        ax.set_xticklabels(ticklabels)
        ax.set_xlabel('Complimentary Determining Region')
        ax.set_ylabel('Amino-Acid Length')
        ax.set_ylim(bottom=0, top=50)

        if save == True:
            plotPath = plotDir + os.sep + patient_id
            if not os.path.exists(plotPath):
                os.makedirs(plotPath)
            fig.savefig(plotPath + os.sep + 'CDRLengthBoxplots.pdf')

    # concatenates DataFrames stored in lists to single DataFrames
    dfAll = pd.concat(dfs)
    dfOld = pd.concat(dfOlds)
    dfYoung = pd.concat(dfYoungs)

    # adds DataFrames to list to loop through
    df_list = [dfAll, dfOld, dfYoung]
    groups = ['AllPatients', 'Old', 'Young']

    # loops through DataFrames and draws boxplots for each DataFrame
    n = 0
    for df in df_list:

        # selects correct group name from groups
        group = groups[n]
        n += 1

        # creates figure and axis
        fig, ax = plt.subplots(figsize=(10,5))

        # draws boxplots on axis
        ax = sns.boxplot(data=df, x='variable', y='value', palette='Paired')

        # sets figure attributes
        ticklabels = ['CDRL1', 'CDRL2', 'CDRL3', 'CDRH1', 'CDRH2', 'CDRH3']
        ax.set_xticklabels(ticklabels)
        ax.set_xlabel('Complimentary Determining Region')
        ax.set_ylabel('Amino-Acid Length')
        ax.set_ylim(bottom=0, top=50)

        if save == True:
            plotPath = plotDir + os.sep + group
            if not os.path.exists(plotPath):
                os.makedirs(plotPath)
            fig.savefig(plotPath + os.sep + 'CDRLengthBoxplots.pdf')


def draw_CDRH3LenCloneSize_Scatterplot(patient_id_list, save=False):

    dfs = []
    for patient_id in patient_id_list:
        df = fl.find_avgCDR3LengthPerClone(patient_id)
        dfs.append(df)
    df = pd.concat(dfs)

    # creates figure and axis
    fig, ax = plt.subplots(figsize=(10,5))

    # draws boxplots on axis
    ax = sns.scatterplot(data=df,
                         x='no_seqs',
                         y='cdr3_aa_Length',
                         palette='mako')

    # sets figure attributes
    ax.set_xlabel('Clone Size (Number of Sequences)')
    ax.set_ylabel('Average CDR H3 Length')

    if save == True:

        # defines directory to store plots
        plotDir = os.getcwd() + os.sep + 'results' + os.sep + 'plots' + os.sep + \
                  'CDR_Lengths'
        if not os.path.exists(plotDir):
            os.makedirs(plotDir)

        # writes plot to png file in plotDir
        fig.savefig(plotDir + os.sep + 'avgCDRH3Len_vs_CloneSize.pdf')
