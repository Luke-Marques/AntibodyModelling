import pandas as pd
import numpy as np
import itertools
import os
import glob
import clean_repertoire_data as crd
import matplotlib.pyplot as plt
import seaborn as sns
import read_pdb as rpdb
import create_sasa_plots as csp


%matplotlib inline
sns.set_theme()
sns.set_context('paper')
sns.set_palette('mako')
sns.set_style('white')

colours = sns.color_palette(palette = 'crest', n_colors = 6)


pd.options.display.html.table_schema = True
pd.options.display.max_rows = None
pd.options.display.max_columns = None


patient_id_list = ['YF189',
                   'YF191',
                   'YF192',
                   'YF199',
                   'YF200',
                   'YF203']
old_patient_id_list = ['YF189',
                       'YF191',
                       'YF192']
young_patient_id_list = ['YF199',
                         'YF200',
                         'YF203']


def write_dataframe_to_csv(df, dir_path, file_name):
    '''Takes pandas dataframe object as argument and writes to .csv file in
    directory'''

    file_path = dir_path + os.sep + file_name
    df.to_csv(file_path, index=False)


def get_aa_type(aa):
    '''Takes a single character string representing an amino acid as an arg
    and returns the classification of that amino acid'''

    if aa in ['R', 'H', 'K']:
        return 'Positively Charged'
    elif aa in ['D', 'E']:
        return 'Negatively Charged'
    elif aa in ['Q', 'N', 'S', 'T', 'Y', 'C']:
        return 'Uncharged'
    elif aa in ['A', 'I', 'L', 'M', 'F', 'V', 'P', 'G', 'W']:
        return 'Non-Polar'


def get_mut_dictionary(patient_id):
    '''Loops through the dataframe created using the 'create_sequence_dataframe'
    function and finds mutations between the germline and sequence amino acid
    code and returns them and their position as a nested dictoinary'''

    # runs 'create_sequence_dataframe' function
    df_sequence = crd.prepare_sequence_dataframe(patient_id)

    # filters 'df_sequence' based on it's 'productive' column
    df_sequence = df_sequence[df_sequence.productive == 'T']

    # removes rows where 'v_sequence_alignment_aa' values contains '*' character
    df_sequence = df_sequence[
        ~df_sequence.v_sequence_alignment_aa.str.contains('*', regex=False)]

    # runs 'create_repertoire_dataframe'
    df_repertoire = crd.prepare_repertoire_dataframe(patient_id) \
                    [['CloneID', 'CellID', 'SeqID', 'UseAsRef']]

    # filters dataframe for 'CloneID' values using 'get_top_heavy_clones'
    # funtion
    df_top_heavy_clones = df_repertoire[df_repertoire.CloneID \
                            .isin(crd.get_top_heavy_clones(patient_id))]

    # filters dataframe for 'CellID' values using 'CellID' values associated
    # with 'CloneID' values from 'df_top_heavy_clones'
    df_top_cells = df_repertoire[df_repertoire.CellID \
                    .isin(list(df_top_heavy_clones.CellID))]

    # list of all 'SeqID' values associated with top 'CloneID' values (H and L)
    top_SeqID_list = list(df_top_cells.SeqID)

    # filters 'df_repertoire' using 'top_SeqID_list'
    df_repertoire = df_repertoire[df_repertoire.SeqID.isin(top_SeqID_list)]

    # dataframe of top 'CloneID' values and 'SeqID' value filtered by 'UseAsRef'
    df_top_clones_ref = df_top_cells[df_top_cells.UseAsRef == True]
    df_top_clones_ref = df_top_clones_ref[['CloneID', 'SeqID']]

    # merges 'df_top_clone_ref' and 'df_sequence' to give dataframe of 'CloneID'
    # values associated with 'SeqID' values that are marked with
    # 'UseAsRef' == True
    df_top_clones_ref = df_top_clones_ref \
        .merge(df_sequence, left_on='SeqID', right_on='sequence_id')
    df_top_clones_ref = df_top_clones_ref.rename(columns={
        'v_germline_alignment_aa': 'v_germline_alignment_aa_ref'})
    df_top_clones_ref = df_top_clones_ref[[
        'CloneID',
        'v_germline_alignment_aa_ref']]

    # merges 'df_repertoire' and 'df_sequence'
    df = df_repertoire \
        .merge(df_sequence, left_on='SeqID', right_on='sequence_id')

    # merges 'df' and 'df_top_clones_ref'
    df = df.merge(df_top_clones_ref, left_on='CloneID', right_on='CloneID')
    df = df[[
        'CloneID',
        'SeqID',
        'v_germline_alignment_aa',
        'v_sequence_alignment_aa',
        'v_germline_alignment_aa_ref' ]]

    # loops through each 'SeqID' value in 'df_sequence' and finds differences
    # between 'v_germline_alignment_aa_ref' and 'sequene_alignment_aa' values
    # stores mutation position, germline aa and sequence aa in dictionary
    nested_dict = {}
    for row in df.itertuples():
        nested_dict[row.SeqID] = {}
        position = 1
        for aa in row.v_germline_alignment_aa_ref:
            if aa != 'X':
                germline_aa = aa
                if position <= len(row.v_sequence_alignment_aa):
                    sequence_aa = row.v_sequence_alignment_aa[position-1]
                    if germline_aa != sequence_aa:
                        nested_dict[row.SeqID][position] = {
                                'germline_aa': germline_aa,
                                'sequence_aa': sequence_aa,
                                'CloneID': row.CloneID}
            position += 1

    return nested_dict


def create_mut_dataframe(patient_id, save=False):
    '''Takes nested dictionary from 'get_mutation_dictionary' as arg and
    converts to multiindex dataframe'''

    mutation_dictionary = get_mut_dictionary(patient_id)

    # converts nested dictionary to dictionary of tuples for later
    dict_tuples = {(k1, k2):v2 for k1,v1 in mutation_dictionary.items() \
                               for k2,v2 in mutation_dictionary[k1].items()}

    # converts dictionary of tuples to multiindex dataframe
    df = pd.DataFrame([dict_tuples[i] for i in sorted(dict_tuples)],
        index=pd.MultiIndex.from_tuples(
            [i for i in sorted(dict_tuples.keys())]))

    # renames columns
    colnames = ['SeqID', 'MutPosition', 'germline_aa', 'sequence_aa', 'CloneID']
    df = df.reset_index()
    df.columns = colnames

    # adds 'Region' columns using conditions based on IMGT numbering scheme
    conditions = [
        (df.MutPosition > 0) & (df.MutPosition <= 26),
        (df.MutPosition > 26) & (df.MutPosition <= 38),
        (df.MutPosition > 38) & (df.MutPosition <= 55),
        (df.MutPosition > 55) & (df.MutPosition <= 65),
        (df.MutPosition > 65) & (df.MutPosition <= 104),
        (df.MutPosition > 104) & (df.MutPosition <= 117),
        (df.MutPosition > 117),]
    values = [
        'FR1',
        'CDR1',
        'FR2',
        'CDR2',
        'FR3',
        'CDR3',
        'FR4']
    df['Region'] = np.select(conditions, values)

    # adds two columns which define the type of amino acid
    # (polar charged, polar uncharged and non-polar)
    df['germline_aa_type'] = df.germline_aa.apply(lambda x: get_aa_type(x))
    df['sequence_aa_type'] = df.sequence_aa.apply(lambda x: get_aa_type(x))

    # creates repertoire dataframe from clean_repertoire_data.py module
    df_repertoire = crd.prepare_repertoire_dataframe(patient_id) \
        [['SeqID', 'Class', 'UseAsRef']]

    # merges df with df_repertoire on SeqID
    df = df.merge(df_repertoire, left_on='SeqID', right_on='SeqID')

    # # uncomment to convert dataframe back to MultiIndex DataFrame
    # df = df.set_index(['SeqID', 'MutPosition'])

    dir_path = os.getcwd() + os.sep + 'results' + os.sep + 'mutation_data'
    file_name = patient_id + '_mutation_data.csv'

    if save == True:
        write_dataframe_to_csv(df, dir_path, file_name)

    return df


def create_mut_dataframe_light(patient_id):
    '''Filters mut_dataframe to only include rows where 'SeqID' value is
    associated with 'Class' values corresonding to light chains'''

    df = create_mut_dataframe(patient_id)

    # filters df for light chain 'Class' values
    df = df[df.Class.isin(['K', 'L'])]

    return df


def create_mut_dataframe_heavy(patient_id):
    '''Filters mut_dataframe to only include rows where 'SeqID' value is
    associated with 'Class' values corresonding to heavy chains'''

    df = create_mut_dataframe(patient_id)

    # filters df for light chain 'Class' values
    df = df[df.Class.isin(['A', 'G', 'M'])]

    return df


def create_SASARef_dataframe(patient_id):
    '''Takes patient_id as arg and creates a dataframe of SeqID values for every
    sequence of which mutations are known using mutation_data.csv'''

    filepath = os.getcwd() + os.sep + 'results' + os.sep + 'mutation_data' \
               + os.sep + patient_id + '_mutation_data.csv'
    df = pd.read_csv(filepath)

    # get list of jobnames which have been modelled successfully via
    # ABodyBuilder script
    pdb_dir = os.getcwd() \
              + os.sep + 'results' \
              + os.sep + 'ABodyBuilder_results' \
              + os.sep + patient_id \
              + os.sep  + 'pdb_structures'
    pdb_file_paths = glob.glob(pdb_dir + os.sep + '*.pdb')
    jobname_list = []
    for pdb in pdb_file_paths:
        pdb_filename = pdb.split(os.sep)[-1]
        pdb_jobname = pdb_filename.split('.')[0]
        jobname_list.append(pdb_jobname)


    df = df[['CloneID', 'SeqID', 'Class']]
    df = df.groupby('SeqID', as_index=False).agg('first').reset_index(drop=True)
    df['jobname'] = df['SeqID'].str.split('_').str[0] + '_' + \
                    df['SeqID'].str.split('_').str[3]
    df['Modelled'] = False
    df['UseAsSASARef'] = False
    df.loc[df.jobname.isin(jobname_list), 'Modelled'] = True

    df_sequence = crd.prepare_sequence_dataframe(patient_id) \
                  [['sequence_id', 'v_identity']]

    df = df.merge(df_sequence, left_on='SeqID', right_on='sequence_id') \
           .drop(columns='sequence_id')

    df = df.sort_values(['CloneID', 'v_identity'], ascending=[True, False]) \
           .reset_index(drop=True)

    # get list of CloneID values in dataframe
    CloneID_list = df.CloneID.unique()

    # slice dataframe for each CloneID in list and check if slice containes any
    # rows with Modelled == True, if not, append CloneID to new list and
    # sets UseAsSASARef column to True for modelled jobname rows with highest
    # v_identity value clone
    df['jobname'] = \
        df.SeqID.str.split('_').str[0] + '_' + \
        df.SeqID.str.split('_').str[3]
    ignore_CloneID_list = []
    SASARef_SeqID_list = []
    counter = 1
    for CloneID in CloneID_list:
        df_slice = df[df.CloneID == CloneID]
        jobnames = df_slice.jobname.unique()
        singleChainPDBs = rpdb.get_singleChainPDBs(patient_id)
        check = any(item in jobnames for item in singleChainPDBs)
        if (df_slice.Modelled.any() == False) | (check == True):
            ignore_CloneID_list.append(CloneID)
        else:
            df_slice = df_slice[df_slice.Modelled == True]
            df_slice = df_slice.head(1).reset_index(drop=True)
            SeqID_value = df_slice.iloc[0].SeqID
            SASARef_SeqID_list.append(SeqID_value)
        print(str(counter) + '/' + str(len(CloneID_list)))
        counter += 1

    df.loc[df.SeqID.isin(SASARef_SeqID_list), 'UseAsSASARef'] = True

    # filters df to not include CloneIDs with no SASARef jobnames
    df = df[~df.CloneID.isin(ignore_CloneID_list)]

    return df


def create_mut_dataframe_with_SASA_data(patient_id, save=False):

    # gets mutation dataframe
    filepath = os.getcwd() + os.sep + 'results' + os.sep + 'mutation_data' \
               + os.sep + patient_id + '_mutation_data.csv'
    df_mut = pd.read_csv(filepath)

    # adds empty SASA columns to df_mut
    df_mut['Q(SASA)'] = np.nan
    df_mut['Exposed'] = ''

    # adds Chain column to df_mut
    df_mut['Chain'] = ''
    df_mut.loc[df_mut.Class.isin(['L', 'K']), 'Chain'] = 'L'
    df_mut.loc[df_mut.Class.isin(['A', 'G', 'M']), 'Chain'] = 'H'

    # gets SASARef dataframe
    df_SASARef = create_SASARef_dataframe(patient_id)

    # gets SASA dataframe
    filepath = os.getcwd() + os.sep + 'results' + os.sep + 'SASA_results' \
               + os.sep + patient_id + '_per_residue_SASA.csv'
    df_SASA = pd.read_csv(filepath)

    # loops through each unique CloneID value in df_SASARef
    CloneID_list = df_SASARef.CloneID.unique()

    # filters df_mut to only include clones from CloneID_list
    df_mut = df_mut[df_mut.CloneID.isin(CloneID_list)]

    for CloneID in CloneID_list:

        # slices dataframe on CloneID value
        df_SASARef_slice = df_SASARef[df_SASARef.CloneID == CloneID]

        # gets jobname of UseAsSASARef == True row
        df_SASARef_slice_jobname = df_SASARef_slice \
            [df_SASARef_slice.UseAsSASARef == True] \
            .reset_index()
        SASARef_jobname = df_SASARef_slice_jobname.iloc[0].jobname

        # gets SASA dataframe and slices with jobname == SASARef_jobname
        filepath = os.getcwd() + os.sep + 'results' + os.sep + 'SASA_results' \
                   + os.sep + patient_id + '_per_residue_SASA.csv'
        df_SASA = pd.read_csv(filepath)
        df_SASA_slice = df_SASA[df_SASA.jobname == SASARef_jobname]

        # slice mutation dataframe on CloneID
        df_mut_slice = df_mut[df_mut.CloneID == CloneID]

        # loops through heavy and light chain types
        for Chain in 'H', 'L':

            # slices SASA dataframe again to only include one chain
            df_SASA_slice_chain = df_SASA_slice[df_SASA_slice.Chain == Chain]

            # slices mutation dataframe again to only include one chain
            df_mut_slice_chain = df_mut_slice[df_mut_slice.Chain == Chain]

            # gets list of mutation positions from mutation dataframe
            position_list = df_mut_slice_chain.MutPosition.unique()

            # loops through list of positions
            if len(position_list) > 0:
                for position in position_list:
                    values = df_SASA_slice_chain.Position_Insertion.values
                    if str(position) in values:

                        # finds SASA value of each position
                        df_SASA_slice_position = \
                            df_SASA_slice_chain \
                            [df_SASA_slice_chain.Position_Insertion \
                                == str(position)].reset_index(drop=True)
                        position_SASA = \
                            df_SASA_slice_position.iloc[0]['Q(SASA)']

                        df_mut.loc[(df_mut.CloneID == CloneID) & \
                                   (df_mut.Chain == Chain) & \
                                   (df_mut.MutPosition == position),
                                   'Q(SASA)'] = position_SASA

    # set Exposed column to True/False depending on Q(SASA) column
    df_mut.loc[(df_mut['Q(SASA)'] != np.nan) & \
               (df_mut['Q(SASA)'] >= 0.15),
               'Exposed'] = True
    df_mut.loc[(df_mut['Q(SASA)'] != np.nan) & \
               (df_mut['Q(SASA)'] < 0.15),
               'Exposed'] = False

    # saves dataframe to .csv file
    if save == True:
        dir_path = os.getcwd() + os.sep + 'results' + os.sep + 'mutation_data'
        file_name = patient_id + '_mutation_data_with_SASA.csv'
        write_dataframe_to_csv(df_mut, dir_path, file_name)

    return df_mut


def find_numMutsPerSASA(df):
    '''Takes mutation dataframe from 'create_mut_dataframe_with_SASA_data'
    function and returns a nested dictionary of number of mutations at exposed
    and buried positions in each chain, H and L.'''

    # creates empty dictionary
    d = {}

    # loops through Chain types
    for Chain in 'H', 'L':

        # slices dataframe by Chain
        df = df[df.Chain == Chain]

        # count all mutations at exposed, buried or unknown positions
        d[Chain] = {'ExposedCount': df.Exposed.value_count().loc[True],
                    'BuriedCount': df.Exposed.value_count().loc[False],
                    'UnknownCount': df.Exposed.isna().sum()}

    return d


def find_numMutsPerRegionPerSASA(df):
    '''Takes mutation dataframe from 'create_mut_dataframe_with_SASA_data'
    function and returns a nested dictionary of number of mutations at exposed,
    buried and unknown positions in each chain and at each region.'''

    # defines list of Region names to loop through
    RegionList = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']

    # creates empty dictionary
    d = {}

    # loops through Chain types
    for Chain in 'H', 'L':

        # slices dataframe by Chain
        df_Chain = df[df.Chain == Chain]

        # creates key in d of Chain
        d[Chain] = {}

        # loops through Regions
        for Region in RegionList:

            # slices dataframe by Region
            df_Region = df_Chain[df_Chain.Region == Region]

            # counts number of exposed, buried and unknown values for the
            # column 'Exposed' in dataframe
            numNaN = df_Region.Exposed.isna().sum()
            numTrue = df_Region.Exposed.sum()
            numFalse = len(df_Region) - numNaN - numTrue

            # count all mutations at exposed, buried or unknown positions
            d[Chain][Region] = {'ExposedCount': numTrue,
                                'BuriedCount': numFalse,
                                'UnknownCount': numNaN}

    return d


def get_avgNumExposedPerRegion(patient_id, chain):

    # get SASA dataframe
    filepath = os.getcwd() + os.sep + 'results' + os.sep + 'SASA_results' + \
               os.sep + 'all_patients_SASA_meta_' + chain + '.csv'
    df = pd.read_csv(filepath)

    # slice dataframe on patient_id
    df = df[df.patient_id == patient_id]

    # adds column 'Position' which does not use insertion codes
    Positions = df.Position_Insertion.values.tolist()
    Positions = [str(x) for x in Positions]
    Positions = [x[:3] for x in Positions]
    df['Position'] = Positions
    df.Position = df.Position.astype(int)

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

    return d


def create_MutSASA_grouped_barplots(patient_id_list):

    # defines Region names and stores in list
    RegionList = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']

    # stores lengths of each region as variables
    FR1_len = 26-0
    FR2_len = 55-38
    FR3_len = 104-65
    CDR1_len = 38-26
    CDR2_len = 65-55
    CDR3_len = 117-104

    # creates empty lists to store dataframes
    df_list_all = []
    df_list_old = []
    df_list_young = []

    d_SASA_list = []

    for patient_id in patient_id_list:

        # gets mutation dataframe with SASA data for patient
        filepath = os.getcwd() + os.sep + 'results' + os.sep + 'mutation_data' \
                   + os.sep + patient_id + '_mutation_data_with_SASA.csv'
        df = pd.read_csv(filepath)

        # generate mutation SASA count dictionary
        d = find_numMutsPerRegionPerSASA(df)

        # loops through Chain types
        for Chain in 'H', 'L':

            # defines chain_name var for use in naming files and titles
            if Chain == 'H':
                chain_name = 'Heavy'
            elif Chain == 'L':
                chain_name = 'Light'

            # gets avg SASA count dictionary for patient
            d_SASA = get_avgNumExposedPerRegion(patient_id, Chain)
            d_SASA_list.append(d_SASA)

            # converts dictionary to dataframe
            df_Chain = pd.DataFrame.from_dict(d[Chain])

            # reshapes dataframe
            df_Chain = df_Chain.reset_index(drop=False) \
                               .rename(columns={'index': 'SASA_cat'})
            df_Chain = pd.melt(df_Chain,
                               id_vars=['SASA_cat'],
                               value_vars=RegionList,
                               value_name='value')
            df_Chain = df_Chain.rename(columns={'variable': 'Region',
                                                'value': 'Count'})
            df_Chain = df_Chain.replace({'ExposedCount': 'Exposed',
                                         'BuriedCount': 'Buried',
                                         'UnknownCount': 'Unknown'})

            # adds columns for patient_id, Chain type, proportional count, and
            # percentage of proportional count
            df_Chain['patient_id'] = patient_id
            df_Chain['Chain'] = Chain
            df_Chain['Prop_Count'] = np.nan
            df_Chain['Prop_pct'] = np.nan

            # creates list of Regions
            Regions = df_Chain.Region.unique()

            # creates column to store avgNumExposed per Region
            df_Chain['avgNumExposed'] = 0
            df_Chain['avgNumBuried'] = 0
            df_Chain['avgNumUnknown'] = 0
            for Region in Regions:
                df_Chain.loc[df_Chain.Region == Region, 'avgNumExposed'] = \
                    d_SASA[Region]['Exposed']
                df_Chain.loc[df_Chain.Region == Region, 'avgNumBuried'] = \
                    d_SASA[Region]['Buried']
                df_Chain.loc[df_Chain.Region == Region, 'avgNumUnknown'] = \
                    d_SASA[Region]['Unknown']

            # sets Prop_Count column to scale count to length of region
            df_Chain.loc[df_Chain.Region == 'FR1', 'Prop_Count'] = \
                df_Chain.Count / FR1_len
            df_Chain.loc[df_Chain.Region == 'FR2', 'Prop_Count'] = \
                df_Chain.Count / FR2_len
            df_Chain.loc[df_Chain.Region == 'FR3', 'Prop_Count'] = \
                df_Chain.Count / FR3_len
            df_Chain.loc[df_Chain.Region == 'CDR1', 'Prop_Count'] = \
                df_Chain.Count / CDR1_len
            df_Chain.loc[df_Chain.Region == 'CDR2', 'Prop_Count'] = \
                df_Chain.Count / CDR2_len
            df_Chain.loc[df_Chain.Region == 'CDR3', 'Prop_Count'] = \
                df_Chain.Count / CDR3_len

            # # scale propcount to avg number of SASA category in region
            # for Region in Regions:
            #     df_Chain.loc[(df_Chain.SASA_cat == 'Exposed') & \
            #                  (df_Chain.Region == Region), 'Prop_Count'] = \
            #         df_Chain.Prop_Count / d_SASA[Region]['Exposed']
            #     df_Chain.loc[(df_Chain.SASA_cat == 'Buried') & \
            #                  (df_Chain.Region == Region), 'Prop_Count'] = \
            #         df_Chain.Prop_Count / d_SASA[Region]['Buried']
            #     df_Chain.loc[(df_Chain.SASA_cat == 'Unknown') & \
            #                  (df_Chain.Region == Region) & \
            #                  (df_Chain.avgNumUnknown != 0), 'Prop_Count'] = \
            #         df_Chain.Prop_Count / d_SASA[Region]['Unknown']

            df_Chain = df_Chain[df_Chain.SASA_cat != 'Unknown']

            # gets total for Prop_Count column
            Prop_Counts_total = df_Chain.Prop_Count.sum()

            # calculates counts as percentage of total counts
            df_Chain.loc[df_Chain.Region == 'FR1', 'Prop_pct'] = \
                df_Chain.Prop_Count / Prop_Counts_total
            df_Chain.loc[df_Chain.Region == 'FR2', 'Prop_pct'] = \
                df_Chain.Prop_Count / Prop_Counts_total
            df_Chain.loc[df_Chain.Region == 'FR3', 'Prop_pct'] = \
                df_Chain.Prop_Count / Prop_Counts_total
            df_Chain.loc[df_Chain.Region == 'CDR1', 'Prop_pct'] = \
                df_Chain.Prop_Count / Prop_Counts_total
            df_Chain.loc[df_Chain.Region == 'CDR2', 'Prop_pct'] = \
                df_Chain.Prop_Count / Prop_Counts_total
            df_Chain.loc[df_Chain.Region == 'CDR3', 'Prop_pct'] = \
                df_Chain.Prop_Count / Prop_Counts_total

            # appends dataframe to df_list
            df_list_all.append(df_Chain)

            # appends dataframe to either df_list_old or df_list_young
            if patient_id in old_patient_id_list:
                df_list_old.append(df_Chain)
            if patient_id in young_patient_id_list:
                df_list_young.append(df_Chain)

            # creates figure for plot
            plt.figure(figsize=(10,8))

            # draws barplot on figure
            sns.barplot(x = 'Region',
                        y = 'Prop_pct',
                        hue = 'SASA_cat',
                        data=df_Chain,
                        palette='mako')

            # sets figure attributes
            ylabel = 'Density of Mutations Relative to Length of Region'
            plt.ylabel(ylabel, size=14)
            xlabel = 'Region'
            plt.xlabel(xlabel, size=14)
            title = 'Density of Mutations within Regions, Grouped by SASA' + \
                    ' Category for Patient ' + patient_id + ', ' + chain_name + \
                    ' Chain'
            plt.title(title, size=18)

            # creates directory for plot if it doesn't already exist
            plot_path = os.getcwd() + os.sep + 'results' + os.sep + \
                        'mutation_data' + os.sep + 'MutGroupedBySASA_Barplots' \
                        + os.sep + chain_name
            if not os.path.exists(plot_path):
                os.makedirs(plot_path)

            # saves plot to directory
            file_name = patient_id + '_' + chain_name + \
                        'Chain_MutsGroupedBySASA_Barplot.svg'
            plt.savefig(plot_path + os.sep + file_name, bbox_inches='tight')

    # creates counter variable to select correct name for plots
    counter = 0
    group_names = ['AllPatients', 'OldPatients', 'YoungPatients']

    for df_list in [df_list_all, df_list_old, df_list_young]:

        # selects correct patient group name
        group_name = group_names[counter]

        # concatenates df_list to single dataframe
        df_all = pd.concat(df_list).reset_index(drop=True)

        # loops through Chain types again
        for Chain in 'H', 'L':

            # defines chain_name var for use in naming files and titles
            if Chain == 'H':
                chain_name = 'Heavy'
            elif Chain == 'L':
                chain_name = 'Light'

            # slices dataframe on Chain
            df_all_Chain = df_all[df_all.Chain == Chain]

            # creates figure for plot
            plt.figure(figsize=(10,5))

            # draws barplot on figure
            sns.barplot(x = 'Region',
                        y = 'Prop_pct',
                        hue = 'SASA_cat',
                        data=df_all_Chain,
                        palette='mako')

            # sets figure attributes
            # ylabel = 'Density of Mutations Relative to Length of Region'
            # plt.ylabel(ylabel, size=14)
            # xlabel = 'Region'
            # plt.xlabel(xlabel, size=14)
            # title = 'Density of Mutations within Regions, Grouped by SASA' + \
            #         ' Category for ' + group_name + ', ' + chain_name + \
            #         ' Chain'
            # plt.title(title, size=18)

            # creates directory for plot if it doesn't already exist
            plot_path = os.getcwd() + os.sep + 'results' + os.sep + \
                        'mutation_data' + os.sep + 'MutGroupedBySASA_Barplots' \
                        + os.sep + group_name
            if not os.path.exists(plot_path):
                os.makedirs(plot_path)

            # saves plot to directory
            file_name = group_name + '_' + chain_name + \
                        'Chain_MutsGroupedBySASA_Barplot.svg'
            plt.savefig(plot_path + os.sep + file_name, bbox_inches='tight')

        # increases counter
        counter += 1


create_MutSASA_grouped_barplots(patient_id_list)
