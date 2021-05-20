import os
import sys
import urllib.request, urllib.error, urllib.parse
import mechanize
from time import sleep
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
import wget
from Bio.SeqIO.FastaIO import SimpleFastaParser
import glob
import shutil
from zipfile import ZipFile

# url for impact on stability (uncomment if stability calculation is needed)
# url = "http://bleoberis.bioc.cam.ac.uk/mcsm/stability"
url = "http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/abodybuilder/"

user_agent = 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)'
headers = {'User-Agent': user_agent}


def generate_df(file):
    '''Takes either a .fasta, .txt or .csv file containing the columns
    'heavychain', 'lightchain', 'jobname' and 'blacklist' and converts the file
    to a dataframe for use in future fucntions'''

    print('Converting ' + file + ' to dataframe...', end=' ')
    if file.endswith('.fasta'):
        # Parses fasta file using biopython module
        with open(file) as fasta_file:
            seq_name = []
            chain_type = []
            seq = []
            for title, sequence in SimpleFastaParser(fasta_file):
                seq_name.append(title.split("|")[0])
                chain_type.append(title.split("|")[1])
                seq.append(sequence)
        # Converts to dataframe
        data = {'seq_name': seq_name, 'chain_type': chain_type, 'seq': seq}
        df = pd.DataFrame(data, columns=['seq_name', 'chain_type', 'seq'])
        df = df.pivot(index='seq_name', columns=['chain_type'], values='seq')
        df = df.rename(columns={'H': 'heavy_chain', 'L': 'light_chain'}).reset_index().rename_axis(None, axis=1)
        df['blacklist'] = np.unique(seq_name)
        print('Done!')
        return df
    elif file.endswith('.csv') or file.endswith('.txt'):
        # Converts file to dataframe
        df = pd.read_csv(file)
        print('Done!')
        return df
    else:
        print('Please enter a valid fasta, txt, or csv file.')


def run_abodybuilder(heavy_chain, light_chain, jobname, blacklist):
    '''Given three strings: one heavy chain AA sequence, one light chain AA
    sequence, and one jobname to act as an id, submit to ABodyBuilder and
    downloads a zip file to the cwd/results/ABodyBuilder_results (will create
    dirs if they do not exist). Returns a list of failed jobs' results urls.
    Optionally take an extra argument, blacklist, (list of pdb files for
    ABodyBuilder to explicitly omit when searching for templates)'''

    # creates result directory
    cwd = os.getcwd()
    results_path = cwd + os.sep + 'results' + os.sep + 'ABodyBuilder_results'
    zip_path = results_path + os.sep + 'zip_files'
    if not os.path.exists(zip_path):
        os.makedirs(zip_path)
    # opens ABodyBuilder url and fills form for job submission
    br = mechanize.Browser(url)
    br.open(url)
    form = br.forms()[0]
    form['hchain'] = heavy_chain
    form['lchain'] = light_chain
    form['jobname'] = jobname
    if blacklist is not None:
        form['blacklist'] = blacklist
    # submits job
    request = form.click()
    # gets results url
    br2 = mechanize.Browser(url)
    br2.open(request)
    response = br2.response().geturl()
    # waits 15 seconds for job to complete
    sleep(15)
    fail_url = None
    # while loop checks if job has completed
    while True:
        htmlout = urllib.request.urlopen(response).read()
        soup = BeautifulSoup(htmlout, 'lxml')
        # if it has failed: saves failed results url to fail_url
        if 'Something went wrong!' in htmlout.decode('utf-8'):
            fail_url = response
            errorMsg = 'Job failed. Jobname and url will be saved to ' + \
                       jobname[:5] + '_failed_jobs.csv in cwd.'
            print(errorMsg)
            break
        # if it has not completed: waits 15 seconds and retries
        elif soup.find('div', id='results') is None:
            sleep(15)
        # if it has completed: downloads results zip file to results directory
        elif soup.find('div', id='results') is not None:
            htmlout = urllib.request.urlopen(response).read()
            soup = BeautifulSoup(htmlout, 'lxml')
            table = soup.find('div', id='results').find('table')
            zip_entry = table.findAll('tr')[-1]
            if 'Zip file containing all output from modelling job' in \
                    str(zip_entry.findAll('td')[0]):
                link = zip_entry.findAll('a')[0]['href']
                download_url = 'http://opig.stats.ox.ac.uk' + link
                destination = zip_path + os.sep + jobname + '_results.zip'
                print('Saving results file to ' + zip_path, end=' ')
                wget.download(download_url, destination)
                print('Done!')
            break
    return fail_url


def get_results(file_path):
    '''Given .csv (or similar) file with columns 'jobname', 'heavy_chain',
    'light_chain', and 'blacklist', creates dataframe from file, using
    create_dataframe function, loops through rows and submits a job to
    ABodyBuilder, using run_ABodyBuilder function, per row unless same job has
    been submitted before'''

    cwd = os.getcwd()
    df = generate_df(file_path)
    fail_dict = {'jobname': [], 'fail_url': []}
    # loops through each row in dataframe (could change to df.itertuples to improve performance)
    for index, row in df.iterrows():
        # prevents running same job twice
        zip_path = cwd + os.sep + 'results' + os.sep + 'ABodyBuilder_results' + os.sep + 'zip_files' + os.sep + row.jobname + '_results.zip'
        if not os.path.exists(zip_path):
            if row.blacklist:
                print('Running ABodyBuilder for ' + row.jobname)
                out = run_abodybuilder(row.heavy_chain, row.light_chain, row.jobname, None)
                if out is not None:
                    fail_dict['jobname'].append(row.jobname)
                    fail_dict['fail_url'].append(out)
            else:
                print('Running ABodyBuilder for ' + row.jobname)
                out = run_abodybuilder(row.heavy_chain, row.light_chain, row.jobname, None)
                if out is not None:
                    fail_dict['jobname'].append(row.jobname)
                    fail_dict['fail_url'].append(out)
            # prevents overloading ABodyBuilder server
            sleep(45)
    # saves dataframe of failed jobs' jobnames and failure urls to check failure reason
    fail_df = pd.DataFrame.from_dict(fail_dict)
    fail_name = file_path[-25:-20] + '_failed_jobs.csv'
    fail_path = os.getcwd() + os.sep + 'results' + os.sep + 'ABodyBuilder_results' + os.sep + fail_name
    fail_df.to_csv(fail_path, index=False)


def extract_pdb_files():
    '''Extracts zip files to patient directory in results directory and then
    copies the pdb structure file from each zip file to pdb_structure directory
    in patient directory'''

    cwd = os.getcwd()
    results_path = cwd + os.sep + 'results' + os.sep + 'ABodyBuilder_results'
    zip_files = glob.glob(cwd + os.sep + 'results' + os.sep + 'ABodyBuilder_results' + os.sep + 'zip_files' + os.sep + '*.zip')
    for zip_path in zip_files:
        zip_file = zip_path.rsplit(os.sep, 1)[-1]
        jobname = zip_file.rsplit('_')[0] + '_' + zip_file.rsplit('_')[1]
        patient_id = jobname[:5]
        if not os.path.exists(results_path + os.sep + patient_id + os.sep + 'all_results'):
            os.makedirs(results_path + os.sep + patient_id + os.sep + 'all_results')
        if not os.path.exists(results_path + os.sep + patient_id + os.sep + 'pdb_structures'):
            os.makedirs(results_path + os.sep + patient_id + os.sep + 'pdb_structures')
        with ZipFile(zip_path, 'r') as zipObj:
            zipObj.extractall(
                results_path + os.sep + patient_id + os.sep + 'all_results' + os.sep + jobname)
    for subdir, dirs, files in os.walk(results_path):
        for file_name in files:
            file_path = subdir + os.sep + file_name
            if file_path.endswith('rank1_imgt_scheme.pdb'):
                patient_id = subdir.rsplit(os.sep)[-5]
                new_file_name = subdir.rsplit(os.sep)[-3] + '.pdb'
                pdb_dir = results_path + os.sep + patient_id + os.sep + 'pdb_structures'
                if not os.path.exists(pdb_dir):
                    os.makedirs(pdb_dir)
                pdb_destination = pdb_dir + os.sep + new_file_name
                shutil.copy(file_path, pdb_destination)


def run_all():
    submission_path = os.getcwd() + os.sep + 'results' + os.sep + 'submission_files'
    for root, dirs, files in os.walk(submission_path):
        for file in files:
            file_path = submission_path + os.sep + file
            get_results(file_path)
    extract_pdb_files()


if __name__ == "__main__":
    if len(sys.argv) == 2:
        get_results(sys.argv[1])
    else:
        print("usage: python run_ABodyBuilder.py input_file ")
