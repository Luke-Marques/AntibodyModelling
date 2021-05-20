import clean_repertoire_data as crd
import run_ABodyBuilder as rab


def main(patient_id_input):
    '''Takes patient_id_input variable as arg (either single patient_id or list
    of patient_ids, these should appear in name of repertoire data file in input
    directory) and generates a CSV file of paired heavy and light chain
    sequences for submission to ABodyBuilder. Then submits these to
    ABodyBuilder and downloads the results zip files to results directory.'''

    crd.save_submission_dataframes(patient_id_input)
    rab.run_all()
