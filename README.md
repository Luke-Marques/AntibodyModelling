## Instructions for use

The script pipeline.py incorporates the key functions from both clean_repertoire_data.py and run_ABodyBuilder.py, to determine paired sequences in the repertoire input file, save those paired sequences to a submission CSV file, and submit each paired sequence to ABodyBuilder for model generation.
These steps can be run by using the main() function in pipeline.py, where the patient_id_input variable must be a string (XX) corresponding to the name of the repertoire input files:

Pacbio_Oct19_lib_XX_withClones_collapsedByCell_annotated_withCellType.csv
XX_sequence_data.tsv

These files must be stored in an sub-directory within the working directory named 'inputs'.
This can be edited by editing the first 4 functions of the clean_repertoire_data.py script.

Running the pipeline.py main() function will call the clean_repertoire_data.py save_submission_dataframes() function and the run_ABodyBuilder.py run_all() function.

## Data Visualisation and Results Analysis Scripts

Also included in the repository is a variety of scripts which were used to generate the results and figures within my paper.