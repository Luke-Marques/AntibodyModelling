# Antibody Modelling and Structural Analysis Project
In this project, a novel pipeline was developed to transform antibody RNA codes, collected from individuals following vaccination trials, to 3-dimensional models which can then be used for structural analysis.

Within the codebase a numerous analyses of these structures and codes, including length analyses as well as a novel clustering method utilising alphabetic encoding, Levenshtein distance, and an unsupervised heirachical clustering approach.


## Disclaimer

This project was my first time implementing the programming and data skills I taught myself during the summer break between my 3rd and 4th years of university. Since starting, and completing this project, I have built upon these skills enormously, completing numerous online data science courses. A number of the functions in this codebase were written in a cumberson manner, often performing many tasks at once. Additionally, some of the methods used were slow (e.g. large Pandas dataframe joins) and could be improved upon.

Overall, a clearer structure to the project is also needed. This would help to improve the modularity of the project.

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