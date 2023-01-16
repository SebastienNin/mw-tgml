#!/usr/bin/env python3

# This script create appropriate libraries.csv file for scRNA-seq multiplexing using the cellplex kit.

import pandas
import os
import numpy as np
import hashlib
from pathlib import Path

if os.path.isfile("../mw-tgml/Sequencing_summary.xlsx"):
    # Read the Sequencing_summary.
    samples = pandas.read_excel("../mw-tgml/Sequencing_summary.xlsx", sheet_name="samples", engine="openpyxl")
    # Get samples to process, i.e those with type == cellplex
    samples_cellplex_1 = samples[samples['Type'].isin(['Cellplex'])]
    samples_cellplex_2 = samples_cellplex_1[samples_cellplex_1['Process'].isin(['yes'])]
    samples_cellplex_3 = samples_cellplex_2[samples_cellplex_2['Analysis_type'].isin(['Demultiplexage_Concatenation_Quantification_QC', 'Concatenation_Quantification_QC'])]
    samples_cellplex = samples_cellplex_3

    experiments_cellplex = samples_cellplex.Accession.unique()
    specie = samples_cellplex.Specie.unique()

    mwconf['config_targets'] = []

    for experiment in experiments_cellplex:
        # Variable to get the loop dataframe
        current_df = samples_cellplex[samples['Accession'] == experiment]
        kit_used = current_df.Kit_index.unique()[0]
        SPECIE = current_df.Specie.unique()[0]
        ACCESSION = current_df.Accession.unique()[0]

        info_to_write = "[gene-expression]\nreference," 

        # For now, on the platform we process only human and mouse.
        if SPECIE in ['human', 'Human', 'Homo_sapiens']:
            scrna_assembly = "GRCh38-2020-A"
        elif SPECIE in ['mouse', 'Mouse', 'Mus_musculus']:
            scrna_assembly = "mm10-2020-A"
        else:
            continue
       
        # Get path to reference
        info_to_write = info_to_write + os.getcwd() + "/out/tar/xvzf_genome_cellranger/wget/https/cf.10xgenomics.com/supp/cell-exp/refdata-gex-" + scrna_assembly + "\n"
        # Add expected cell number
        #expected_cell_nb_list = str(int(current_df.Expected_cell_number.unique()[0]))
        #print(expected_cell_nb_list)
        #info_to_write = info_to_write + "expect-cells," + str(int(current_df.Expected_cell_number.unique()[0])) + "\n\n" 
        info_to_write = info_to_write + "[libraries]\nfastq_id,fastqs,feature_types\n"
        
        for index, row in current_df.iterrows():
            info_to_write = info_to_write + row['Sample_Name'] + "," + os.getcwd() + "/out/cellranger/mkfastq" + row['Accession'] + "/" + "_".join(str(row['Run_Name']).split("_")[0:2]) + "/outs/fastq_path/" + str(row['Run_Name']) + "/" + str(row['Sample_ID']) + "," + row['Cellplex_feature_type'] + "\n"
            #info_to_write = info_to_write + row['Sample_Name'] + "," + os.getcwd() + "/out/cellranger/mkfastq" + row['Accession'] + "/" + "_".join(str(row['Run_Name']).split("_")[0:2]) + "/outs/fastq_path/" + "," + row['Cellplex_feature_type'] + "\n"

        info_to_write = info_to_write + "\n[samples]\nsample_id,cmo_ids,description\n" 

        for index, row in current_df.iterrows():
            if(row['Cellplex_feature_type'] == "Multiplexing Capture"):
                CMO_list = row['CMO_information'].split(";")
                for CMO in CMO_list:
                    CMO_split = CMO.split(",")
                    SAMPLE = CMO_split[0]
                    CMO_NB = CMO_split[1]
                    info_to_write = info_to_write + SAMPLE + "," + CMO_NB + "," + SAMPLE + "\n"

        #print(info_to_write)

        cellranger_multi_target = "out/cellranger/multi_" + scrna_assembly + "/cellranger/mkfastq" + ACCESSION + "/process_done"
        # Test to add ln for cellranger
        #cellranger_multi_prefix = "/cellranger/multi_" + scrna_assembly + "/cellranger/mkfastq" + ACCESSION
        #ln_cellranger_multi = "out/ln/alias/sst/by_run/run" + RUN + cellranger_multi_prefix

        filename = "out/csv/cellranger/multi/cellranger/mkfastq" + ACCESSION  + "/libraries.csv"
        tmp_filename = "out/csv/cellranger/multi/cellranger/mkfastq" + ACCESSION + "/libraries.csv.tmp"

        output_prefix = "out/csv/cellranger/multi/cellranger/mkfastq" + ACCESSION

        if current_df[['Process']].iloc[0, 0] == "yes":
            os.makedirs(output_prefix, exist_ok=True)
            myfile = open(filename, 'w')
            myfile.write(info_to_write)
            myfile.close()
            mwconf['config_targets'].append(cellranger_multi_target)
