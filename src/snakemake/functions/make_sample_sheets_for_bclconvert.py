#!/usr/bin/env python3

# This script produce appropriate SampleSheets for all samples requiring it from Sequencing_summary.xlsx
#if not os.path.isfile("../mw-sst/Sequencing_summary.csv"):

# Modified 03/09/2020
# if the data are single-cell RNA-seq bcl, we must not add 'I5_Index_ID' and 'index2' in the samplesheet or cellranger mkfastq throws an error.

# Modified 09/10/2020
# Merged scrna_bcl and bcl in the sequencing_summary after discussion with biologists.

import pandas
import os
import numpy as np
import hashlib
import re
import subprocess
from pathlib import Path

if os.path.isfile("../mw-tgml/Sequencing_summary.xlsx"):
    samples = pandas.read_excel("../mw-tgml/Sequencing_summary.xlsx", sheet_name="samples", engine="openpyxl")
    samples_process_yes = samples[samples['Process'] == "yes"]
    samples_from_bcl = samples_process_yes[samples_process_yes['Origin'].isin(['bcl_NS2000_p1p2', 'bcl_NS2000_p3', 'bcl_NS2000'])]
    # Add adapter file reading and dictionnary creation to add adapter sequences in the samplesheet.
    if os.path.isfile("../mw-tgml/src/snakemake/tables/adapter.tsv"):
        adapter_dict = {}
        with open("../mw-tgml/src/snakemake/tables/adapter.tsv", "r") as adapter_file:
            # Skip first line containing header
            next(adapter_file)
            # Fill the dictionnary with all kits information:   Kit Adapter1    Adapter2
            for line in adapter_file:
                splitted_line = line.split("\t")
                adapter_dict[splitted_line[0]] = [splitted_line[1], splitted_line[2].replace("\n", "")]

    experiments_from_bcl = samples_from_bcl.Accession.unique()

    mwconf['bclconvert_targets'] = []

    for experiment in experiments_from_bcl:
        # Variable to get the loop dataframe
        current_df = samples_from_bcl.loc[samples['Accession'] == experiment]
        current_df.rename({'Run_Name' : 'Sample_Project'}, axis=1, inplace=True)
        kit_used = current_df.Kit_index.unique()[0]

        # Condition to determine wether data are single-cell and single or double indexed
        if(['scRNA-seq', 'scRNA_HTO', 'Cellplex'] in current_df.Type.unique()):
            bclconvert_prefix = "out/cellranger/mkfastq/" + experiment
            bclconvert_target = bclconvert_prefix + "/Reports/html/tree.html"
            # 8/10/2021 Add fillna to correct an error where I5_Index_ID = NaN but is not detected by str(row['I5_Index_ID']) != "NaN"
            # fillna replace all na by ""
            #data_tmp = current_df[['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_10X','I7_Index_ID','Index', 'I5_Index_ID', 'Index2', 'Run_Name','Description']].fillna("")
            data_tmp = current_df[['sample_name','Sample_ID','Sample_Name','Index_10X','I7_Index_ID','Index', 'I5_Index_ID', 'Index2', 'Sample_Project','Description']].fillna("")
            # 2021 04 07 Update to enable dual indexing demultiplexing. Cellranger can read index name and will replace it by the index sequence. It also detect dual index or single index.
            # 2021-04-23 Update the conditions to detect if samples have double index or not.
            for index, row in data_tmp.iterrows():
                if row['Index_10X'] == "yes":
                    mRNA_index=index
                    index_name=row["I7_Index_ID"]
                    data_tmp.at[mRNA_index, "Index"] = index_name
                    if(str(row['I5_Index_ID']) != ""):
                        data_tmp.at[mRNA_index, "I5_Index_ID"] = index_name
                        data_tmp.at[mRNA_index, "Index2"] = index_name
                        data_dropped = data_tmp.drop(columns=["Index_10X"])
                    else:
                        data_dropped = data_tmp.drop(columns=["Index_10X", "I5_Index_ID", "Index2"])
                elif(str(row['I5_Index_ID']) == ''):
                    data_dropped = data_tmp.drop(columns=["Index_10X", "I5_Index_ID", "Index2"])
                else:
                    data_dropped = data_tmp.drop(columns=["Index_10X"])
            data = data_dropped
        #elif(['bcl_no_mismatch'] in current_df.Origin.unique()):
            #bcl2fastq_prefix = "out/bcl2fastq/_--no-lane-splitting_--barcode-mismatches_0/" + experiment
            #bcl2fastq_target = bcl2fastq_prefix + "/Reports/html/tree.html"
            #data = current_df[['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','I7_Index_ID','Index','I5_Index_ID','Index2','Run_Name','Description']]
            #data = current_df[['Sample_ID','Sample_Name','I7_Index_ID','Index','I5_Index_ID','Index2','Run_Name','Description']]
            #data = current_df[['Sample_ID','Sample_Name','I7_Index_ID','Index','I5_Index_ID','Index2','Sample_Project','Description']]
        else:
            # By default, create the fastq for index reads! Important to use when debugging!
            bclconvert_prefix = "out/bcl-convert/_--force/" + experiment
            bclconvert_target = bclconvert_prefix + "/Logs/FastqComplete.txt"
            #data = current_df[['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','I7_Index_ID','Index','I5_Index_ID','Index2','Run_Name','Description']]
            #data = current_df[['Sample_ID','Sample_Name','I7_Index_ID','Index','I5_Index_ID','Index2','Run_Name','Description']]
            data = current_df[['sample_name','Sample_ID','Sample_Name','I7_Index_ID','Index','I5_Index_ID','Index2', 'Sample_Project','Description']]
        # Sort by index to prevent trouble with mix of ID and numbers
        data_to_write = data[['sample_name','Index','Index2']]
        data_to_write = data_to_write.rename(columns={'sample_name':'Sample_ID', 'Index':'Index', 'Index2':'Index2'})
        data_to_write = data_to_write.sort_index()
        filename = bclconvert_prefix + "/SampleSheet.csv"
        tmp_filename = bclconvert_prefix + "/SampleSheet.csv.tmp"

        filename = filename.replace('//', '/')
        tmp_filename = tmp_filename.replace('//', '/')

        # Condition to detect SE or PE and create adapter string to add to the file
        if len(str(current_df[['Reads']].iloc[0,0]).split("x")) == 1:
            readCycles = 'Read1Cycles,' + str(current_df[['Reads']].iloc[0,0]).split("x") + "\n"
            adapterString = 'AdapterRead1,' + adapter_dict[kit_used][0] + '\n'
        elif len(str(current_df[['Reads']].iloc[0,0]).split("x")) == 0:
            eprint("Sample " + current_df[['Sample_Name']] + " does not have read size in column 'Reads', please, check the Sequencing_Sumarry")
        else:
            readCycles = 'Read1Cycles,' + str(current_df[['Reads']].iloc[0,0]).split("x")[0] + "\n"
            readCycles = readCycles + 'Read2Cycles,' + str(current_df[['Reads']].iloc[0,0]).split("x")[1] + "\n"
            adapterString = 'AdapterRead1,' + adapter_dict[kit_used][0] + '\n' + 'AdapterRead2,' + adapter_dict[kit_used][1] + '\n\n'

        # Condition to detect single or double index
        if len(str(current_df[['Index']].iloc[0,0])) != 0:
            index1length = 'Index1Cycles,' + str(len(str(current_df[['Index']].iloc[0,0]))) + "\n"
        if len(str(current_df[['Index2']].iloc[0,0])) != 0:
            index2length = 'Index2Cycles,' + str(len(str(current_df[['Index2']].iloc[0,0]))) + "\n"

        # Get the bcl-convert version
        # Exécuter la commande bcl-convert -V et récupérer la sortie
        output = subprocess.check_output(['bcl-convert', '-V'], stderr=subprocess.STDOUT)

        # Convertir la sortie en une chaîne de caractères
        output_str = output.decode('utf-8').strip()

        cmd_split = output_str.split()
        version_full = cmd_split[2]
        version_split = version_full.split(".")
        # Finally, get the version
        version = ".".join(version_split[3:6])

        ## SampleSheet should be created only if it does not exist
        # to ensure bcl2fastq is not redone only because the SampleSheet is touched
        # by the process.
        # Maybe replace this with a hash check to still write a new SampleSheet if the
        # infos from Sequencing_summary are different.
        #if not os.path.isfile(filename):
        os.makedirs(bclconvert_prefix, exist_ok = True)
        header=('[Header]\n'
                'FileFormatVersion,2\n'
                'RunName,' + current_df[['Sample_Project']].iloc[0,0]+'\n'
                'InstrumentPlatform,NextSeq1k2k\n'
                '[Reads]\n' +
                readCycles +
                index1length +
                index2length + "\n" + 
                '[BCLConvert_Settings]\n'
                'SoftwareVersion,' + version + '\n' +  
                adapterString + 
                '[BCLConvert_Data]\n')

        if current_df[['Process']].iloc[0,0] == "yes":
            myfile = open(filename, 'w')
            myfile.write(header)
            myfile.write(data_to_write.to_csv(index=False))
            myfile.close()

        #if(['scRNA-seq', 'snRNA-seq', 'Cellplex'] not in current_df.Type.unique()):
        #    adapter_content=('[Settings]\n'
        # 'Adapter,' + adapter_dict[kit_used][0] + '\n'
        # 'AdapterRead2,' + adapter_dict[kit_used][1] + '\n\n'
        # '[Data]\n')
        #else:
        #    adapter_content=("\n"
        #    '[Data]\n')

        # To avoid id, double slashes "//" are replaced by single one.
        bclconvert_target = bclconvert_target.replace('//', '/')
        # Add a condition to read the process column of the summary.
        # Sometimes, we may not want to run bcl2fastq.
        if current_df[['Process']].iloc[0,0] == "yes":
            mwconf['bclconvert_targets'].append(bclconvert_target)
