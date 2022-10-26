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
from pathlib import Path

if os.path.isfile("../mw-tgml/Sequencing_summary.xlsx"):
    samples = pandas.read_excel("../mw-tgml/Sequencing_summary.xlsx", sheet_name="samples", engine="openpyxl")
    samples_from_bcl = samples[samples['Origin'].isin(['bcl', 'bcl_no_mismatch', 'bcl_index_generation'])]
    # Add adapter file reading and dictionnary creation to add adapter sequences in the samplesheet.
    if os.path.isfile("../mw-lib/src/snakemake/tables/adapter.tsv"):
        adapter_dict = {}
        with open("../mw-lib/src/snakemake/tables/adapter.tsv", "r") as adapter_file:
            # Skip first line containing header
            next(adapter_file)
            # Fill the dictionnary with all kits information:   Kit Adapter1    Adapter2
            for line in adapter_file:
                splitted_line = line.split("\t")
                adapter_dict[splitted_line[0]] = [splitted_line[1], splitted_line[2].replace("\n", "")]

    experiments_from_bcl = samples_from_bcl.Accession.unique()

    mwconf['bcl2fastq_targets'] = []

    for experiment in experiments_from_bcl:
        # Variable to get the loop dataframe
        current_df = samples_from_bcl[samples['Accession'] == experiment]
        kit_used = current_df.Kit_index.unique()[0]

        # Condition to determine wether data are single or double indexed
        if(['scRNA-seq', 'scRNA_HTO', 'Cellplex'] in current_df.Type.unique()):
            bcl2fastq_prefix = "out/cellranger/mkfastq/" + experiment
            bcl2fastq_target = bcl2fastq_prefix + "/Reports/html/tree.html"
            # 8/10/2021 Add fillna to correct an error where I5_Index_ID = NaN but is not detected by str(row['I5_Index_ID']) != "NaN"
            # fillna replace all na by ""
            data_tmp = current_df[['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','Index_10X','I7_Index_ID','Index', 'I5_Index_ID', 'Index2', 'Sample_Project','Description']].fillna("")
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
        elif(['bcl_no_mismatch'] in current_df.Origin.unique()):
            bcl2fastq_prefix = "out/bcl2fastq/_--no-lane-splitting_--barcode-mismatches_0/" + experiment
            bcl2fastq_target = bcl2fastq_prefix + "/Reports/html/tree.html"
            data = current_df[['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','I7_Index_ID','Index','I5_Index_ID','Index2','Sample_Project','Description']]
        else:
            # By default, create the fastq for index reads! Important to use when debugging!
            bcl2fastq_prefix = "out/bcl2fastq/_--no-lane-splitting_--create-fastq-for-index-reads/" + experiment
            bcl2fastq_target = bcl2fastq_prefix + "/Reports/html/tree.html"
            data = current_df[['Sample_ID','Sample_Name','Sample_Plate','Sample_Well','I7_Index_ID','Index','I5_Index_ID','Index2','Sample_Project','Description']]
        
        # Sort by index to prevent trouble with mix of ID and numbers
        data = data.sort_index()
        filename = bcl2fastq_prefix + "/SampleSheet.csv"
        tmp_filename = bcl2fastq_prefix + "/SampleSheet.csv.tmp"

        filename = filename.replace('//', '/')
        tmp_filename = tmp_filename.replace('//', '/')

        # SampleSheet should be created only if it does not exist
        # to ensure bcl2fastq is not redone only because the SampleSheet is touched
        # by the process.
        # Maybe replace this with a hash check to still write a new SampleSheet if the
        # infos from Sequencing_summary are different.
        #if not os.path.isfile(filename):
        os.makedirs(bcl2fastq_prefix, exist_ok = True)
        header=('[Header]\n'
                'IEMFileVersion,4\n'
                'Investigator Name,' + current_df[['Customer']].iloc[0,0]+'\n'
                'Experiment Name,' + current_df[['Sample_Project']].iloc[0,0]+'\n'
                'Date,' + str(current_df[['Date_run']]) + '\n'
                'Workflow,GenerateFASTQ\n'
                'Application,NextSeq FASTQ Only\n'
                'Description,' + current_df[['Description']].iloc[0,0]+'\n'
                'Chemistry,Default\n\n'
                '[Reads]\n'
                #+current_df[['Reads']].iloc[0,0].replace('x','\n')+'\n\n'
                # Add str convertion to the read value in order to prevent error is SE (i.e x is absent)
                +str(current_df[['Reads']].iloc[0,0]).replace('x','\n')+'\n\n')
                # Add settings section containing adapter information
                #'[Settings]\n'
                #'Adapter,' + adapter_dict[kit_used][0] + '\n'
                #'AdapterRead2,' + adapter_dict[kit_used][1] + '\n\n'
                #'[Data]\n')

        if(['scRNA-seq', 'scRNA_HTO', 'Cellplex'] not in current_df.Type.unique()):
            adapter_content=('[Settings]\n'
         'Adapter,' + adapter_dict[kit_used][0] + '\n'
         'AdapterRead2,' + adapter_dict[kit_used][1] + '\n\n'
         '[Data]\n')
        else:
            adapter_content=("\n"
            '[Data]\n')
        if current_df[['Process']].iloc[0,0] == "yes":
            myfile = open(filename, 'w')
            myfile.write(header)
            myfile.write(adapter_content)
            myfile.write(data.to_csv(index=False))
            myfile.close()

            #print(str(os.path.isfile(filename)) + " SampleSheet.csv " + filename)
            #print(str(os.path.isfile(tmp_filename)) + " SampleSheet.csv.tmp " + tmp_filename )
            #if(os.path.isfile(filename)):
                #print(hashlib.md5(open(filename,'rb').read()).hexdigest())
                #print(hashlib.md5(open(tmp_filename,'rb').read()).hexdigest())
                #print(hashlib.md5(open(filename,'rb').read()).hexdigest() == hashlib.md5(open(tmp_filename,'rb').read()).hexdigest())
                #print(hashlib.md5(open(filename,'rb').read()).hexdigest() != hashlib.md5(open(tmp_filename,'rb').read()).hexdigest())
                #if(hashlib.md5(open(filename,'rb').read()).hexdigest() != hashlib.md5(open(tmp_filename,'rb').read()).hexdigest()):
                    #os.replace(filename, tmp_filename)
                #else:
                    #if(os.path.isfile(tmp_filename)):
                        #print(tmp_filename)
                        #print(os.path.isfile(tmp_filename))
                        #os.remove(tmp_filename)

            #else:
                #if(os.path.isfile(tmp_filename)):
                    #os.rename(tmp_filename, filename)

        # To avoid id, double slashes "//" are replaced by single one.
        bcl2fastq_target = bcl2fastq_target.replace('//', '/')
        # Add a condition to read the process column of the summary.
        # Sometimes, we may not want to run bcl2fastq.
        if current_df[['Process']].iloc[0,0] == "yes":
            mwconf['bcl2fastq_targets'].append(bcl2fastq_target)
