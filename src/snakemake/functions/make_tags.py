#!/usr/bin/env python3

# This script produces appropriate tags file required to process hashtagged oligos in cell hashing experiments.
# It gets require informations from the Sequencing_summary.xlsx

# Created 25/02/2021

import pandas
import os
import numpy as np
import hashlib
import re
from pathlib import Path

if os.path.isfile("../mw-tgml/Sequencing_summary.xlsx"):
    samples = pandas.read_excel(
        "../mw-tgml/Sequencing_summary.xlsx", sheet_name="samples", engine="openpyxl")

    samples_from_bcl = samples[samples['type'].isin(['scRNA'])]

    experiments_from_bcl = samples_from_bcl.accession.unique()

    #mwconf['bcl2fastq_targets'] = []

    for experiment in experiments_from_bcl:
        # Variable to get the loop dataframe
        current_df = samples_from_bcl[samples['accession'] == experiment]
        if(current_df.analysis_type.unique() == "Demultiplexage_Concatenation_Quantification_QC"):
            # print(current_df)
            current_HTO_list = list(current_df.HTO_information)
            current_ADT_list = list(current_df.ADT_information)

            #print(current_HTO_list)
            #print(current_ADT_list)

            # Identify data where HTO or ADT are presents
            # Naive search of the letter A in string. Maybe find an other way to do it.
            if len(current_HTO_list) > 1:
                HTO_matching = ["A" in str(f) for f in current_HTO_list]
                #print(HTO_matching)
            else:
                continue

            if len(current_ADT_list) > 1:
                ADT_matching = ["A" in str(f) for f in current_ADT_list]
                #print(ADT_matching)
            else:
                continue

            matching_list = []
            if(True in HTO_matching):
                matching_list.append("HTO")
            if(True in ADT_matching):
                matching_list.append("ADT")

            if(len(matching_list)>0):
                # Loop over available information for HTO and/or ADT
                for match in matching_list:
                    info_list = [match, "matching"]
                    info_joined = "_".join(info_list)

                    # Get the line where HTO information is present
                    res = [i for i, val in enumerate(eval(info_joined)) if val]
                    # print(res)
                    info_df = current_df.iloc[res[0]]
                    # Get output folder name
                    # Barcode sequences are different depending on the kit.
                    # TotalSeq A have the Barcode sequence on the first position of the read, whereas TotalseqB having Barcode sequence beginning after the 10 first nucleotides
                    # print(info_df.Kit_barcode_scRNA)
                    if(info_df.Kit == "Total_seq_A"):
                        output_prefix = "out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_-o_Results_" + match + "/" + info_df.accession
                    elif(info_df.Kit == "Total_seq_B"):
                        output_prefix = "out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_-trim_10_-o_Results_" + match + "/" + info_df.accession
                    else:
                        output_prefix = "out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_-o_Results_" + match + "/" + info_df.accession
                    output_prefix = output_prefix.replace("//", "/")

                    output_target = output_prefix + "/Results_" + match + "/run_report.yaml"
                    #print(output_target)

                    splitted_info = eval("info_df." + match + "_information.split(';')")
                    info_to_write = []

                    for info in splitted_info:
                        splitted_info = info.split(",")
                        info_to_write.append(splitted_info[1] + "," + splitted_info[0])

                    """                 # Get the HTO information
                    splitted_HTO_info = info_df.HTO_information.split(";")
                    HTO_to_write = []
                    for info in splitted_HTO_info:
                        splitted_info = info.split(",")
                        HTO_to_write.append(
                            splitted_info[1] + "," + splitted_info[0])

                    # Get the ADT information
                    splitted_ADT_info = info_df.ADT_information.split(";")
                    ADT_to_write = []
                    for info in splitted_ADT_info:
                        splitted_info = info.split(",")
                        ADT_to_write.append(
                            splitted_info[1] + "," + splitted_info[0]) """
                    #tmp_filename = "/data/nin/Workspace/dev/tags.tmp"
                    #filename = "/data/nin/Workspace/dev/tags.csv"
                    tmp_filename = output_prefix + "/tags.tmp"
                    filename = output_prefix + "/tags.csv"

                    # Add 

                    if current_df[['process']].iloc[0, 0] == "yes":
                        os.makedirs(output_prefix, exist_ok=True)
                        myfile = open(tmp_filename, 'w')
                        for el in info_to_write:
                            myfile.write(el + "\n")
                        myfile.close()
                        if(os.path.isfile(filename)):
                            if(hashlib.md5(open(filename, 'rb').read()).hexdigest() != hashlib.md5(open(tmp_filename, 'rb').read()).hexdigest()):
                                os.replace(filename, tmp_filename)
                            else:
                                if(os.path.isfile(tmp_filename)):
                                    os.remove(tmp_filename)
                                # 2021-05-27: Try to correct missing tmp file error by adding a touch if the tmp file doesn't exist.
                                if(os.path.isfile(tmp_filename)):
                                    Path(tmp_filename).touch()
                                    os.remove(tmp_filename)
                                else:
                                    Path(tmp_filename).touch()
                                    os.remove(tmp_filename)
                    else:
                        if(os.path.isfile(tmp_filename)):
                            os.rename(tmp_filename, filename)
