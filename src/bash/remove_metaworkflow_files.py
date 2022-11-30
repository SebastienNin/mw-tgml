# Sequencing summary reader
# This script aims is to create all the needed target for the TGML workflow.
# It is based on the mapper from Salvatore Spicuglia's mapper

import sys, os, argparse, pandas, glob, itertools

import functools
import operator
regular_list = []

import socket

import shutil

current_cluster = socket.gethostname()

if current_cluster == "sacapus":
    mw_path = "/gpfs/tgml/"
elif current_cluster == "bigmemorix":
    mw_path = "/data/TGML"
else:  
    mw_path = "/gpfs/tgml"

# Transform irregular 2D list into a regular one.
def transform(nested_list):
    for ele in nested_list:
        if type(ele) is list:
            regular_list.append(ele)
        else:
            regular_list.append([ele])
    return regular_list

parser = argparse.ArgumentParser(description='This script remove principal targets (bcl2fastq, fastq, fastQC and multiQC output for specified run numbers. You can give this script a list of Runs. For example "python remove_metaworkflow_files_for_a_run.py 200 203 499 301"')
parser.add_argument('run_number', metavar='N', type=int, nargs='+',
                    help='List of Run numbers to process')

args = parser.parse_args()

#if not os.path.isfile("../new_mw-tgml/Sequencing_summary.xlsx"):
if not os.path.isfile(mw_path + "mw-tgml/Sequencing_summary.xlsx"):
    print("No Sequencing_summary found. Please check that the Sequencing_summary.xlsx file exists in mw-sst folder.")
else:
    samples = pandas.read_excel(mw_path + "mw-tgml/Sequencing_summary.xlsx", sheet_name="samples", engine='openpyxl')
    #samples = pandas.read_excel("../new_mw-tgml/Sequencing_summary.xlsx", sheet_name="samples", engine='openpyxl')

FILE_TO_REMOVE = []

for run in args.run_number:
    for index, row in samples.iterrows():
        SAMPLE_NAME = str(row['sample_name'])
        PROCESS = str(row['Process'])
        if PROCESS in ['yes','done']:            
            if pandas.isna(row['Run']):
                continue
            if row['Run'] == run:
                RUN = str(int(row['Run']))
                TYPE = str(row['Type'])
                SE_OR_PE = str(row['Se_or_Pe'])
                SAMPLE_NAME = str(row['sample_name'])
                EXP = str(row['Exp'])
                PROJECT = str(row['Sample_Project'])
                CELL_TYPE = str(row['Cell_Type'])
                CUSTOMER = str(row['Customer'])
                ACCESSION = str(row['Accession'])
                ANALYSIS_TYPE = str(row['Analysis_type'])

                # New output from NextSeq500 after Windows 10 update (September 2020)
                if str(row['Origin']) == 'NS500_W10':
                    if SE_OR_PE == 'se':
                        fq_to_cat = [ACCESSION + "/" + SAMPLE_NAME.replace("_", "-") + "_S" + str(int(row['Sample_Well'])) + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                        #fq_to_cat = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                        id_cat = "merge-nexsteq500-se/" + SAMPLE_NAME + '.fastq.gz'
                        fq_to_rename = "out/cat/" + id_cat

                    if SE_OR_PE == 'pe':
                        #fq_to_cat_1 = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                        #fq_to_cat_2 = [ACCESSION + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                        fq_to_cat_1 = [ACCESSION + "/" + SAMPLE_NAME.replace("_", "-") + "_S" + str(int(row['Sample_Well'])) + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                        fq_to_cat_2 = [ACCESSION + "/" + SAMPLE_NAME.replace("_", "-") + "_S" + str(int(row['Sample_Well'])) + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                        id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                        id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                        fq_to_rename_1 = "out/cat/" + id_cat_1
                        fq_to_rename_2 = "out/cat/" + id_cat_2

                # Old output from NextSeq500. Keep it in case it is needed
                elif str(row['Origin']) == 'NextSeq500':
                    if SE_OR_PE == 'se':
                        fq_to_cat = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                        id_cat = "merge-nexsteq500-se/" + SAMPLE_NAME + '.fastq.gz'
                        fq_to_rename = "out/cat/" + id_cat

                    if SE_OR_PE == 'pe':
                        fq_to_cat_1 = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                        fq_to_cat_2 = [ACCESSION + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                        id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                        id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                        fq_to_rename_1 = "out/cat/" + id_cat_1
                        fq_to_rename_2 = "out/cat/" + id_cat_2

                # Process bcl files from NextSeq500. Output are NOT splitted by lane and fastq files for indexes are generated
                elif(str(row['Origin']) in ['bcl', 'bcl_no_mismatch'] and str(row['type']) not in ['scRNA-seq', 'cellplex', 'snRNA-seq']):
                    if(str(row['Origin']) in ['bcl']):
                        bcl_prefix = "out/bcl2fastq/_--no-lane-splitting_--create-fastq-for-index-reads/" + ACCESSION + "/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                    else:
                        bcl_prefix = "out/bcl2fastq/_--no-lane-splitting_--barcode-mismatches_0/" + ACCESSION + "/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                    bcl_prefix = bcl_prefix.replace('//','/')

                    if SE_OR_PE == 'se':
                        fq_to_rename = bcl_prefix + "_R1_001.fastq.gz"
                    elif SE_OR_PE == 'pe':
                        fq_to_rename_1 = bcl_prefix + "_R1_001.fastq.gz"
                        fq_to_rename_2 = bcl_prefix + "_R2_001.fastq.gz"

                    FILE_TO_REMOVE.append(bcl_prefix)
                    FILE_TO_REMOVE.append(os.path.dirname(os.path.dirname(os.path.dirname(bcl_prefix))))

                # Add case for scrna_bcl
                elif(str(row['Origin']) == 'bcl' and str(row['type']) in ['scRNA','cellplex', 'snRNA-seq']):
                    bcl_prefix = "out/cellranger/mkfastq/" + ACCESSION + "/" + "_".join(str(row['Sample_Project']).split("_")[0:2]) + "/outs/fastq_path/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"])  + "_S" + str(int(row["Sample_Well"]))
                    bcl_prefix = bcl_prefix.replace('//','/')

                    # scRNA-seq can't be se, if SE_OR_PE == 'se', throw an error
                    if SE_OR_PE == 'se':
                        sys.exit("scRNA-seq can't be single-end! Please check the excel file.")
                    elif SE_OR_PE == 'pe':
                        fq_to_cat_1 = [bcl_prefix + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                        fq_to_cat_2 = [bcl_prefix + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                        id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                        id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                        FILE_TO_REMOVE.append(fq_to_cat_1)
                        FILE_TO_REMOVE.append(fq_to_cat_2)
                        fq_to_rename_1 = "out/cat/" + id_cat_1
                        fq_to_rename_2 = "out/cat/" + id_cat_2

                    FILE_TO_REMOVE.append(bcl_prefix)

                elif str(row['Origin']) == 'sra':
                    INPUT_FASTQ_PREFIX = "out/sra-tools/fastq-dump_" + SE_OR_PE + "/" + ACCESSION
                    if SE_OR_PE == 'se':
                        fq_to_rename = INPUT_FASTQ_PREFIX + ".fastq.gz"

                    elif SE_OR_PE == 'pe':
                        fq_to_rename_1 = INPUT_FASTQ_PREFIX + "_1.fastq.gz"
                        fq_to_rename_2 = INPUT_FASTQ_PREFIX + "_2.fastq.gz"

                # 1)
                # (2
                # aliases are created for fastq in order to gather all samples in the same directory
                # as well as other trees defined in the base_stem_dict dict below
                base_stem_dict = {
                        "all"           : "sst/all_samples",
                        "by_type_and_run": "sst/by_type_and_run/" + TYPE + "/run" + RUN,
                        "by_type_and_exp": "sst/by_type_and_exp/" + TYPE + "/" + EXP,
                        "by_project"     : "sst/by_project/" + PROJECT,
                        "by_cell_type"   : "sst/by_cell_type/" + CELL_TYPE,
                        "by_customer"   : "sst/by_customer/" + CUSTOMER,
                        "by_run" : "sst/by_run/run" + RUN
                        }

                fq_stem_dict = {}
                for k in base_stem_dict.keys():
                    fq_stem_dict[k] = base_stem_dict[k] + '/fastq/' + SAMPLE_NAME

                fq_stem = fq_stem_dict['all']
                trim_stem = "sickle/" + SE_OR_PE + "_-t_sanger_-q_20/ln/alias/" + fq_stem

                for k in fq_stem_dict.keys():
                    if SE_OR_PE == 'se':
                        fq_suffix = fq_stem_dict[k] + ".fastq.gz"
                        FILE_TO_REMOVE.append(fq_to_rename)
                        fq_path = "out/ln/alias/" + fq_suffix
                        if PROCESS == 'yes':
                            FILE_TO_REMOVE.append(fq_path)

                    elif SE_OR_PE == 'pe':
                        fq_suffix_1 = fq_stem_dict[k] + "_1.fastq.gz"
                        fq_suffix_2 = fq_stem_dict[k] + "_2.fastq.gz"
                        FILE_TO_REMOVE.append(fq_to_rename_1)
                        FILE_TO_REMOVE.append(fq_to_rename_2)
                        fq_path_1 = "out/ln/alias/" + fq_suffix_1
                        fq_path_2 = "out/ln/alias/" + fq_suffix_2

                        if PROCESS == 'yes':
                            FILE_TO_REMOVE.append(fq_path_1)
                            FILE_TO_REMOVE.append(fq_path_2)
                # 2)

                # (3
                # Here could be a good spot to run fastqc and fastq_screen on se or pe from all_samples tree
                #trim_stem
                if SE_OR_PE == "se":
                    fastqc_path =  "out/fastqc/fastq.gz/ln/alias/" + fq_stem + "_fastqc.zip"
                    fastq_screen_path = "out/fastq_screen/filter/" + trim_stem + "_screen.txt"
                    qc_paths = [fastqc_path, fastq_screen_path]
                elif SE_OR_PE == "pe":
                    fastqc_path_1 = "out/fastqc/fastq.gz/ln/alias/" + fq_stem + "_1_fastqc.zip"
                    fastqc_path_2 = "out/fastqc/fastq.gz/ln/alias/" + fq_stem + "_2_fastqc.zip"
                    fastq_screen_path_1 = "out/fastq_screen/filter/" + trim_stem + "_1_screen.txt"
                    fastq_screen_path_2 = "out/fastq_screen/filter/" + trim_stem + "_2_screen.txt"
                    qc_paths = [fastqc_path_1, fastqc_path_2, fastq_screen_path_1, fastq_screen_path_2]

                if PROCESS == 'yes':
                    FILE_TO_REMOVE.append(qc_paths)
                # 3)

                # (4
                # QC log are symlinked in the alternative trees in order
                # to provide a mechanism to get multiQC reports with only
                # subsets of samples for each subtree.
                for k in base_stem_dict.keys():
                    if SE_OR_PE == "se":
                        fastqc_suffix = base_stem_dict[k] + "/logs/fastqc/" + SAMPLE_NAME + "_fastqc.zip"
                        FILE_TO_REMOVE.append(fastqc_path)
                        ln_fastqc_path = "out/ln/alias/" + fastqc_suffix

                        fastq_screen_suffix = base_stem_dict[k] + "/logs/fastq_screen/" + SAMPLE_NAME + "_screen.txt"
                        FILE_TO_REMOVE.append(fastq_screen_path)
                        ln_fastq_screen_path = "out/ln/alias/" + fastq_screen_suffix

                        ln_qc_paths = [ln_fastqc_path, ln_fastq_screen_path]

                    if SE_OR_PE == "pe":
                        fastqc_suffix_1 = base_stem_dict[k] + "/logs/fastqc/" + SAMPLE_NAME + "_1_fastqc.zip"
                        FILE_TO_REMOVE.append(fastqc_path_1)
                        ln_fastqc_path_1 = "out/ln/alias/" + fastqc_suffix_1

                        fastqc_suffix_2 = base_stem_dict[k] + "/logs/fastqc/" + SAMPLE_NAME + "_2_fastqc.zip"
                        FILE_TO_REMOVE.append(fastqc_path_2)
                        ln_fastqc_path_2 = "out/ln/alias/" + fastqc_suffix_2

                        fastq_screen_suffix_1 = base_stem_dict[k] + "/logs/fastq_screen/" + SAMPLE_NAME + "_1_screen.txt"
                        FILE_TO_REMOVE.append(fastq_screen_path_1)
                        ln_fastq_screen_path_1 = "out/ln/alias/" + fastq_screen_suffix_1

                        fastq_screen_suffix_2 = base_stem_dict[k] + "/logs/fastq_screen/" + SAMPLE_NAME + "_2_screen.txt"
                        FILE_TO_REMOVE.append(fastq_screen_path_2)
                        ln_fastq_screen_path_2 = "out/ln/alias/" + fastq_screen_suffix_2

                        ln_qc_paths = [ln_fastqc_path_1, ln_fastqc_path_2, ln_fastq_screen_path_1, ln_fastq_screen_path_2]

                    if PROCESS == 'yes':
                        FILE_TO_REMOVE.append(ln_qc_paths)

                # (5
                # Assemblies for each specie keyword should be defined here
                SPECIE = str(row['Specie'])
                if SPECIE in ['human', 'Human', 'Homo_sapiens']:
                    assembly_list = ["GRCh38", "hg19"]
                    gsize = "hs"
                    scrna_assembly = "GRCh38-2020-A"
                elif SPECIE in ['mouse', 'Mouse', 'Mus_musculus']:
                    #assembly_list = ["GRCm38", "mm9", "GRCm38-merge-attr", "GRCm38-merge-attr-retrieve"]
                    #2021-03-31 Remove GRCm3[89]-merge-attr, cause redondant with GRCm3[89]-merge-attr-retrieve. The results are the same if using gtftk to retrieve the gtf from ensembl or getting it from the ensembl ftp
                    #2021-04-20: Remove GRCm39 and GRCm39-merge-attr-retriev. Not used by the scientific community for now.
                    #assembly_list = ["GRCm38", "mm9", "GRCm38-merge-attr-retrieve", "GRCm39", "GRCm39-merge-attr-retrieve"]
                    assembly_list = ["GRCm38", "mm9"]
                    gsize = "mm"
                    scrna_assembly = "mm10-2020-A"
                elif SPECIE in ['drosophila', 'Fruit_fly', 'Drosophila_melanogaster']:
                    assembly_list = ["BDGP6"]
                    gsize = "dm"
                elif SPECIE in ['Yeast', 'Saccharomyces_cerevisiae']:
                    assembly_list = ["R64-1-1"]
                    gsize = "12e6"
                elif SPECIE in ['Rat']:
                    assembly_list = ["Rnor6"]
                elif not pandas.isna(SPECIE):
                    assembly_list = [SPECIE]
                else:
                    continue
                # 5)

                # Processing of single-cell RNA-seq
                scRNA_samples = samples[(samples['Process'].isin(['yes','done'])) & (samples['Type'].isin(['scRNA-seq', 'Cellplex', 'snRNA-seq'])].Sample_Project.unique()
                for scRNA_project in scRNA_samples:
                    project_samples = samples[(samples['Process'].isin(['yes','done'])) & (samples['Sample_Project'] == scRNA_project)]
                    # Run cellranger count on cellranger mkfastq output
                    if project_samples['Analysis_type'].all() in ['Demultiplexage_Concatenation_Quantification_QC']:
                        # Loop on sample to run cellranger count on each sample
                        SAMPLES = project_samples['Sample_Name']
                        protected_underscore_line = "cellranger/count\t"
                        samples_to_protect = []
                        for SAMPLE in SAMPLES:
                            samples_to_protect.append(SAMPLE)
                            SPECIE = str(project_samples['Specie'].unique()[0])
                            ACCESSION = str(project_samples['Accession'].unique()[0])
                            PROCESS = str(project_samples['Process'].unique()[0])
                            RUN = str(int(project_samples['Run'].unique()[0]))
                            # Checking specie to process.
                            # For now, on the platform we process only human and mouse.
                            if SPECIE in ['human', 'Human', 'Homo_sapiens']:
                                scrna_assembly = "GRCh38-2020-A"
                            elif SPECIE in ['mouse', 'Mouse', 'Mus_musculus']:
                                scrna_assembly = "mm10-2020-A"
                            else:
                                continue
                            cellranger_count_target = "out/cellranger/count_--sample_" + SAMPLE + "_" + scrna_assembly + "/cellranger/mkfastq" + ACCESSION + "/process_done"
                            # Test to add ln for cellranger
                            cellranger_count_prefix = "/cellranger/count_--sample_" + SAMPLE + "_" + scrna_assembly + "/cellranger/mkfastq" + ACCESSION
                            ln_cellranger_count = "out/ln/alias/sst/by_run/run" + RUN + cellranger_count_prefix
                            if PROCESS == 'yes':
                                FILE_TO_REMOVE.append(cellranger_count_target)

                sst_subtrees = ["out/ln/alias/sst/all_samples/"] +\
                glob.glob("out/ln/alias/sst/by_project/*/") +\
                glob.glob("out/ln/alias/sst/by_customer/*/") +\
                glob.glob("out/ln/alias/sst/by_cell_type/*/") +\
                glob.glob("out/ln/alias/sst/by_type_and_exp/*/*/") +\
                glob.glob("out/ln/alias/sst/by_type_and_run/*/*/") +\
                glob.glob("out/ln/alias/sst/by_run/*/")

                multiqc_targets = [stem + "multiqc_report.html" for stem in sst_subtrees]
                multiqc_targets_interactive = [sub.replace("out/", "out/multiqc/dir_--interactive/") for sub in multiqc_targets]
                multiqc_targets_flat = [sub.replace("out/", "out/multiqc/dir_--flat/") for sub in multiqc_targets]

                md5sum_targets = [sub.replace("out/", "out/find/md5sum/") for sub in sst_subtrees]
                md5sum_targets = [stem + "md5sum.txt" for stem in md5sum_targets]
            
                FILE_TO_REMOVE.append("out/ln/alias/sst/by_run/run" + RUN)
                FILE_TO_REMOVE.append("out/ln/alias/sst/by_type_and_run/" + TYPE + "/run" + RUN)
                FILE_TO_REMOVE.append("out/ln/alias/sst/by_type_and_exp/" + TYPE + "/" + EXP)
                FILE_TO_REMOVE.append("out/ln/alias/sst/by_project/" + PROJECT)

LIST = transform(FILE_TO_REMOVE)
FINAL_LIST = list(set(functools.reduce(operator.iconcat, LIST, [])))

for el in FINAL_LIST:
    el_to_remove = mw_path + "/mw-tgml/" + el
    try:
        os.remove(el_to_remove)
    except OSError:
        pass
    if(os.path.exists(el_to_remove)):
        try:
            shutil.rmtree(el_to_remove)
        except OSError:
            pass
