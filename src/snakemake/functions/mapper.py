#if not os.path.isfile("../mw-sst/Sequencing_summary.csv"):
if not os.path.isfile("../mw-sst/Sequencing_summary.xlsx"):
    #eprint('No Sequencing_summary.xlsx')
    mwconf['targets'] = []
    mwconf['bcl2fastq_targets'] = []
    mwconf['qc_targets'] = []
else:
    samples = pandas.read_excel("../mw-sst/Sequencing_summary.xlsx", sheet_name="samples", engine='openpyxl')
    #genomes = pandas.read_excel("../mw-sst/Sequencing_summary.xlsx", sheet_name="genomes")
    #samples = pandas.read_csv("../mw-sst/Sequencing_summary.csv")
    samples.type.fillna("Unknown", inplace=True)

    # TGML users have this column in their sequencing summary
    # but SST users do not.
    if "analysis_type" not in samples:
        samples['analysis_type'] = "default"

    mwconf['targets'] = []
    if 'ids' not in mwconf:
        mwconf['ids'] = {}

    id_chip_qc_to_cat = []
    id_chip_qc_fingerprint_to_cat = []
    id_multiqc_sst = []
    # TO DO: implement a per run, per experiment and per project id_multiqc_* so I can provide report for each level.
    # e.g. id_multiqc_run*, id_multiqc_exp*, id_multiqc_proj*.
    ln_bam_list_all_chip_atac = []

    for index, row in samples.iterrows():
        SAMPLE_NAME = str(row['sample_name'])
        PROCESS = str(row['process'])

        #if PROCESS == 'no':
        #eprint('Skipping ' + SAMPLE_NAME + ' because of process trigger.')
        #else:
        if PROCESS in ['yes','done']:
            #print(row['sample_merge_list'])
            # sample who have are not a merge of other samples are defined here:
            #if pandas.isna(row['sample_merge_list']) and not pandas.isna(row['tgml_fastq_prefix']):
            #if pandas.isna(row['sample_merge_list']):
            TYPE = str(row['type'])
            SE_OR_PE = str(row['se_or_pe'])
            SAMPLE_NAME = str(row['sample_name'])
            EXP = str(row['exp'])
            PROJECT = str(row['Sample_Project'])
            CELL_TYPE = str(row['cell_type'])
            CUSTOMER = str(row['customer'])
            ACCESSION = str(row['accession'])
            ANALYSIS_TYPE = str(row['analysis_type'])


            if numpy.isnan(row['run']):
                RUN = 'nan'
            else:
                RUN = str(int(row['run']))

            if pandas.isna(ACCESSION):
                eprint('Sample ' + SAMPLE_NAME + ' should have a valid accession to be processed.')

            # (1
            # Standardize fastq depending on library layout
            # fq_to_rename is created for single-end
            # fq_to_rename_1 and _2 are created for paired-end
            # for any sources
            if str(row['origin']) == 'tgml':
                TGML_FASTQ_PREFIX = "out/ln/updir/mw/" + ACCESSION
                MERGED_OR_UNMERGED = str(row['merged_or_unmerged'])

            # print(row['sample_name'], row['tgml_fastq_prefix'])
            # path = prefix + stem + suffix
            #base_stem = "sst/" + TYPE + "/run" + RUN
                if SE_OR_PE == 'se':
                    if MERGED_OR_UNMERGED == 'unmerged':
                        fq_to_cat = [TGML_FASTQ_PREFIX + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                        id_cat = "merge-nexsteq500-se/" + SAMPLE_NAME + '.fastq.gz'
                        mwconf['ids'][id_cat] = str(fq_to_cat)
                        fq_to_rename = "out/cat/" + id_cat

                    elif row['merged_or_unmerged'] == 'merged':
                        fq_to_rename = TGML_FASTQ_PREFIX + ".fastq.gz"

                elif SE_OR_PE == 'pe':
                    if MERGED_OR_UNMERGED == 'unmerged':
                        fq_to_cat_1 = [TGML_FASTQ_PREFIX + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                        fq_to_cat_2 = [TGML_FASTQ_PREFIX + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                        id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                        id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                        mwconf['ids'][id_cat_1] = str(fq_to_cat_1)
                        mwconf['ids'][id_cat_2] = str(fq_to_cat_2)
                        fq_to_rename_1 = "out/cat/" + id_cat_1
                        fq_to_rename_2 = "out/cat/" + id_cat_2

                    elif MERGED_OR_UNMERGED == 'merged':
                        fq_to_rename_1 = TGML_FASTQ_PREFIX + "_R1_001.fastq.gz"
                        fq_to_rename_2 = TGML_FASTQ_PREFIX + "_R2_001.fastq.gz"

            # New output form sequencer after Windows 10 update
            elif str(row['origin']) == 'NS500_W10':
                if SE_OR_PE == 'se':
                    fq_to_cat = [ACCESSION + "/" + SAMPLE_NAME.replace("_", "-") + "_S" + str(int(row['Sample_Well'])) + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                    #fq_to_cat = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                    id_cat = "merge-nexsteq500-se/" + SAMPLE_NAME + '.fastq.gz'
                    mwconf['ids'][id_cat] = str(fq_to_cat)
                    fq_to_rename = "out/cat/" + id_cat

                if SE_OR_PE == 'pe':
                    #fq_to_cat_1 = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                    #fq_to_cat_2 = [ACCESSION + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                    fq_to_cat_1 = [ACCESSION + "/" + SAMPLE_NAME.replace("_", "-") + "_S" + str(int(row['Sample_Well'])) + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                    fq_to_cat_2 = [ACCESSION + "/" + SAMPLE_NAME.replace("_", "-") + "_S" + str(int(row['Sample_Well'])) + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                    id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                    id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                    mwconf['ids'][id_cat_1] = str(fq_to_cat_1)
                    mwconf['ids'][id_cat_2] = str(fq_to_cat_2)
                    fq_to_rename_1 = "out/cat/" + id_cat_1
                    fq_to_rename_2 = "out/cat/" + id_cat_2

            elif str(row['origin']) == 'NextSeq500':
                if SE_OR_PE == 'se':
                    fq_to_cat = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                    id_cat = "merge-nexsteq500-se/" + SAMPLE_NAME + '.fastq.gz'
                    mwconf['ids'][id_cat] = str(fq_to_cat)
                    fq_to_rename = "out/cat/" + id_cat

                if SE_OR_PE == 'pe':
                    fq_to_cat_1 = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                    fq_to_cat_2 = [ACCESSION + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                    id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                    id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                    mwconf['ids'][id_cat_1] = str(fq_to_cat_1)
                    mwconf['ids'][id_cat_2] = str(fq_to_cat_2)
                    fq_to_rename_1 = "out/cat/" + id_cat_1
                    fq_to_rename_2 = "out/cat/" + id_cat_2

            elif str(row['origin']) == 'fastq_absolute_path_agilent_XT_HS2':
                    fq_to_cat_1 = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,3)]
                    fq_to_cat_2 = [ACCESSION + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,3)]
                    id_cat_1 = "merge-illumina-2-lanes-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                    id_cat_2 = "merge-illumina-2-lanes-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                    mwconf['ids'][id_cat_1] = str(fq_to_cat_1)
                    mwconf['ids'][id_cat_2] = str(fq_to_cat_2)
                    fq_to_rename_1 = "out/cat/" + id_cat_1
                    fq_to_rename_2 = "out/cat/" + id_cat_2

            # A special case from files from Necker Institute where paired-end reads correspond to R1 and R3
            elif str(row['origin']) == 'NextSeq500_R1_R3':
                fq_to_cat_1 = [ACCESSION + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                fq_to_cat_2 = [ACCESSION + "_L00" + str(n) + "_R3_001.fastq.gz" for n in range(1,5)]
                id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                mwconf['ids'][id_cat_1] = str(fq_to_cat_1)
                mwconf['ids'][id_cat_2] = str(fq_to_cat_2)
                fq_to_rename_1 = "out/cat/" + id_cat_1
                fq_to_rename_2 = "out/cat/" + id_cat_2

            elif str(row['origin']) == 'bcl_NextSeq500':
                # Add /_/ after out/bcl2fastq and to prevent decimal in the sample Well, add int() in the str conversion to first convert into int then into str
                bcl_prefix = "out/bcl2fastq/_/" + ACCESSION + "/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                bcl_prefix = bcl_prefix.replace('//','/')

                if SE_OR_PE == 'se':
                    fq_to_cat = [bcl_prefix + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                    id_cat = "merge-nexsteq500-se/" + SAMPLE_NAME + '.fastq.gz'
                    mwconf['ids'][id_cat] = str(fq_to_cat)
                    fq_to_rename = "out/cat/" + id_cat
                    #mw/out/bcl2fastq/Run_310_NS500-217_25_02_2020_SK/S003693_9_ChIP_H4K5Cro_G5/9_ChIP_H4K5Cro_G5_S9_L001_R1_001.fastq.gz

                if SE_OR_PE == 'pe':
                    fq_to_cat_1 = [bcl_prefix + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                    fq_to_cat_2 = [bcl_prefix + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                    id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                    id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                    mwconf['ids'][id_cat_1] = str(fq_to_cat_1)
                    mwconf['ids'][id_cat_2] = str(fq_to_cat_2)
                    fq_to_rename_1 = "out/cat/" + id_cat_1
                    fq_to_rename_2 = "out/cat/" + id_cat_2

            # Add case for bcl_no_split_lane
            elif(str(row['origin']) in ['bcl', 'bcl_no_mismatch', 'bcl_index_generation'] and str(row['type']) not in ['scRNA', 'scRNA_HTO', 'cellplex']):
                # Add /_/ after out/bcl2fastq and to prevent decimal in the sample Well, add int() in the str conversion to first convert into int then into str
                if(str(row['origin']) == 'bcl'):
                    bcl_prefix = "out/bcl2fastq/_--no-lane-splitting/" + ACCESSION + "/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                elif(str(row['origin']) == 'bcl_index_generation'):
                    bcl_prefix = "out/bcl2fastq/_--no-lane-splitting_--create-fastq-for-index-reads/" + ACCESSION + "/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                    #bcl_prefix = "out/bcl2fastq/_--no-lane-splitting_--create-fastq-for-index-reads_--ignore-missing-bcls/" + ACCESSION + "/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                else:
                    bcl_prefix = "out/bcl2fastq/_--no-lane-splitting_--barcode-mismatches_0/" + ACCESSION + "/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                bcl_prefix = bcl_prefix.replace('//','/')

                if SE_OR_PE == 'se':
                    fq_to_rename = bcl_prefix + "_R1_001.fastq.gz"
                elif SE_OR_PE == 'pe':
                    fq_to_rename_1 = bcl_prefix + "_R1_001.fastq.gz"
                    fq_to_rename_2 = bcl_prefix + "_R2_001.fastq.gz"
                    #id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                    #id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                    #fq_to_cat_1 = "out/ln/alias/sst/all_samples/fastq/" + SAMPLE_NAME + "_1.fastq.gz"
                    #fq_to_cat_2 = "out/ln/alias/sst/all_samples/fastq/" + SAMPLE_NAME + "_2.fastq.gz"
                    #mwconf['ids'][id_cat_1] = str(fq_to_cat_1)
                    #mwconf['ids'][id_cat_2] = str(fq_to_cat_2)


            # Add case for scrna_bcl
            elif(str(row['origin']) == 'bcl' and str(row['type']) in ['scRNA', 'scRNA_HTO', 'cellplex']):
                #bcl_prefix = "out/cellranger/mkfastq/" + ACCESSION + "/" + "_".join(str(row['Sample_Project']).split("_")[0:2]) + "/outs/fastq_path/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                bcl_prefix = "out/cellranger/mkfastq/" + ACCESSION + "/" + "_".join(str(row['Sample_Project']).split("_")[0:2]) + "/outs/fastq_path/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"])  + "_S" + str(int(row["Sample_Well"]))
                #bcl_prefix = "out/cellranger/mkfastq/" + ACCESSION + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                #bcl_prefix = "out/cellranger/mkfastq/" + ACCESSION + "/" + str(row['Sample_Project']) + "/" + str(row["Sample_ID"]) + "/" + str(row["Sample_Name"]) + "_S" + str(int(row["Sample_Well"]))
                bcl_prefix = bcl_prefix.replace('//','/')

                # scRNA-seq can't be se, if SE_OR_PE == 'se', throw an error
                if SE_OR_PE == 'se':
                    sys.exit("scRNA-seq can't be single-end! Please check the excel file.")
                elif SE_OR_PE == 'pe':
                    fq_to_cat_1 = [bcl_prefix + "_L00" + str(n) + "_R1_001.fastq.gz" for n in range(1,5)]
                    fq_to_cat_2 = [bcl_prefix + "_L00" + str(n) + "_R2_001.fastq.gz" for n in range(1,5)]
                    id_cat_1 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_1.fastq.gz'
                    id_cat_2 = "merge-nexsteq500-pe/" + SAMPLE_NAME + '_2.fastq.gz'
                    mwconf['ids'][id_cat_1] = str(fq_to_cat_1)
                    mwconf['ids'][id_cat_2] = str(fq_to_cat_2)
                    fq_to_rename_1 = "out/cat/" + id_cat_1
                    fq_to_rename_2 = "out/cat/" + id_cat_2
                    #fq_to_rename_1 = bcl_prefix + "_L001_R1_001.fastq.gz"
                    #fq_to_rename_2 = bcl_prefix + "_L001_R2_001.fastq.gz"


            elif str(row['origin']) == 'sra':
                INPUT_FASTQ_PREFIX = "out/sra-tools/fastq-dump_" + SE_OR_PE + "/" + ACCESSION
                if SE_OR_PE == 'se':
                    fq_to_rename = INPUT_FASTQ_PREFIX + ".fastq.gz"

                elif SE_OR_PE == 'pe':
                    fq_to_rename_1 = INPUT_FASTQ_PREFIX + "_1.fastq.gz"
                    fq_to_rename_2 = INPUT_FASTQ_PREFIX + "_2.fastq.gz"

            elif str(row['origin']) == 'blueprint':
                INPUT_FASTQ_PREFIX = "out/ln/updir/mw/" + ACCESSION
                if SE_OR_PE == 'se':
                    fq_to_rename = INPUT_FASTQ_PREFIX + ".fastq.gz"

                elif SE_OR_PE == 'pe':
                    fq_to_rename_1 = INPUT_FASTQ_PREFIX + "_1.fastq.gz"
                    fq_to_rename_2 = INPUT_FASTQ_PREFIX + "_2.fastq.gz"

            elif str(row['origin']) == 'fastq_absolute_path':
                if SE_OR_PE == 'se':
                    fq_to_rename = ACCESSION + ".fastq.gz"

                elif SE_OR_PE == 'pe':
                    fq_to_rename_1 = ACCESSION + "_1.fastq.gz"
                    fq_to_rename_2 = ACCESSION + "_2.fastq.gz"

            elif  str(row['origin']) == 'illumina_fastq_sets':
                if SE_OR_PE == 'se':
                    samples_to_cat = str(glob.glob(ACCESSION + "_[0-9][0-9][0-9].fastq.gz"))
                    concat_sample = "merge-illumina-fastq-sets/" + SAMPLE_NAME + ".fastq.gz"
                    mwconf['ids'][concat_sample] = samples_to_cat
                    fq_to_rename = "out/cat/merge-illumina-fastq-sets/" + SAMPLE_NAME + ".fastq.gz"
                    # TODO : add PE HERE LATER


            elif str(row['origin']) == 'merge_fastq':
                SAMPLES_TO_MERGE = ACCESSION.split(",")
                if SE_OR_PE == 'se':
                    samples_to_cat = str(["out/ln/alias/sst/all_samples/fastq/" + SAMPLE + ".fastq.gz" for SAMPLE in SAMPLES_TO_MERGE])
                    concat_sample = "merge-fastq-samples/" + SAMPLE_NAME + ".fastq.gz"
                    mwconf['ids'][concat_sample] = samples_to_cat
                    fq_to_rename = "out/cat/merge-fastq-samples/" + SAMPLE_NAME + ".fastq.gz"

                elif SE_OR_PE == 'pe':
                    samples_to_cat_1 = str(["out/ln/alias/sst/all_samples/fastq/" + SAMPLE + "_1.fastq.gz" for SAMPLE in SAMPLES_TO_MERGE])
                    samples_to_cat_2 = str(["out/ln/alias/sst/all_samples/fastq/" + SAMPLE + "_2.fastq.gz" for SAMPLE in SAMPLES_TO_MERGE])
                    concat_sample_1 = "merge-fastq-samples/" + SAMPLE_NAME + "_1.fastq.gz"
                    concat_sample_2 = "merge-fastq-samples/" + SAMPLE_NAME + "_2.fastq.gz"
                    mwconf['ids'][concat_sample_1] = samples_to_cat_1
                    mwconf['ids'][concat_sample_2] = samples_to_cat_2
                    fq_to_rename_1 = "out/cat/merge-fastq-samples/" + SAMPLE_NAME + "_1.fastq.gz"
                    fq_to_rename_2 = "out/cat/merge-fastq-samples/" + SAMPLE_NAME + "_2.fastq.gz"
            
            elif str(row['origin']) == 'subsample_fastq':
                SUBSAMPLE_ARGS = ACCESSION.split(",")
                SAMPLE=SUBSAMPLE_ARGS.pop(0)
                N_READS_OR_FRACTION = SUBSAMPLE_ARGS.pop(0)
                if len(SUBSAMPLE_ARGS) == 0:
                    SUBSAMPLE_SEED="123"
                else:
                    SUBSAMPLE_SEED = SUBSAMPLE_ARGS.pop(0)
                
                fq_to_rename_prefix = "out/seqtk/sample_-s_" + \
                    SUBSAMPLE_SEED + \
                    "_" + \
                    N_READS_OR_FRACTION + \
                    "/ln/alias/sst/all_samples/fastq/" + \
                    SAMPLE

                if SE_OR_PE == 'se':
                    fq_to_rename = fq_to_rename_prefix + ".fastq.gz"
                
                elif SE_OR_PE == 'pe':
                    fq_to_rename_1 = fq_to_rename_prefix + "_1.fastq.gz"
                    fq_to_rename_2 = fq_to_rename_prefix + "_2.fastq.gz"

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
                    mwconf['ids'][fq_suffix] = fq_to_rename
                    fq_path = "out/ln/alias/" + fq_suffix
                    if PROCESS == 'yes':
                        mwconf['targets'].append(fq_path)

                elif SE_OR_PE == 'pe':
                    fq_suffix_1 = fq_stem_dict[k] + "_1.fastq.gz"
                    fq_suffix_2 = fq_stem_dict[k] + "_2.fastq.gz"
                    mwconf['ids'][fq_suffix_1] = fq_to_rename_1
                    mwconf['ids'][fq_suffix_2] = fq_to_rename_2
                    fq_path_1 = "out/ln/alias/" + fq_suffix_1
                    fq_path_2 = "out/ln/alias/" + fq_suffix_2

                    if PROCESS == 'yes':
                        mwconf['targets'].append(fq_path_1)
                        mwconf['targets'].append(fq_path_2)
            # 2)

            # (3
            # Here could be a good spot to run fastqc and fastq_screen on se or pe from all_samples tree
            #trim_stem
            #fastqc_path_R1 = "out/fastqc/fastq.gz/ln/alias/" + fq_stem + "_1_fastqc.zip"
            #fastqc_path_R2 = "out/fastqc/fastq.gz/ln/alias/" + fq_stem + "_2_fastqc.zip"
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
                mwconf['targets'].append(qc_paths)
            # 3)

            # (4
            # QC log are symlinked in the alternative trees in order
            # to provide a mechanism to get multiQC reports with only
            # subsets of samples for each subtree.
            for k in base_stem_dict.keys():
                if SE_OR_PE == "se":
                    fastqc_suffix = base_stem_dict[k] + "/logs/fastqc/" + SAMPLE_NAME + "_fastqc.zip"
                    mwconf['ids'][fastqc_suffix] = fastqc_path
                    ln_fastqc_path = "out/ln/alias/" + fastqc_suffix

                    fastq_screen_suffix = base_stem_dict[k] + "/logs/fastq_screen/" + SAMPLE_NAME + "_screen.txt"
                    mwconf['ids'][fastq_screen_suffix] = fastq_screen_path
                    ln_fastq_screen_path = "out/ln/alias/" + fastq_screen_suffix

                    ln_qc_paths = [ln_fastqc_path, ln_fastq_screen_path]

                if SE_OR_PE == "pe":
                    fastqc_suffix_1 = base_stem_dict[k] + "/logs/fastqc/" + SAMPLE_NAME + "_1_fastqc.zip"
                    mwconf['ids'][fastqc_suffix_1] = fastqc_path_1
                    ln_fastqc_path_1 = "out/ln/alias/" + fastqc_suffix_1

                    fastqc_suffix_2 = base_stem_dict[k] + "/logs/fastqc/" + SAMPLE_NAME + "_2_fastqc.zip"
                    mwconf['ids'][fastqc_suffix_2] = fastqc_path_2
                    ln_fastqc_path_2 = "out/ln/alias/" + fastqc_suffix_2

                    fastq_screen_suffix_1 = base_stem_dict[k] + "/logs/fastq_screen/" + SAMPLE_NAME + "_1_screen.txt"
                    mwconf['ids'][fastq_screen_suffix_1] = fastq_screen_path_1
                    ln_fastq_screen_path_1 = "out/ln/alias/" + fastq_screen_suffix_1

                    fastq_screen_suffix_2 = base_stem_dict[k] + "/logs/fastq_screen/" + SAMPLE_NAME + "_2_screen.txt"
                    mwconf['ids'][fastq_screen_suffix_2] = fastq_screen_path_2
                    ln_fastq_screen_path_2 = "out/ln/alias/" + fastq_screen_suffix_2

                    ln_qc_paths = [ln_fastqc_path_1, ln_fastqc_path_2, ln_fastq_screen_path_1, ln_fastq_screen_path_2]

                if PROCESS == 'yes':
                    mwconf['targets'].append(ln_qc_paths)


            #for k in base_stem_dict.keys():
            #   for qc_path in qc_paths:
            #       qc_suffix = qc_path.replace('out/', base_stem_dict[k] + "/logs/")
            #       mwconf['ids'][qc_suffix] = qc_path
            #       print('qc_suffix')
            #       print(qc_suffix)
            #       print('qc_path')
            #       print(qc_path)
            #       print("mwconf['ids'][qc_suffix]")
            #       print(mwconf['ids'][qc_suffix])
            #       qc_k_path = "out/ln/alias/" + qc_suffix
            #       if PROCESS == 'yes':
            #           mwconf['targets'].append(qc_k_path)
            # 4)
            #sst/all_samples/logs/fastqc/fastq.gz/ln/alias/sst/all_samples/fastq/RPMI_H3K4me3_fastqc.zip

            # (5
            # Assemblies for each specie keyword should be defined here
            SPECIE = str(row['specie'])
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
            elif SPECIE in ['drosophila', 'Fruit_fly', 'Drosophila_melanogaster']:
                scrna_assembly = "mm10-2020-A"
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

            # (6
            # Targets for post-alignment, experiment-type-unspecific files are produced
            # I do not want to throw an error if specie is not referred and just skip sample instead.
            if 'assembly_list' in locals():
                for assembly in assembly_list:
                    aligned_stem_dict = {}
                    for k in base_stem_dict.keys():
                        aligned_stem_dict[k] = base_stem_dict[k] + "/" + assembly

                    if row['analysis_type'] in ['Demultiplexage_Concatenation_QC', 'Concatenation_QC']:
                        continue
                    elif row['type'] == 'RNA' and row['origin'] != "fastq_absolute_path_agilent_XT_HS2":
                        aligned_stem = "samtools/index/star/" + row['se_or_pe'] + "_fastq.gz_to_bam_standard_staridx-" + assembly + "-merge-attr-retrieve-ensembl_gtf-" + assembly +"-merge-attr-retrieve-ensembl/sickle/" + row['se_or_pe'] + "_-t_sanger_-q_20/ln/alias/" + fq_stem
                    elif row['type'] == 'RNA' and row['origin'] == "fastq_absolute_path_agilent_XT_HS2":
                        aligned_stem = "samtools/index/samtools/sort/agent/locatit_mbc_-i_-R/samtools/sort_-n/star/" + row['se_or_pe'] + "_fastq.gz_to_bam_standard_staridx-" + assembly + "-merge-attr-retrieve-ensembl_gtf-" + assembly +"-merge-attr-retrieve-ensembl/agent/trim_-v2/ln/alias/" + fq_stem
                    else:
                        bowtie2_stem = "bowtie2/" + row['se_or_pe'] + "_" + assembly + "/sickle/" + row['se_or_pe'] + "_-t_sanger_-q_20/ln/alias/" + fq_stem
                        bowtie2_log_path = "out/" + bowtie2_stem + ".log"
                        aligned_stem = "samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/" + bowtie2_stem

                    mw_bam_path = "out/" + aligned_stem + ".bam"
                    mw_bai_path = mw_bam_path + ".bai"
                    mw_idxstat_path = "out/samtools/idxstats/" + aligned_stem + ".idxstat.tsv"
                    mw_bw_path = "out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM/" + aligned_stem + ".bw"

                    if row['type'] == 'RNA':
                        mw_bw_fwd_path = "out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM_--filterRNAstrand_forward/" + aligned_stem + ".bw"
                        mw_bw_rev_path = "out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM_--filterRNAstrand_reverse/" + aligned_stem + ".bw"
                        # 2021-03-25 Add rseqc geneBody_coverage
                        mw_rseqc_path = "out/rseqc/geneBody_coverage_bed-housekeeping-genes-" + assembly + "/" + aligned_stem + ".geneBodyCoverage.curves.pdf"

                    for k in aligned_stem_dict.keys():
                        bam_suffix = aligned_stem_dict[k] + "/bam/" + SAMPLE_NAME + ".bam"
                        mwconf['ids'][bam_suffix] = mw_bam_path
                        ln_bam_path = "out/ln/alias/" + bam_suffix

                        bai_suffix = aligned_stem_dict[k] + "/bam/" + SAMPLE_NAME + ".bam.bai"
                        mwconf['ids'][bai_suffix] = mw_bai_path
                        ln_bai_path = "out/ln/alias/" + bai_suffix

                        idxstat_suffix = aligned_stem_dict[k] + "/logs/samtools_idxstats/" + SAMPLE_NAME + ".idxstat.tsv"
                        mwconf['ids'][idxstat_suffix] = mw_idxstat_path
                        ln_idxstat_path = "out/ln/alias/" + idxstat_suffix

                        bw_suffix = aligned_stem_dict[k] + "/bw/" + SAMPLE_NAME + ".bw"
                        mwconf['ids'][bw_suffix] = mw_bw_path
                        ln_bw_path = "out/ln/alias/" + bw_suffix

                        ln_aligned_unspecif_paths = [ln_bam_path, ln_bai_path, ln_idxstat_path, ln_bw_path]

                        if row['type'] == 'RNA':
                            bw_fwd_suffix = aligned_stem_dict[k] + "/bw/stranded/" + SAMPLE_NAME + "_fwd.bw"
                            mwconf['ids'][bw_fwd_suffix] = mw_bw_fwd_path
                            ln_bw_fwd_path = "out/ln/alias/" + bw_fwd_suffix

                            bw_rev_suffix = aligned_stem_dict[k] + "/bw/stranded/" + SAMPLE_NAME + "_rev.bw"
                            mwconf['ids'][bw_rev_suffix] = mw_bw_rev_path
                            ln_bw_rev_path = "out/ln/alias/" + bw_rev_suffix

                            ln_aligned_unspecif_paths.append(ln_bw_fwd_path)
                            ln_aligned_unspecif_paths.append(ln_bw_rev_path)

                            #rseqc_suffix = aligned_stem_dict[k] + "/pdf/" + SAMPLE_NAME + ".pdf"
                            #mwconf['ids'][rseqc_suffix] = mw_rseqc_path
                            #ln_rseqc_path = "out/ln/alias/" + rseqc_suffix

                            #ln_aligned_unspecif_paths.append(ln_rseqc_path)

                        if PROCESS == 'yes' and row['type'] not in ['RNA_fq_only', 'ChIP_fq_only', 'scRNA', 'scRNA_HTO', 'cellplex', 'Demultiplexage_Concatenation_QC']:
                            mwconf['targets'].append(ln_aligned_unspecif_paths)

                    # (7
                    # Targets for files specific of ChIP-like approaches.
                    if row['type'] in ['ChIP','ATAC','FAIRE','DNASE','MNase'] and row['analysis_type'] not in ['Demultiplexage_Concatenation_QC', 'Concatenation_QC']:
                        mw_chip_qc_fingerprint_prefix = "out/deepTools/plotFingerprint/" + aligned_stem
                        mw_chip_qc_fingerprint_metrics = mw_chip_qc_fingerprint_prefix + ".metrics.tsv"
                        mw_chip_qc_fingerprint_counts = mw_chip_qc_fingerprint_prefix + ".counts.tsv"
                        mw_chip_qc_phantompeakqualtools = "out/phantompeakqualtools/bam_noctrl_-savp/" + aligned_stem + '.spp.out'
                        mw_bed_broad = "out/macs2/noctrl_callpeak_--broad_--gsize_" + gsize + "/" + aligned_stem + "_peaks.bed"
                        mw_xls_broad = "out/macs2/noctrl_callpeak_--broad_--gsize_" + gsize + "/" + aligned_stem + "_peaks.xls"
                        mw_bed_narrow = "out/macs2/noctrl_callpeak_--gsize_" + gsize + "/" + aligned_stem + "_peaks.bed"
                        mw_xls_narrow = "out/macs2/noctrl_callpeak_--gsize_" + gsize + "/" + aligned_stem + "_peaks.xls"

                        for k in aligned_stem_dict.keys():
                            fingerprint_metrics_suffix = aligned_stem_dict[k] + "/logs/fingerprint/" + SAMPLE_NAME + ".metrics.tsv"
                            mwconf['ids'][fingerprint_metrics_suffix] = mw_chip_qc_fingerprint_metrics
                            ln_fingerprint_metrics_path = "out/ln/alias/" + fingerprint_metrics_suffix

                            fingerprint_counts_suffix = aligned_stem_dict[k] + "/logs/fingerprint/" + SAMPLE_NAME + ".counts.tsv"
                            mwconf['ids'][fingerprint_counts_suffix] = mw_chip_qc_fingerprint_counts
                            ln_fingerprint_counts_path = "out/ln/alias/" + fingerprint_counts_suffix

                            phantompeakqualtools_suffix = aligned_stem_dict[k] + "/logs/phantompeakqualtools/" + SAMPLE_NAME + ".spp.out"
                            mwconf['ids'][phantompeakqualtools_suffix] = mw_chip_qc_phantompeakqualtools
                            ln_phantompeakqualtools_path = "out/ln/alias/" + phantompeakqualtools_suffix

                            bed_broad_suffix = aligned_stem_dict[k] + "/bed/broad/" + SAMPLE_NAME + "_peaks.bed"
                            mwconf['ids'][bed_broad_suffix] = mw_bed_broad
                            ln_bed_broad_path = "out/ln/alias/" + bed_broad_suffix

                            bed_narrow_suffix = aligned_stem_dict[k] + "/bed/narrow/" + SAMPLE_NAME + "_peaks.bed"
                            mwconf['ids'][bed_narrow_suffix] = mw_bed_narrow
                            ln_bed_narrow_path = "out/ln/alias/" + bed_narrow_suffix

                            xls_broad_suffix = aligned_stem_dict[k] + "/bed/broad/" + SAMPLE_NAME + "_peaks.xls"
                            mwconf['ids'][xls_broad_suffix] = mw_xls_broad
                            ln_xls_broad_path = "out/ln/alias/" + xls_broad_suffix

                            xls_narrow_suffix = aligned_stem_dict[k] + "/bed/narrow/" + SAMPLE_NAME + "_peaks.xls"
                            mwconf['ids'][xls_narrow_suffix] = mw_xls_narrow
                            ln_xls_narrow_path = "out/ln/alias/" + xls_narrow_suffix

                            ln_aligned_chip_specif_paths = [
                                    ln_fingerprint_metrics_path,
                                    ln_fingerprint_counts_path,
                                    ln_phantompeakqualtools_path,
                                    ln_bed_broad_path,
                                    ln_bed_narrow_path,
                                    ln_xls_broad_path,
                                    ln_xls_narrow_path
                                    ]

                            if PROCESS == 'yes':
                                mwconf['targets'].append(ln_aligned_chip_specif_paths)

                        # 7)

                        # (8
                        # Peak-calling with control:
                        if not pandas.isna(row['control_name']):
                            mw_bed_broad = "out/macs2/callpeak_--broad_--gsize_" + gsize + "/" + aligned_stem + "_over_" + row['control_name'] + "_peaks.bed"
                            mw_bed_narrow = "out/macs2/callpeak_--gsize_" + gsize + "/" + aligned_stem + "_over_" + row['control_name'] + "_peaks.bed"

                            for k in aligned_stem_dict.keys():
                                bed_broad_suffix = aligned_stem_dict[k] + "/bed/broad/" + SAMPLE_NAME + "_over_" + row['control_name'] + "_peaks.bed"
                                mwconf['ids'][bed_broad_suffix] = mw_bed_broad
                                ln_bed_broad_path = "out/ln/alias/" + bed_broad_suffix

                                bed_narrow_suffix = aligned_stem_dict[k] + "/bed/narrow/" + SAMPLE_NAME + "_over_" + row['control_name'] + "_peaks.bed"
                                mwconf['ids'][bed_narrow_suffix] = mw_bed_narrow
                                ln_bed_narrow_path = "out/ln/alias/" + bed_narrow_suffix

                                if PROCESS == 'yes':
                                    mwconf['targets'].append([ln_bed_broad_path, ln_bed_narrow_path])
                        # 8)

                        # (9
                        # Quantile normalization
                        if not pandas.isna(row['quantile_normalization_name']):
                            ##
                            # TODO: Check if sample and quantile_normalization share the same library type : se or pe.
                            # Test: bam2wig as a replacement to danpos dtriple for the wig generation
                            # bioconductor-epigenomix may be a substitution too.
                            #if row['se_or_pe'] == "pe":
                            #   PAIRED="1"
                            #elif row['se_or_pe'] == "se":
                            #   PAIRED="0"
                            #else:
                            #   print('se_or_pe should be either pe or se')
                            #mw_danpos_wiq_path = "out/ucsc/wigToBigWig_-clip_chrominfo-" + assembly + "/danpos/wiq_chrominfo-" + assembly + "/danpos/dtriple_-m_" + PAIRED + "/" + aligned_stem + "_qnorVS_" + row['quantile_normalization_name'] + ".bw"
                            # the code above is working for most usecases but can not work if the sample and reference are from different library types se VS pe.
                            # The quickest workaround is to always treat samples as single end for danpos, although it comes with cost to accuracy:
                            # Note this is important to start from aliased bam, instead of fastq, to have se and pe samples in the same directory
                            mw_danpos_wiq_path = "out/ucsc/wigToBigWig_-clip_chrominfo-" + assembly + "/danpos/wiq_chrominfo-" + assembly + "/danpos/dtriple/ln/alias/sst/all_samples/" + assembly + "/bam/" + SAMPLE_NAME + "_qnorVS_" + row['quantile_normalization_name'] + ".bw"

                            #mw_danpos_wiq_path = "out/ucsc/wigToBigWig_-clip_chrominfo-hg38/danpos/wiq_chrominfo-hg38/jvarkit/bam2wig/" + aligned_stem + "_qnorVS_" + row['quantile_normalization_name'] + ".bw"
                            #mw_danpos_wiq_path = "out/ucsc/wigToBigWig_-clip_chrominfo-" + assembly + "/danpos/wiq_chrominfo-" + assembly + "/danpos/dtriple/" + aligned_stem + "_qnorVS_" + row['quantile_normalization_name'] + ".bw"
                            #mw_danpos_wiq_path = "out/ucsc/wigToBigWig_-clip_chrominfo-hg38/danpos/wiq_chrominfo-hg38/danpos/dtriple_--width_10_--distance_100_--edge_1_--paired_0/ln/alias/" + aligned_stem + "/bam/" + row["sample_name"] + "_qnorVS_" + row['quantile_normalization_name'] + ".bw"
                            #mw_danpos_wiq_path = "out/ucsc/wigToBigWig_-clip_chrominfo-hg38/danpos/wiq_chrominfo-hg38/jvarkit/bam2wig/ln/alias/" + aligned_stem + "/bam/" + row['sample_name'] + "_qnorVS_" + row['quantile_normalization_name'] + ".bw"

                            for k in aligned_stem_dict.keys():
                                id_suffix = aligned_stem_dict[k] + "/bw/quantile_normalized/" + row['sample_name'] + "_over_" + row['quantile_normalization_name'] + ".bw"
                                mwconf['ids'][id_suffix] = mw_danpos_wiq_path
                                ln_path = "out/ln/alias/" + id_suffix

                                if PROCESS == 'yes':
                                    mwconf['targets'].append(ln_path)
                        # 9)
    # (10
    # Add here treatment by exp for RNA-Seq
    #samples_rna_exp = samples[(samples['type'] == 'RNA') & (samples['process'] in ['yes','done'])]
    # First
    #rna_exps = samples[(samples['type'] == 'RNA') & (samples['process'] == 'yes') & (samples['exp'] != '')].exp.unique()

    rna_exps = samples[(samples['process'].isin(['yes','done'])) & samples['type'].isin(['RNA']) & (samples['exp'] != '')].exp.unique()

    rna_exps_to_process = samples[(samples['process'].isin(['yes'])) & samples['type'].isin(['RNA']) & (samples['exp'] != '')].exp.unique()

    for rna_exp in rna_exps:
        rna_exp_samples = samples[(samples['process'].isin(['yes','done'])) & (samples['exp'] == rna_exp)]

        if len(rna_exp_samples.specie.unique()) != 1:
            eprint('More than one specie for this RNA experiment ' + str(rna_exp) + '. Analysis steps involving all the samples from this experiment are skipped. There is likely an error in your Sequencing_summary.xlsx')
        else:
            SPECIE = rna_exp_samples.specie.unique()
            if SPECIE in ['human', 'Human', 'Homo_sapiens']:
                assemblies = ["GRCh38", "hg19"]
            elif SPECIE in ['mouse', 'Mouse', 'Mus_musculus']:
                #assemblies = ["GRCm38", "mm9", "GRCm38-merge-attr", "GRCm38-merge-attr-retrieve"]
                # 2021-03-31 Remove GRCm3[89]-merge-attr, cause redondant with GRCm3[89]-merge-attr-retrieve. The results are the same if using gtftk to retrieve the gtf from ensembl or getting it from the ensembl ftp.

                #2021-04-20: Remove GRCm39 and GRCm39-merge-attr-retriev. Not used by the scientific community for now.
                #assemblies = ["GRCm38", "mm9", "GRCm38-merge-attr-retrieve", "GRCm39", "GRCm39-merge-attr-retrieve"]
                assemblies = ["GRCm38", "mm9"]
            elif SPECIE in ['drosophila', 'Fruit_fly', 'Drosophila_melanogaster']:
                assemblies = ["BDGP6"]
            elif SPECIE in ['Yeast', 'Saccharomyces_cerevisiae']:
                assemblies = ["R64-1-1"]
            elif SPECIE in ['Rat']:
                assemblies = ["Rnor6"]
            elif not pandas.isna(SPECIE):
                assemblies = [SPECIE]

            for assembly in assemblies:
                bam_id = bam_list_id = "bam-" + assembly + "-" + rna_exp
                # Not sure if '_' should absolutely by replaced by '-'
                bam_id = bam_id.replace('_','-')
                bam_paths = str(["out/ln/alias/sst/all_samples/" + assembly + "/bam/" + sample + ".bam" for sample in rna_exp_samples.sample_name])
                mwconf['ids'][bam_id] = bam_paths
                mwconf['ids'][bam_list_id] = bam_paths
                #sst/by_type_and_exp/" + TYPE + "/exp" + EXP
                #id_suffix =
                # UPDATE (2021-03-15: Add a condition to process RNA-seq using merged attributes gtf thanks to gtftk.)
                # UPDATE (2021-04-21: Remove the condition to process RNA-seq using merged attributes, now all RNA-seq are processed using merged attributes.)
                #if("merge-attr" in assembly):
                    #g_param = "merge_gene_id_name"
                #else:
                    #g_param = "gene_id"
                # TODO: FIND WAY TO FACTORISE aligned_stem_dict so I do not have to redo it again here.
                for norm in ["raw", "rpkm"]:
                    #mw_path = "out/r/tidy_featureCounts/subread/featureCounts_-O_-t_exon_-g_gene_id_gtf-" + assembly + "-ensembl_" + bam_id + "_" + norm + ".tsv"
                    mw_path = "out/r/tidy_featureCounts/subread/featureCounts_-O_-t_exon_-g_merge_gene_id_name_gtf-" + assembly + "-merge-attr-retrieve-ensembl_" + bam_id + "_" + norm + ".tsv"
                    for stem in ["sst/all_samples/", "sst/by_type_and_exp/RNA/" + rna_exp + "/"]:
                        id_suffix = stem + assembly + "-merge-attr-retrieve/counts/" + rna_exp + "_" + norm + ".tsv"
                        mwconf['ids'][id_suffix] = mw_path

                        if rna_exp in rna_exps_to_process:
                            ln_path = "out/ln/alias/" + id_suffix
                            mwconf['targets'].append(ln_path)

                # 2021-03-30 Add a condition to know if rseqc needs to be run
                # 2021-04-21 Run only on GRCh38 and GRCm38
                if assembly in ['GRCh38', 'GRCm38']:
                    # 2021-03-29 Add rseqc output for bam list
                    bam_list_id = bam_list_id.replace('_', '-')
                    rseqc_genecov_path = "rseqc/geneBody_coverage_bam_list_bed-housekeeping-genes-" + assembly + "/" + bam_list_id + ".geneBodyCoverage.curves.pdf"
                    ln_rseqc_genecov_path = "out/" + rseqc_genecov_path
                    # 2021-03-30 Add rseqc tin.py
                    rseqc_tin_path = "rseqc/tin_bam_list_bed-housekeeping-genes-" + assembly + "/" + bam_list_id + ".summary.txt"
                    ln_rseqc_tin_path = "out/" + rseqc_tin_path

                    if rna_exp in rna_exps_to_process:
                        mwconf["targets"].append(ln_rseqc_genecov_path)
                        mwconf["targets"].append(ln_rseqc_tin_path)

    projects = samples[(samples['process'].isin(['yes','done'])) & (samples['project'] != '')].project.unique()
    for project in projects:
        project_samples = samples[(samples['process'].isin(['yes','done'])) & (samples['project'] == project)]

        if len(project_samples.specie.unique()) != 1:
            eprint('More than one specie for this project ' + str(project_samples) + '. There is likely an error in your Sequencing_summary.xlsx.')
        else:
            SPECIE = str(row['specie'])
            if SPECIE in ['human', 'Human', 'Homo_sapiens']:
                assemblies = ["GRCh38", "hg19"]
            elif SPECIE in ['mouse', 'Mouse', 'Mus_musculus']:
                assemblies = ["GRCm38", "mm9"]
            elif SPECIE in ['drosophila', 'Fruit_fly', 'Drosophila_melanogaster']:
                assemblies = ["BDGP6"]
            elif SPECIE in ['Yeast', 'Saccharomyces_cerevisiae']:
                assemblies = ["R64-1-1"]
            elif not pandas.isna(SPECIE):
                assemblies = [SPECIE]

            for assembly in assemblies:
                bam_id = "bam-" + assembly + "-project-" + project
                # Not sure if '_' should absolutely by replaced by '-'
                bam_id = bam_id.replace('_','-')
                bam_paths = str(["out/ln/alias/sst/all_samples/" + assembly + "/bam/" + sample + ".bam" for sample in project_samples.sample_name])
                mwconf['ids'][bam_id] = bam_paths

                bed_broad_id = "bed-broad-" + assembly + "-project-" + project
                bed_broad_id_with_ext = bed_broad_id + ".bed"
                bed_broad_paths = str(["out/ln/alias/sst/all_samples/" + assembly + "/bed/broad/" + sample + "_peaks.bed" for sample in project_samples.sample_name])
                mwconf['ids'][bed_broad_id_with_ext] = bed_broad_paths

                merged_bed_broad_path = str(["out/bedtools/merge/sort/_-k1,1_-k2,2n/cat/" + bed_broad_id_with_ext])
                #print(merged_bed_broad_path)
                merged_bed_broad_id = "bed-merged-broad-" + assembly + "-project-" + project
                mwconf['ids'][merged_bed_broad_id] = merged_bed_broad_path




                #mw_path = "out/subread/featureCounts_-O_-t_exon_-g_gene_id_gtf-" + assembly + "-ensembl_" + bam_id
                #id_suffix =
                # TODO: FIND WAY TO FACTORISE aligned_stem_dict so I do not have to redo it again here.
                #for k in aligned_stem_dict.keys():
                #                id_suffix = aligned_stem_dict[k] + "/bw/quantile_normalized/" + row['sample_name'] + "_over_" + row['quantile_normalization_name'] + ".bw"
                #                mwconf['ids'][id_suffix] = mw_danpos_wiq_path
                #                ln_path = "out/ln/alias/" + id_suffix





    #experiments_from_bcl = samples_from_bcl.accession.unique()

    #mwconf['bcl2fastq_targets'] = []

    #for experiment in experiments_from_bcl:

    # 10)

            # (7 Specific process for scRNA-seq.
            # This part must contains process for hto also!

            # First, merge all fastq from scRNA. Should only merge mRNA data. HTO data are processed later.

    #scRNA_samples = samples[(samples['process'].isin(['yes','done'])) & samples['type'] == 'scRNA'].Sample_Project.unique()
    scRNA_samples = samples[(samples['process'].isin(['yes','done'])) & (samples['type'] == 'scRNA')].Sample_Project.unique()
    for scRNA_project in scRNA_samples:
        project_samples = samples[(samples['process'].isin(['yes','done'])) & (samples['Sample_Project'] == scRNA_project)]

        # Run cellranger count on cellranger mkfastq output
        if project_samples['analysis_type'].all() in ['Demultiplexage_Concatenation_Quantification_QC']:
            # Loop on sample to run cellranger count on each sample
            SAMPLES = project_samples['Sample_Name']
            #print(SAMPLES)
            protected_underscore_line = "cellranger/count\t"
            samples_to_protect = []
            for SAMPLE in SAMPLES:
                samples_to_protect.append(SAMPLE)
                #print(SAMPLE)
                SPECIE = str(project_samples['specie'].unique()[0])
                ACCESSION = str(project_samples['accession'].unique()[0])
                PROCESS = str(project_samples['process'].unique()[0])
                RUN = str(int(project_samples['run'].unique()[0]))
                # Checking specie to process.
                # For now, on the platform we process only human and mouse.
                if SPECIE in ['human', 'Human', 'Homo_sapiens']:
                    scrna_assembly = "GRCh38-2020-A"
                elif SPECIE in ['mouse', 'Mouse', 'Mus_musculus']:
                    scrna_assembly = "mm10-2020-A"
                else:
                    continue
                #cellranger_count_target = "out/ln/alias/sst/by_run/run" + RUN +  "/cellranger/count_--sample_"+ SAMPLE + "_" + scrna_assembly + "/cellranger/mkfastq" + ACCESSION + "/process_done"
                cellranger_count_target = "out/cellranger/count_--sample_" + SAMPLE + "_" + scrna_assembly + "/cellranger/mkfastq" + ACCESSION + "/process_done"
                #cellranger_count_target = "out/cellranger/count_" + scrna_assembly + "/cellranger/mkfastq" + ACCESSION + "web_summary.html"
                #cellranger_count_target = cellranger_count_target.replace("//", "/")
                # Test to add ln for cellranger
                cellranger_count_prefix = "/cellranger/count_--sample_" + SAMPLE + "_" + scrna_assembly + "/cellranger/mkfastq" + ACCESSION
                ln_cellranger_count = "out/ln/alias/sst/by_run/run" + RUN + cellranger_count_prefix
                if PROCESS == 'yes':
                    mwconf['targets'].append(cellranger_count_target)
                    #mwconf['targets'].append(ln_cellranger_count)
                    #print(cellranger_count_target)
            # Add the protected_underscore_line to the scRNA_protectedunderscores.tsv file, which is created on each scRNA-seq processing.
            # 07-07-2021: doesn't work when multiple scRNA are processed.
            if(str(project_samples['process'].unique()[0]) == 'yes'):
                #print(protected_underscore_line)
                protected_underscore_line = protected_underscore_line + ",".join(samples_to_protect)
                #print(protected_underscore_line)
                protected_underscore_file = open("out/mw/Run_" + RUN  + "_protectedunderscores.tsv", 'w')
                protected_underscore_file.write("tool\tprotectedunderscore")
                protected_underscore_file.write("\n")
                protected_underscore_file.write(protected_underscore_line)
                protected_underscore_file.write("\n")
                protected_underscore_file.close()


            
            # 2021-07-05: Add CITE-seq and cellhashing processing
            ADT = list(project_samples['ADT_information'])
            HTO = list(project_samples['HTO_information'])
            SCRNA_KIT = str(project_samples['Kit_barcode_scRNA'].unique()[0])

            # Identify data where HTO or ADT are presents
            # Naive search of the letter A in string. Maybe find an other way to do it.
            if len(HTO) > 1:
                HTO_matching = ["A" in str(f) for f in HTO]
                #print(HTO_matching)
            else:
                continue

            if len(ADT) > 1:
                ADT_matching = ["A" in str(f) for f in ADT]
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

                    # Get output folder name
                    # Barcode sequences are different depending on the kit.
                    # TotalSeq A have the Barcode sequence on the first position of the read, whereas TotalseqB having Barcode sequence beginning after the 10 first nucleotides
                    # print(info_df.Kit_barcode_scRNA)
                    if(SCRNA_KIT == "Total_seq_A"):
                        citeseq_output_prefix = "out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_-o_Results_" + match + "/cat/merge-nextseq500-pe/" 
                    elif(SCRNA_KIT == "Total_seq_B"):
                        citeseq_output_prefix = "out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_-trim_10_-o_Results_" + match + "/cat/merge-nextseq500-pe/"
                    else:
                        citeseq_output_prefix = "out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_-o_Results_" + match + "/cat/merge-nextseq500-pe/"
                    citeseq_output_prefix = citeseq_output_prefix.replace("//", "/")

                    #mwconf['ids'][citeseq_output_prefix] = "out/" + cellranger_count_prefix + "/Run_" + str(int(project_samples['run'].unique())) + "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" 

                    #print(mwconf['ids'][citeseq_output_prefix])

                    citeseq_output_target = citeseq_output_prefix + "/Results_" + match + "/run_report.yaml"
                    #print(citeseq_output_target)
                    #SAMPLE1 = "/cat/merge-nexsteq500-pe/" + SAMPLE + "_1.fastq.gz"
                    #SAMPLE2 = "/cat/merge-nexsteq500-pe/" + SAMPLE + "_2.fastq.gz"
                    #mwconf['ids'][citeseq_output_prefix] = SAMPLE1
                    #mwconf['ids'][citeseq_output_prefix] = SAMPLE2
                    #print(str(project_samples['Sample_Name']))

                    #if PROCESS == 'yes':
                        #mwconf['targets'].append(citeseq_output_target)
    
    ##      # Dealing with samples and merge of samples gathered into an experiment:
    ##      # Currently working on it for Capstarrseq:
    ##      if row['exp'] != '':
    ##          print('sample_name "' + row['sample_name'] + '" in exp "' + row['exp'] +'"')
    ##          if row['type'] == 'CapStarr' and row['specie'] == 'mouse':
    ##              path_bam_in = "out/ln/sst_exp/" + row['exp'] + "/mm9/" + row['sample_name'] + ".bam"
    ##              path_bam_ln = "sst/experiments/" + row['exp'] + "/mm9/" + row['sample_name'] + ".bam"
    ##              d_paths_in_to_ln[path_bam_in] = path_bam_ln
    ##
    ##              path_bai_in = path_bam_in + ".bai"
    ##              path_bai_ln = path_bam_ln + ".bai"
    ##              d_paths_in_to_ln[path_bai_in] = path_bai_ln
    ##
    ##              path_bw_in="out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM_--extendReads_314/" + in_aligned_suffix + ".bw"
    ##              path_bw_ln = "sst/experiments/" + row['exp'] + "/mm9/" + row['sample_name'] + ".bw"
    ##              d_paths_in_to_ln[path_bw_in] = path_bw_ln
    ##
    ##              if row['control_name'] not in ['','irrelevant']:
    ##                  path_tsv_in = 'out/capstarrseq/merge_all_data_mTDHS/awk/extend_reads_314/awk/keep_first_mate_for_pe_bedtools_bamtobed/bedtools/bamtobed/ln/sst_exp/' + row['exp'] + "/mm9/" + row['sample_name'] + "_over_" + row['control_name'] + ".allData.tsv"
    ##                  path_tsv_ln = "sst/experiments/" + row['exp'] + "/mm9/" + row['sample_name'] + "_over_" + row['control_name'] + ".allData.tsv"
    ##                  d_paths_in_to_ln[path_tsv_in] = path_tsv_ln
    ##
    ##                  path_pdf_in = 'out/capstarrseq/grouping_crms/capstarrseq/fold_change_mTDHS/awk/extend_reads_314/awk/keep_first_mate_for_pe_bedtools_bamtobed/bedtools/bamtobed/ln/sst_exp/' + row['exp'] + "/mm9/" + row['sample_name'] + "_over_" + row['control_name'] + ".inflexionPointGroups.pdf"
    ##                  path_pdf_ln = "sst/experiments/" + row['exp'] + "/mm9/" + row['sample_name'] + "_over_" + row['control_name'] + ".inflexionPointGroups.pdf"
    ##                  d_paths_in_to_ln[path_pdf_in] = path_pdf_ln
    ##
    ##  # This file will log links done between mw and sst
    ##  if usage == 'run':
    ##      file = open("out/ln/sst_table.tsv","w")
    ##      for path_in in d_paths_in_to_ln.keys():
    ##          path_ln = d_paths_in_to_ln[path_in]
    ##          shell("mkdir -p `dirname {path_ln}`; ln -f {path_in} {path_ln}")
    ##          file.write(path_in + "  " + path_ln)
    ##      file.close()
    ##
    ##  return list(d_paths_in_to_ln.keys())
    #
    #### function to create an ID table for bam to merge:
    ###for row in csv.DictReader(open(TSV),delimiter='  '):
    ### paths = list()
    ### print('Creating IDs for samples to merge')
    ### if row['sample_merge_list'] != "":
    ###     samples = row['sample_merge_list'].split(",")
    ###     for sample in samples:
    ###         for row2 in csv.DictReader(open(TSV),delimiter='    '):
    ###             if row2['sample_name'] == sample:
    ###             paths.append(mwconf[sample]['path_bam_in'])
    ### return(paths)
    ##
    ##def input_bam_samtools_merge_samples_sst(wildcards):
    ##  """
    ##  Created:
    ##      2018-01-01 22:29:32
    ##  Aim:
    ##  """
    ##  paths = list()
    ##  assembly = wildcards['assembly']
    ##
    ##  for row in csv.DictReader(open(TSV),delimiter=' '):
    ##
    ##      # Only look for samples that have available files
    ##      #if ',' in row['sample_merge_list']:
    ##      #if row['exp'] == wildcards['exp']:···
    ##      if row['sample_name'] == wildcards['sample_name']:
    ##          if row['exp'] != wildcards['exp']:
    ##              raise Exception("This sample can't be in this experiment based on metadata.")
    ##          #print("sample_name : " + row['sample_name'])
    ##          samples = row['sample_merge_list'].split(",")
    ##          for sample in samples:
    ##              #print('sample: ' + sample)
    ##              for row2 in csv.DictReader(open(TSV),delimiter='    '):
    ##                  if row2['sample_name'] == sample:
    ##                      path_bam_in = "out/samtools/sort/samtools/view_bSh/bowtie2/" + row2['se_or_pe'] + "_" + assembly +"/sickle/" + row2['se_or_pe'] + "_-t_sanger_-q_20/" + "gunzip/merge_lanes_nextseq500_" + row2['se_or_pe'] + "_raw/ln/updir/mw/" + row2['tgml_fastq_prefix'] + ".bam"
    ##                      print('bam: ' + path_bam_in)
    ##                      paths.append(path_bam_in)
    ##
    ##  return(paths)
    ##
    ##def input_bam_ln_sst_exp(wildcards):
    ##  """
    ##  Created:
    ##      2018-01-29 14:54:54
    ##  Aim:
    ##  """
    ##  assembly = wildcards['assembly']
    ##  exp = wildcards['exp']
    ##  sample_name = wildcards['sample_name']
    ##
    ##  for row in csv.DictReader(open(TSV),delimiter=' '):
    ##      if row['sample_name'] == wildcards['sample_name']:
    ##          if exp not in [row['exp'],'ChIP']:
    ##              # 'ChIP' added here because I want all ChIP-Seq samples to be put together in a 'ChIP' type experiment, mainly because it is easier then to deal with peak callin with control samples.
    ##              raise Exception("This sample can't be in this experiment based on metadata.")
    ##          if ',' in row['sample_merge_list']:
    ##              path = "out/samtools/merge_samples_sst/" + exp + "/" + assembly + "/" + sample_name + ".bam"
    ##          else:
    ##              path = "out/samtools/sort/samtools/view_bSh/bowtie2/" + row['se_or_pe'] + "_" + assembly +"/sickle/" + row['se_or_pe'] + "_-t_sanger_-q_20/" + "gunzip/merge_lanes_nextseq500_" + row['se_or_pe'] + "_raw/ln/updir/mw/" + row['tgml_fastq_prefix'] + ".bam"
    ##
    ##  return(path)
    ##
    ##print(ln_bam_list_all_chip_atac)
    ##mwconf['ids']['bam-sst-all-chip-atac'] = str(ln_bam_list_all_chip_atac)
    ##mwconf['targets'].append('out/deepTools/plotFingerprint_bam-sst-all-chip-atac.txt')
    #
    #mwconf['ids']['multiqc-sst'] = id_multiqc_sst
    #mwconf['targets'].append('out/multiqc/req_multiqc-sst/multiqc_report.html')
    #mwconf['targets'].append(id_multiqc_sst)

    #with open('out/multiQC_test_file_list.txt', 'w+') as f:
    #   for item in id_multiqc_sst:
    #       f.write("%s\n" % item)

    #mwconf['ids']['sst/fingerprint.tsv'] = str(id_chip_qc_fingerprint_to_cat)
    #mwconf['targets'].append('out/sort/_-u/cat/sst/fingerprint.tsv')
    #
    #mwconf['ids']['sst/phantompeak.tsv'] = str(id_chip_qc_to_cat)
    #mwconf['targets'].append('out/cat/sst/phantompeak.tsv')
    #print(mwconf['targets'])
    #print(id_chip_qc_to_cat)
    #print(str(id_chip_qc_to_cat))

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

    mwconf['qc_targets'] = multiqc_targets_interactive + multiqc_targets_flat + md5sum_targets
    #eprint(mwconf['targets'])
