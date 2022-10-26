# Sequencing summary reader
# This script aims is to create all the needed target for the TGML workflow.
# It is based on the mapper from Salvatore Spicuglia's mapper

if not os.path.isfile("../mw-tgml/Sequencing_summary.xlsx"):
    #eprint('No Sequencing_summary.xlsx')
    mwconf['targets'] = []
    mwconf['bcl2fastq_targets'] = []
    mwconf['qc_targets'] = []
else:
    samples = pandas.read_excel("../mw-tgml/Sequencing_summary.xlsx", sheet_name="samples", engine='openpyxl')
    samples.Type.fillna("Unknown", inplace=True)

    if "Analysis_type" not in samples:
        samples['Analysis_type'] = "default"

    mwconf['targets'] = []
    if 'ids' not in mwconf:
        mwconf['ids'] = {}

    id_chip_qc_to_cat = []
    id_chip_qc_fingerprint_to_cat = []
    id_multiqc_tgml = []
    ln_bam_list_all_chip_atac = []

    for index, row in samples.iterrows():
        SAMPLE_NAME = str(row['Sample_Name'])
        PROCESS = str(row['Process'])

        # If samples have to be processed, get all the needed info from summary
        if PROCESS in ['yes']:
            TYPE = str(row['Type']).replace("-seq", "")
            SE_OR_PE = str(row['Se_or_Pe'])
            SAMPLE_NAME = str(row['sample_name'])
            EXP = str(row['Exp'])
            PROJECT = str(row['Sample_Project'])
            CELL_TYPE = str(row['Cell_type'])
            CUSTOMER = str(row['Customer'])
            ACCESSION = str(row['Accession'])
            ANALYSIS_TYPE = str(row['Analysis_type'])

            if numpy.isnan(row['Run']):
                RUN = 'nan'
            else:
                RUN = str(int(row['Run']))

            if pandas.isna(ACCESSION):
                eprint('Sample ' + SAMPLE_NAME + ' should have a valid accession to be processed.')

            # (1
            # Standardize fastq depending on library layout
            # fq_to_rename is created for single-end
            # fq_to_rename_1 and _2 are created for paired-end
            # for any sources

            # New output from NextSeq500 after Windows 10 update (September 2020)
            elif str(row['Origin']) == 'NS500_W10':
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

            # Old output from NextSeq500. Keep it in case it is needed
            elif str(row['Origin']) == 'NextSeq500':
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

            # Process bcl files from NextSeq500. Output are NOT splitted by lane and fastq files for indexes are generated
            elif(str(row['Origin']) in ['bcl', 'bcl_no_mismatch'] and str(row['Type']) not in ['scRNA', 'scRNA_HTO', 'Cellplex']):
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

            # Add case for scrna_bcl
            elif(str(row['Origin']) == 'bcl' and str(row['Type']) in ['scRNA', 'scRNA_HTO', 'Cellplex']):
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
                    mwconf['ids'][id_cat_1] = str(fq_to_cat_1)
                    mwconf['ids'][id_cat_2] = str(fq_to_cat_2)
                    fq_to_rename_1 = "out/cat/" + id_cat_1
                    fq_to_rename_2 = "out/cat/" + id_cat_2

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
                    "all"           : "tgml/all_samples",
                    "by_type_and_run": "tgml/by_type_and_run/" + TYPE + "/run" + RUN,
                    "by_type_and_exp": "tgml/by_type_and_exp/" + TYPE + "/" + EXP,
                    "by_project"     : "tgml/by_project/" + PROJECT,
                    "by_cell_type"   : "tgml/by_cell_type/" + CELL_TYPE,
                    "by_customer"   : "tgml/by_customer/" + CUSTOMER,
                    "by_run" : "tgml/by_run/run" + RUN
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

            # 4)

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
            elif SPECIE in ['Rat', 'Rattus_norvegicus']:
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

                    # Only run this part if analysis_type is Demultiplexage_Concatenation_Quantification_QC
                    if row['Analysis_type'] in ['Demultiplexage_Concatenation_QC', 'Concatenation_QC']:
                        continue
                    else:
                        bowtie2_stem = "bowtie2/" + row['Se_or_Pe'] + "_" + assembly + "/sickle/" + row['Se_or_Pe'] + "_-t_sanger_-q_20/ln/alias/" + fq_stem
                        bowtie2_log_path = "out/" + bowtie2_stem + ".log"
                        aligned_stem = "samtools/index/samtools/sort/samtools/view_sam_to_bam_-q_30/" + bowtie2_stem

                    mw_bam_path = "out/" + aligned_stem + ".bam"
                    mw_bai_path = mw_bam_path + ".bai"
                    mw_idxstat_path = "out/samtools/idxstats/" + aligned_stem + ".idxstat.tsv"
                    mw_bw_path = "out/deepTools/bamCoverage_--binSize_20_--minMappingQuality_0_--normalizeUsing_RPKM/" + aligned_stem + ".bw"

                    if row['Type'] == 'RNA':
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

                        if row['Type'] == 'RNA':
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

                        if PROCESS == 'yes' and row['Type'] not in ['RNA_fq_only', 'ChIP_fq_only', 'scRNA', 'scRNA_HTO', 'Cellplex', 'Demultiplexage_Concatenation_QC']:
                            mwconf['targets'].append(ln_aligned_unspecif_paths)

                    # (7
                    # Targets for files specific of ChIP-like approaches.
                    if row['Type'] in ['ChIP','ATAC','FAIRE','DNASE','MNase'] and row['Analysis_type'] not in ['Demultiplexage_Concatenation_QC', 'Concatenation_QC']:
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
                        if not pandas.isna(row['Input_chip']):
                            mw_bed_broad = "out/macs2/callpeak_--broad_--gsize_" + gsize + "/" + aligned_stem + "_over_" + row['Input_chip'] + "_peaks.bed"
                            mw_bed_narrow = "out/macs2/callpeak_--gsize_" + gsize + "/" + aligned_stem + "_over_" + row['Input_chip'] + "_peaks.bed"

                            for k in aligned_stem_dict.keys():
                                bed_broad_suffix = aligned_stem_dict[k] + "/bed/broad/" + SAMPLE_NAME + "_over_" + row['Input_chip'] + "_peaks.bed"
                                mwconf['ids'][bed_broad_suffix] = mw_bed_broad
                                ln_bed_broad_path = "out/ln/alias/" + bed_broad_suffix

                                bed_narrow_suffix = aligned_stem_dict[k] + "/bed/narrow/" + SAMPLE_NAME + "_over_" + row['Input_chip'] + "_peaks.bed"
                                mwconf['ids'][bed_narrow_suffix] = mw_bed_narrow
                                ln_bed_narrow_path = "out/ln/alias/" + bed_narrow_suffix

                                if PROCESS == 'yes':
                                    mwconf['targets'].append([ln_bed_broad_path, ln_bed_narrow_path])
                        # 8)

                        # (9
                        # Quantile normalization
                        if not pandas.isna(row['Quantile_normalization_name']):
                            mw_danpos_wiq_path = "out/ucsc/wigToBigWig_-clip_chrominfo-" + assembly + "/danpos/wiq_chrominfo-" + assembly + "/danpos/dtriple/ln/alias/tgml/all_samples/" + assembly + "/bam/" + SAMPLE_NAME + "_qnorVS_" + row['Quantile_normalization_name'] + ".bw"
                            for k in aligned_stem_dict.keys():
                                id_suffix = aligned_stem_dict[k] + "/bw/quantile_normalized/" + row['Sample_Name'] + "_over_" + row['Quantile_normalization_name'] + ".bw"
                                mwconf['ids'][id_suffix] = mw_danpos_wiq_path
                                ln_path = "out/ln/alias/" + id_suffix

                                if PROCESS == 'yes':
                                    mwconf['targets'].append(ln_path)
                        # 9)
    # (10
    # Add here treatment by exp for RNA-Seq
    # Only works if exp (column Q) as been specified
    rna_exps = samples[(samples['Process'].isin(['yes'])) & samples['Type'].isin(['RNA']) & (samples['Exp'] != '')].Exp.unique()
    rna_exps_to_process = samples[(samples['Process'].isin(['yes'])) & samples['Type'].isin(['RNA']) & (samples['Exp'] != '')].Exp.unique()

    for rna_exp in rna_exps:
        rna_exp_samples = samples[(samples['Process'].isin(['yes','done'])) & (samples['Exp'] == rna_exp)]

        if len(rna_exp_samples.specie.unique()) != 1:
            eprint('More than one specie for this RNA experiment ' + str(rna_exp) + '. Analysis steps involving all the samples from this experiment are skipped. There is likely an error in your Sequencing_summary.xlsx')
        else:
            SPECIE = rna_exp_samples.specie.unique()
            if SPECIE in ['human', 'Human', 'Homo_sapiens']:
                assemblies = ["GRCh38", "hg19"]
            elif SPECIE in ['mouse', 'Mouse', 'Mus_musculus']:
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
                bam_paths = str(["out/ln/alias/tgml/all_samples/" + assembly + "/bam/" + sample + ".bam" for sample in rna_exp_samples.sample_name])
                mwconf['ids'][bam_id] = bam_paths
                mwconf['ids'][bam_list_id] = bam_paths
                
                for norm in ["raw", "rpkm"]:
                    mw_path = "out/r/tidy_featureCounts/subread/featureCounts_-O_-t_exon_-g_merge_gene_id_name_gtf-" + assembly + "-merge-attr-retrieve-ensembl_" + bam_id + "_" + norm + ".tsv"
                    for stem in ["tgml/all_samples/", "tgml/by_type_and_exp/RNA/" + rna_exp + "/"]:
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

    # Processing of single-cell RNA-seq
    scRNA_samples = samples[(samples['Process'].isin(['yes'])) & (samples['Type'] == 'scRNA')].Sample_Project.unique()
    for scRNA_project in scRNA_samples:
        project_samples = samples[(samples['Process'].isin(['yes'])) & (samples['Sample_Project'] == scRNA_project)]
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
                ln_cellranger_count = "out/ln/alias/tgml/by_run/run" + RUN + cellranger_count_prefix
                if PROCESS == 'yes':
                    mwconf['targets'].append(cellranger_count_target)
                    #mwconf['targets'].append(ln_cellranger_count)
            # Add the protected_underscore_line to the scRNA_protectedunderscores.tsv file, which is created on each scRNA-seq processing. 
            # This prevent error when having sample name with space
            # 07-07-2021: doesn't work when multiple scRNA are processed.
            if(str(project_samples['Process'].unique()[0]) == 'yes'):
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
            SCRNA_KIT = str(project_samples['Kit_HTO'].unique()[0])

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
                    if(SCRNA_KIT == "Total_seq_A"):
                        citeseq_output_prefix = "out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_-o_Results_" + match + "/cat/merge-nextseq500-pe/" 
                    elif(SCRNA_KIT == "Total_seq_B"):
                        citeseq_output_prefix = "out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_-trim_10_-o_Results_" + match + "/cat/merge-nextseq500-pe/"
                    else:
                        citeseq_output_prefix = "out/cite-seq_count/_-cbf_1_-cbl_16_-umif_17_-umil_28_-o_Results_" + match + "/cat/merge-nextseq500-pe/"
                    citeseq_output_prefix = citeseq_output_prefix.replace("//", "/")
                    citeseq_output_target = citeseq_output_prefix + "/Results_" + match + "/run_report.yaml"


    tgml_subtrees = ["out/ln/alias/tgml/all_samples/"] +\
    glob.glob("out/ln/alias/tgml/by_project/*/") +\
    glob.glob("out/ln/alias/tgml/by_customer/*/") +\
    glob.glob("out/ln/alias/tgml/by_cell_type/*/") +\
    glob.glob("out/ln/alias/tgml/by_type_and_exp/*/*/") +\
    glob.glob("out/ln/alias/tgml/by_type_and_run/*/*/") +\
    glob.glob("out/ln/alias/tgml/by_run/*/")

    multiqc_targets = [stem + "multiqc_report.html" for stem in tgml_subtrees]
    multiqc_targets_interactive = [sub.replace("out/", "out/multiqc/dir_--interactive/") for sub in multiqc_targets]
    multiqc_targets_flat = [sub.replace("out/", "out/multiqc/dir_--flat/") for sub in multiqc_targets]

    md5sum_targets = [sub.replace("out/", "out/find/md5sum/") for sub in tgml_subtrees]
    md5sum_targets = [stem + "md5sum.txt" for stem in md5sum_targets]

    mwconf['qc_targets'] = multiqc_targets_interactive + multiqc_targets_flat + md5sum_targets
