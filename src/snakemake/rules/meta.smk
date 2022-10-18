localrules: config_targets, config_qc_targets, config_blc2fastq_targets
rule config_targets:
    input: mwconf['targets']

rule config_qc_targets:
    input: mwconf['qc_targets']

rule config_blc2fastq_targets:
    input: mwconf['bcl2fastq_targets']

