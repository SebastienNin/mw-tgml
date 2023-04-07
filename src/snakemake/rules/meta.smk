localrules: config_targets, config_qc_targets, config_bcl2fastq_targets
rule config_targets:
    input: mwconf['targets']

rule config_qc_targets:
    input: mwconf['qc_targets']

rule config_bcl2fastq_targets:
    input: mwconf['bcl2fastq_targets']

rule config_bclconvert_targets:
    input: mwconf['bclconvert_targets']
