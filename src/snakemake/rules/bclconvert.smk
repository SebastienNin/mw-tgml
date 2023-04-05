rule bclconvert:
    """
    Created:
        2020-05-18 00:37:53
    Aim:
        Produce fastq from bcl
        Note the fastq are not explicitely defined here so you should first have a
        Snakemake workflow requiring the tree.html file for all your runs, which will also produce all fastq files.
        Then use a second workflow that knows where the produced fastq are.
        Please open a github issue if you know how to properly explicitely define fastq in this case.
        Maybe checkpoints ?

    Test:
        out/bcl2fastq/_--barcode-mismatches_1/gpfs/projects/spicuglia/mw/inp/bcl/Run_310_200225_NS500637_0217_AH2JLTBGXF/Reports/html/tree.html out/bcl2fastq/_--barcode-mismatches_1_--no-lane-splitting/gpfs/projects/spicuglia/mw/inp/bcl/Run_310_200225_NS500637_0217_AH2JLTBGXF/Reports/html/tree.html
    """
    input:
        xml="/{filler}/RunInfo.xml", 
        csv="out/{tool}{extra}/{filler}/SampleSheet.csv",
    output:
        xml="out/{tool}{extra}/{filler}/RunInfo.xml",
        txt="out/{tool}{extra}/{filler}/Logs/FastqComplete.txt"
    params:
        extra = params_extra
    wildcard_constraints:
        tool="bcl-convert/"
    threads:
        8
    log:
        bclconvert_log = "out/{tool}{extra}/{filler}/log",
    shell:
        """
        cp {input.xml} {output.xml}
        INDIR=`dirname {input.xml}`
        OUTDIR=`dirname {output.xml}`
        (bcl-convert --bcl-input-directory $INDIR --output-directory $OUTDIR --sample-sheet {input.csv} {params.extra}) &> {log.bclconvert_log}
        """
