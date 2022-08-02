#!/bin/bash

# Activate a controlled environment for all users.
. /gpfs/tgml/nin/anaconda3/etc/profile.d/conda.sh
conda activate /gpfs/tgml/nin/anaconda3/envs/dev

# Download the configuration table.
cd /gpfs/tgml/mw-tgml
#cp -f /mnt/wintgml/TGML/Sequencing_summary.xlsx .
cp -f /mnt/panoramix/tgml/TGML/Sequencing_summary.xlsx .
#cp -f /gpfs/tgml/nin/Sequencing_summary.xlsx .
chmod 664 ./Sequencing_summary.xlsx
#wget -N https://www.dropbox.com/s/vqffgdxvc9u6end/Sequencing_summary.xlsx


#SNAKEMAKE_ARGS='-prk --rerun-incomplete --cores 16 --resources wget_limit=2 conda_token=1 --cluster "qsub -V -q lifescope -o log/qsub -e log/qsub -l nodes=1:ppn={threads},walltime=99:00:00" --use-conda --max-jobs-per-second 3'
SNAKEMAKE_ARGS='-prk --rerun-incomplete --cores 16 --resources wget_limit=2 conda_token=1 --cluster "qsub -V -q lifescope -o log/qsub -e log/qsub -l nodes=1:ppn={threads},walltime=99:00:00" --use-conda --max-jobs-per-second 3'

# Run continuous integration tests.
# Send an error by mail if one file can not be generated from inputs files and current workflow state.
#cd ../mw-sst-rulegraph
#snakemake --use-conda -f --rerun-incomplete -prk --verbose out/snakemake/stdout_--rulegraph_ci-sst.dot
#echo "${SNAKEMAKE_ARGS}" | xargs snakemake out/snakemake/stdout_--rulegraph_ci-sst.dot
#cd ../mw-sst

## Produces fastq from bcl
echo "${SNAKEMAKE_ARGS}" | xargs snakemake config_blc2fastq_targets 
echo "${SNAKEMAKE_ARGS}" | xargs snakemake config_targets 
echo "${SNAKEMAKE_ARGS}" | xargs snakemake config_qc_targets  

# Create output directories for the edge case where workflow is first executed with no target (maybe not required anymore, to test)
#mkdir -p out/ln/alias/sst/ out/config_targets/out/ln/alias/sst/

# Create a dereferenced tree with multiple hardlinks for the same files.
# This is required to get rsync sync only one time files that are duplicated in multiple trees.
#find out/ln/alias/sst/ -type l -print0 | while read -r -d $'\0' file; do mkdir -p `dirname out/config_targets/$file` ; ln -Lf $file out/config_targets/$file ; done
#find out/ln/alias/sst/ -type l -print0 | while IFS= read -r -d $'\\0' file; do mkdir -p `dirname out/config_targets/$file` ; ln -Lf $file out/config_targets/$file ; done

# Rsync can now send changed files to the NAS
# rsync --recursive --times -H -vPh --stats out/config_targets/out/ln/alias/sst/ /mnt/SynologySpicuglia/guillaume/illumina_sequencing_data
# rsync --recursive out/multiqc/out/ /mnt/SynologySpicuglia/guillaume/illumina_sequencing_data/multiQC
#rsync -zarmvP out/multiqc/ sebastiennin@gcbio.info:~ 

# Huge processed files may now be removed from the output directory.
# Small ones should be conserved in order to have multiqc find the log files of previous runs,
# thus allowing to still have "done" samples in the QC report.
#find out -type f -size +100000k -delete

#chmod -R 774 out
python ../scripts/utilitaires/symbolic_link_from_metaworkflow.py
