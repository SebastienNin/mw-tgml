#!/bin/bash

# Download the configuration table.
cd /mnt/data/TGML/mw-tgml
#chmod 664 ./Sequencing_summary.xlsx

#SNAKEMAKE_ARGS='-prk --rerun-incomplete --cores 16 --resources wget_limit=2 conda_token=1 --cluster "qsub -V -q lifescope -o log/qsub -e log/qsub -l nodes=1:ppn={threads},walltime=99:00:00" --use-conda --max-jobs-per-second 3'
SNAKEMAKE_ARGS='-prk --rerun-incomplete --cores 16 --resources wget_limit=2 conda_token=1 --use-conda --max-jobs-per-second 3'

## Produces fastq from bcl
echo "${SNAKEMAKE_ARGS}" | xargs snakemake config_blc2fastq_targets 
echo "${SNAKEMAKE_ARGS}" | xargs snakemake config_targets 
echo "${SNAKEMAKE_ARGS}" | xargs snakemake config_qc_targets  

#chmod -R 774 out
#python ../scripts/utilitaires/symbolic_link_from_metaworkflow.py
find /mnt/data/TGML/mw-tgml/out/ -user $USER -type d -exec chmod 775 {} \; -exec chgrp tgml {} \; -exec chmod g+s {} \; -exec setfacl -d -m g::rwx {} \; -exec setfacl -d -m u::rwx {} \; -exec setfacl -d -m o::rx {} \; 
find /mnt/data/TGML/mw-tgml/out/ -user $USER -type f -exec chmod u+rw,g+rw,o+r {} \; -exec chgrp tgml {} \;

