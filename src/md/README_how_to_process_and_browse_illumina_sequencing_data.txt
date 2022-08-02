# Browsing Illumina sequencing results.

Everything you need from an already processed run can be found following 'illumina_sequencing_data' link.

# How to process sequencing data from a new TGML run.

1. Copy the raw fastq.gz files from TGML into 'mw/inp/fastq/runXXX' directory. TGML files are shared through a SFTP server. You can download them directly from Sacapus using sftp command line entering the server using "sftp USERNAME@SERVERIP", then "get -r DIRECTORY". Sacapus sftp version is too old to have the '-r' argument, thus you should use a conda environment with the "openssh" package installed.

2. Update the file "Sequencing summary.xlsx" in Dropbox appending metadata from your experiment. tgml_prefix_path should not include the "_L001_R1_001" part. Here is an example of how you can get a list of prefixes ready to paste in this column. "cd /gpfs/projects/spicuglia/mw;  find inp/fastq/runXXX/ -name '*_L001_R1_001.fastq.gz' | sed 's/_L001_R1_001.fastq.gz//'"

3. Launch the analysis: "./process_illumina_sequencing_data.sh".  Error "Directory cannot be locked." is expected if you or someone else already started the analysis. Just wait until the other run is complete and launch it again.

4. Wait until the workflow finishes.

5. Enjoy your nicely ordered processed files in 'mw/sst'. A symbolic link here 'illumina_sequencing_data' will lead you there.

# Side notes.

- Do NOT add anything in "mw/sst" directory. It contains only files produced by the Snakemake workflow and ordered according to Salvatore Spicuglia Tree needs. The tree is removed at each execution of the workflow to keep it tidy up after updates, so any files put here will be lost forever.

- Files in mw/inp are tracked by git-annex (https://git-annex.branchable.com/) so they are only symlinks to annex hidden directory. Common linux tools (cp, ln,...) can take "-L" argument to handle these symlinks like they were files, e.g. if you need to copy input files in your own directory.

- Errors like "chgrp: changing group of `.': Operation not permitted" are expected when the targeted file is not your own. Do not worry and contact the owner if these files are not already in the thymus group and you need write access to them. If the owner is not reachable, make an Astec ticket.

- Please send bug reports or suggestions for improvement to guillaume.charbonnier@outlook.com starting mail subject with "[mw] Bug report:" or "[mw] Suggest:"
