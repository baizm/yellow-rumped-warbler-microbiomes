#!/bin/bash

# loops through a list of individual sample names (indivs.txt) and substitutes a string in a template script (pear-bwa_template.sh), creates a new slurm job script with the substitution, and submits the job.

for i in `cat indivs.txt`; do sed "s/indiv/$i/g" pear-bwa-template.sh > $i.pear-bwa.sh; sbatch $i.pear-bwa.sh; done