#!/bin/bash

# loops through a list of bedfiles (chunks.txt) and substitutes a string in a template script (mpileup_template.sh), creates a new slurm job script with the substitution, and submits the job.

for i in `cat chunks.txt`; do sed "s/chunk/$i/g" mpileup_template.sh > $i.mpile.sh; sbatch $i.mpile.sh; done