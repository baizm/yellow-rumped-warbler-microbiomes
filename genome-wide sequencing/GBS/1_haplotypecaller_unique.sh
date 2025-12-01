#pick up here with sorted bams and index files
#to re-do joint genotyping after removing replicates 

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/
cd $WORKDIR

#load the necessary modules
module load gcccore/11.2.0
module load gatk/4.3.0.0-Java-11.0.16

readpair_array=( $( find /projects/academic/mbaiz/YRWA_GBS_2023/bams_sorted2 -name "*.sorted.bam" ))
for j in ${readpair_array[@]}
do
i=$( echo $j | sed 's/.sorted.bam//')
sample=$(basename $i)

echo "#!/bin/sh
#
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=50GB
#SBATCH --job-name="slurm_haplotype_caller"
#SBATCH --output=slurm_haplotype_caller.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

gatk HaplotypeCaller -R /projects/academic/mbaiz/mywa_genome/mywagenomev2.1.fa \
  -I /projects/academic/mbaiz/YRWA_GBS_2023/bams_sorted2/$sample.sorted.bam \
  -ERC GVCF \
  -O /vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/$sample.g.vcf >& /vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/$sample.gvcf.log" >> $sample.snpcall.sbatch
sbatch $sample.snpcall.sbatch
done

