module use /storage/group/dut374/default/sw/modules
module load all

##or /storage/group/dut374/default/sw/common to get other softwares
## make sure echo looks exactly as below to prevent echo from skipping the first line & bash "event not found" error

cd /storage/group/dut374/default/YRWA_GBS_2023/

##Demultiplex using GBSX## 
#make sure: -Xmx___g is asking for enough memory, -kc is true (retain cut site), -t is number of threads#
echo "#! /bin/bash
#SBATCH --account=dut374_c
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=36:00:00
#SBATCH --mem=100GB
java -Xmx150g -jar /storage/group/dut374/default/sw/common/bin/GBSX_v1.3.jar --Demultiplexer -f1 /storage/group/dut374/default/YRWA_GBS_2023/fastq/YRWA_GBS_S1_R1_001.fastq.gz -f2 /storage/group/dut374/default/YRWA_GBS_2023/fastq/YRWA_GBS_S1_R2_001.fastq.gz -i /storage/group/dut374/default/YRWA_GBS_2023/4plate_barcode.txt -rad false -gzip true -kc true -t 15 -o /storage/group/dut374/default/YRWA_GBS_2023/demultiplex_new/" >> gbs_demultiplex.sbatch
sbatch gbs_demultiplex.sbatch

##AdapterRemoval## 
#create the /collapsed folder first, then run script#
cd /storage/group/dut374/default/YRWA_GBS_2023/collapsed
echo "#! /bin/bash
#SBATCH --account=dut374_c
#SBATCH --nodes=1
#SBATCH --ntasks=4 
#SBATCH --time=1-5:00:00
for gfile in /storage/group/dut374/default/YRWA_GBS_2023/demultiplex_new/*.R1.fastq.gz;
	do /storage/group/dut374/default/sw/adapterremoval-2.1.7/bin/AdapterRemoval --file1 $gfile --file2 ${gfile/.R1.fastq.gz/.R2.fastq.gz} --collapse --trimns --minlength 20 --qualitybase 33 --basename /storage/group/dut374/default/YRWA_GBS_2023/collapsed/${gfile##*/}
done" >> adapterremoval.sbatch
sbatch adapterremoval.sbatch

##BowTie2## 
module use /storage/group/dut374/default/sw/modules
module load all
##Align## remember to make /BAM/SAM folder#
cd /storage/group/dut374/default/YRWA_GBS_2023/
readpair_array=( $( find collapsed -name "*.fastq.gz.collapsed" ))
for j in ${readpair_array[@]}
do
i=$( echo $j | sed 's/.fastq.gz.collapsed//')
sample=$(basename $i)

echo "#! /bin/bash
#SBATCH --account=dut374_c
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=1-5:00:00
bowtie2 -p 1 --very-sensitive-local --local -N 0 --phred33 -x /storage/group/dut374/default/mywa_genome_2/final_assembly/mywagenomev2.1 --rg-id $sample --rg SM:$sample -1 /storage/group/dut374/default/YRWA_GBS_2023/collapsed/${sample}.fastq.gz.pair1.truncated -2 /storage/group/dut374/default/YRWA_GBS_2023/collapsed/${sample}.fastq.gz.pair2.truncated -U /storage/group/dut374/default/YRWA_GBS_2023/collapsed/${sample}.fastq.gz.collapsed -X 700 -S /storage/group/dut374/default/YRWA_GBS_2023/BAM/SAM/${sample}.sam >& /storage/group/dut374/default/YRWA_GBS_2023/BAM/SAM/bowtie.${sample}.log" >> $sample.sbatch
sbatch $sample.sbatch
done
#move scripts to /scripts after job is finished; otherwise it will break the for loop & can't queue all using $sample.sbatch
#removed undetermined file, size ~ 120GB

###NOTE 10/25/2023: SINCE WE DON'T MARK DUPPLICATES, MOVE EVERYTHING TO /BAM/SAM AND DELETE /unmarked_bam BEFORE PROCEEDING###
##Rename file paths above for reproducible scripts

##SamTools##
module use /storage/group/dut374/default/sw/modules
module load all
##Make BAM from SAM##
cd /storage/group/dut374/default/YRWA_GBS_2023/BAM/
readpair_array=( $( find SAM -name "*.sam" ))
for j in ${readpair_array[@]}
do
i=$( echo $j | sed 's/.sam//')
sample=$(basename $i)

echo "#! /bin/bash
#SBATCH --account=dut374_c
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=25:00:00

samtools view -bS /storage/group/dut374/default/YRWA_GBS_2023/BAM/SAM/$sample.sam > /storage/group/dut374/default/YRWA_GBS_2023/BAM/$sample.bam" >> $sample.bam.sbatch
sbatch $sample.bam.sbatch
done
#move script to /scripts

##Sort BAM##
module use /storage/group/dut374/default/sw/modules
module load all
#don't try to change readpair_array into ".bam"--it will duplicate the file names and doesn't work
cd /storage/group/dut374/default/YRWA_GBS_2023/BAM/
readpair_array=( $( find SAM -name "*.sam" ))
for j in ${readpair_array[@]}
do
i=$( echo $j | sed 's/.sam//')
sample=$(basename $i)

echo "#! /bin/bash
#SBATCH --account=dut374_c
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=25:00:00
#SBATCH --mem=50GB

samtools sort /storage/group/dut374/default/YRWA_GBS_2023/BAM/$sample.bam -T /storage/group/dut374/default/YRWA_GBS_2023/BAM/$sample.temp.bam -o /storage/group/dut374/default/YRWA_GBS_2023/BAM/$sample.sorted.bam" >> $sample.sorted.sbatch
sbatch $sample.sorted.sbatch
done
#move script to /scripts

##To check if BAM files are sorted## can add -n X after 'head' which specify how many lines can be viewed
cd /storage/group/dut374/default/YRWA_GBS_2023/BAM/
samtools view -H CF26D06.R1.sorted.bam | head

##Index BAM##
module use /storage/group/dut374/default/sw/modules
module load all

cd /storage/group/dut374/default/YRWA_GBS_2023/BAM/
readpair_array=( $( find SAM -name "*.sam" ))
for j in ${readpair_array[@]}
do
i=$( echo $j | sed 's/.sam//')
sample=$(basename $i)

echo "#! /bin/bash
#SBATCH --account=dut374_c
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00

samtools index /storage/group/dut374/default/YRWA_GBS_2023/BAM/$sample.sorted.bam /storage/group/dut374/default/YRWA_GBS_2023/BAM/$sample.sorted.bai" >> $sample.index.sbatch
sbatch $sample.index.sbatch
done
#move script to /scripts
