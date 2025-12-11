# used to merge bams for sample replicates

for indiv in `cat indivs.txt`; do grep $indiv yrwabams.txt > temp.bams.txt; \
varA=`awk "NR==1" temp.bams.txt`;
varB=`awk "NR==2" temp.bams.txt`;
sed "s/fileA/$varA/g" mergetemplate.sh | \
sed "s/fileB/$varB/g" | \
sed "s/indiv/$indiv/g" > ./mergescripts/$indiv.merge.sh; \
done


#merge the two bams
samtools merge ./mergedbams/indiv.merged.bam fileA fileB -@ 24

#index the final bam
samtools index ./mergedbams/indiv.merged.bam
