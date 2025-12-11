#!/bin/bash

# break up the genome into intervals 
bedtools makewindows -g ~/shared/bwaDb/mywagenomev2.1.fa -w 20000000 > genome_chunks.bed

for i in {1..69}; do
    sed -n "${i}p" genome_chunks.bed > chunk${i}.bed
done
