#!/bin/bash

# I need to get the largest Hi-C phaseset per chr

# first compress and index vcf
bgzip -c hapcut2.concatenated.phased.vcf > hapcut2.concatenated.phased.vcf.gz && bcftools index hapcut2.concatenated.phased.vcf.gz

# now extract the largest phaseset per chr
for x in chr1 chr2 chr3 chr4 chr5 chr6 chr7
do
    bcftools view -r ${x} hapcut2.concatenated.phased.vcf.gz | bcftools view -i "PS=$(sort -k 6 -nr hapcut2.concatenated.phased.block.list | grep ${x} | cut -f3 | head -n 1)" > ${x}.hapcut.phased.largest.vcf
done

# concatenate vcf to one
bcftools concat chr*.hapcut.phased.largest.vcf > all.hapcut.phased.largest.vcf

# get the chr pos information
grep -v '^#' all.hapcut.phased.largest.vcf | cut -f1,2 > all.hapcut.phased.largest.pos.txt

# make bed file for the position
paste all.hapcut.phased.largest.pos.txt <(cut -f2 all.hapcut.phased.largest.pos.txt) > all.hapcut.phased.largest.pos.bed
