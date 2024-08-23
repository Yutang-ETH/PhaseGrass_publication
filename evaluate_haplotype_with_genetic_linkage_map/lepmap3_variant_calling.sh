#!/bin/bash

# try lep_map3 variant calling pipeline

# path to the Lempam3 package
mypath=/home/yutachen/public/Yutangchen/Lepmap3_20240208/bin

# generate new_sample_order
# map is the folder where the *.bam files are stored
cat <(ls -1 map | grep 'bam$' | xargs basename -s '.bam' | grep -v 'lmGbsJP') <(ls -1 map | grep 'bam$' | xargs basename -s '.bam' | grep 'lmGbsJP') > new_sample_order

# prepare the sorted_bams file
sed 's/^/map\//' new_sample_order | sed 's/$/.bam/' | tr '\n' '\t' > sorted_bams

# prepare the mapping.txt file required by Pileup2Likelihoods
cat new_sample_order > mapping.txt

# run variant calling in parallel
for i in chr1_h2 chr2_h2 chr3_h2 chr4_h2 chr5_h2 chr6_h2 chr7_h2
do
    echo "samtools mpileup -r \"$i\" -q 10 -Q 10 -s \$(cat sorted_bams) | java -cp $mypath Pileup2Likelihoods > \"$i\".post"
done >SNP_calling.txt

# run samtools mpileup in parallel
parallel --jobs 48 < SNP_calling.txt

# cat all the post files
cat *.post | awk '(NR==1 || ($1!="CHR"))' | gzip > all_post.gz
