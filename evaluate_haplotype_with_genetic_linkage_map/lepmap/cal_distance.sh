#!/bin/bash

# calculate distance with given order

mypath=/home/yutachen/public/Yutangchen/Lepmap3_20240208/bin
post=/home/yutachen/public/Yutangchen/Sikem_hifi/validate_phase/genetic_linkage_map/all_post.gz

# import data
zcat ${post} | java -cp $mypath ParentCall2 data=pedigree.txt posteriorFile=- > p.call

# filter data
java -cp $mypath Filtering2 data=p.call dataTolerance=0.0001 removeNonInformative=1 > p_f.call

# make physical order
grep 'chr' p_f.call | cut -f1 > chr.txt
seq 1 $(wc -l chr.txt | cut -f1 -d " ") > physical_order.txt
paste physical_order.txt chr.txt > po.txt
rm physical_order.txt chr.txt

# calculate the distance
for i in $(seq 1 7)
do
    grep "chr$i" po.txt | cut -f1 > mypo_$i.txt 
    java -cp $mypath OrderMarkers2 evaluateOrder=mypo_$i.txt data=p_f.call improveOrder=0 numThreads=40 chromosome=$i informativeMask=1 >order_$i.txt 2>order_$i.err 
    echo "chr$i finished"
done

# cat chroms
cat order_*.txt > order_x.txt

# convert map to R/qtl format and get the phased genotypes
awk -f map2genotypes_phased.awk order_x.txt > genotypes_x.txt

# extract the row number of each SNP in the p.call file
cut -f 1,2 p_f.call | awk '(NR>=7)' > snps.txt

# match marker number to position
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt order_x.txt > order_x.mapped

# combine the maker position and the phased genotypes
paste <(grep -v "#" order_x.mapped | cut -f 1-4) <(grep -v 'CHR' genotypes_x.txt | cut -f5 | tr ' ' '\t') > scaffold_position_map_phased.txt

# make a bed file to get the reference allele
paste <(cut -f1 scaffold_position_map_phased.txt) <(awk '{ print $2 - 1 }' scaffold_position_map_phased.txt) <(cut -f2 scaffold_position_map_phased.txt) > scaffold_position_map_phased.bed

# use bedtools getfasta to get the reference allele
bedtools getfasta -fi ../genome/sikem_ont_phased_hifi_h1_final.fasta -bed scaffold_position_map_phased.bed | grep -v '>' > scaffold_position_map_reference_allele

# add ref allele to scaffold_position_map_phased.txt
paste scaffold_position_map_phased.txt scaffold_position_map_reference_allele > scaffold_position_map_phased_with_ref_allele.txt
