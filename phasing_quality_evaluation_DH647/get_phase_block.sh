#!/bin/bash

# get phased SNPs in phase blocks

bcftools view -H -p whatshap_long_read_hic_phased.vcf.gz | cut -f1,2 > chr_pos.txt
bcftools view -H -p whatshap_long_read_hic_phased.vcf.gz | cut -f10 | cut -f9 -d':' > phase_block.txt
paste chr_pos.txt phase_block.txt > DH647_phase_block.txt
rm chr_pos.txt phase_block.txt

