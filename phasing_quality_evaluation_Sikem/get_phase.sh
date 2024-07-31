#!/bin/bash

bcftools view -H -p whatshap_long_read_hic_phased.vcf.gz | cut -f1,2 > chr_pos.txt
bcftools view -H -p whatshap_long_read_hic_phased.vcf.gz | cut -f10 | cut -f1 -d':' | tr '|' '\t' > genotype.txt
paste chr_pos.txt genotype.txt > sikem_phase.txt
rm chr_pos.txt genotype.txt 


