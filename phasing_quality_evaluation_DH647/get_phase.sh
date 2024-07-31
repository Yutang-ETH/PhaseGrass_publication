#!/bin/bash

bcftools view -H -p whatshap_long_read_hic_phased.vcf.gz | cut -f1,2 > chr_pos.txt
bcftools view -H -p whatshap_long_read_hic_phased.vcf.gz | cut -f10 | cut -f1 -d':' | cut -f1 -d'|' > genotype.txt
paste chr_pos.txt genotype.txt > DH647_phase.txt
rm chr_pos.txt genotype.txt 


