#!/bin/bash

# 05.06.2024

asm_h1=/scratch/yutang/sikem/chr_2_chr/asm/Sikem_h1_HiFi_ONT_final.fasta
asm_h2=/scratch/yutang/sikem/chr_2_chr/asm/Sikem_h2_HiFi_ONT_final.fasta
asm_dip=/scratch/yutang/sikem/accuracy_depth/sikem_dip/sikem_dip_final.fasta

short_read=/scratch/yutang/sikem/accuracy_depth/wgs_short

kat comp -t 30 -m 23 -h -H 50000000000 -o ./h1 "${short_read}/sikem_R?.noN.fastq.gz" ${asm_h1}
kat comp -t 30 -m 23 -h -H 50000000000 -o ./h2 "${short_read}/sikem_R?.noN.fastq.gz" ${asm_h2}
kat comp -t 30 -m 23 -h -H 50000000000 -o ./dip "${short_read}/sikem_R?.noN.fastq.gz" ${asm_dip}

kat plot spectra-cn -x 200 -o ./sikem_h1.png ./h1-main.mx
kat plot spectra-cn -x 200 -o ./sikem_h2.png ./h2-main.mx
kat plot spectra-cn -x 200 -o ./sikem_dip.png ./dip-main.mx

