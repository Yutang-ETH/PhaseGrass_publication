#!/bin/bash

hifiasm -o hifiasm/redclover.asm -t 48 -u 0 --h1 HIC/SRR15433790_1.fastq,HIC/SRR15433791_1.fastq,HIC/SRR15433792_1.fastq --h2 HIC/SRR15433790_2.fastq,HIC/SRR15433791_2.fastq,HIC/SRR15433792_2.fastq HIFI/SRR15433789.fastq

awk '/^S/{print ">"$2;print $3}' hifiasm/redclover.asm.hic.hap1.p_ctg.gfa > hifiasm/redclover.hap1.ctg.fa
awk '/^S/{print ">"$2;print $3}' hifiasm/redclover.asm.hic.hap2.p_ctg.gfa > hifiasm/redclover.hap2.ctg.fa
awk '/^S/{print ">"$2;print $3}' hifiasm/redclover.asm.hic.p_utg.gfa > hifiasm/redclover.utg.fa

assembly-stats -t hifiasm/redclover.hap1.ctg.fa hifiasm/redclover.hap2.ctg.fa hifiasm/redclover.utg.fa > hifiasm/hifiasm.asm.stats 
