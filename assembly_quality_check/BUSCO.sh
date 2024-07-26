#!/bin/bash

# 25.01.2024

# busco for sikem hifiasm hic and ul assemblies

database="/scratch/yutang/busco_database/embryophyta_odb10"

for asm in $(ls ../hifiasm_hic_ul_ctg)
do
    echo ${asm}
    mkdir busco_${asm}
    busco -i ../hifiasm_hic_ul_ctg/${asm} \
          -l ${database} \
          -m genome \
          -o busco_${asm} \
          -f \
          -c 60 \
          --offline
    echo 'Finished'
done
