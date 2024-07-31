This folder contains scripts that were used to map rep and norep seqs to haplomes

align_unique_and_shared_rep is a snakemake pipeline maping norep and shared rep seqs between haplomes in one haplome to the other

bed_cigar_parse.py is a python script to parse cigar string in a bed file. This script is needed by align_unique_and_shared_rep

nonrep_rep_window_cov is a snakemake pipeline calculates the rep and norep seq coverage per 1 Mb 
