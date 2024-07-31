#!/usr/bin/env python3

# I want to parse cigar string in the bed file
# I want to know the following:
# total length of matched bases
# total length of SNVs
# total length of consecutive SNVs, namely multiple nucleotide substitution
# total length of small insertions >=1 & <50
# total length of large insertions >50
# total length of smal deletions >=1 & <50
# total length of large deletions >50
# total length of soft clipping

# import library
import re
import sys

# only one positional argument, the input bed
myinput = sys.argv[1]

# standard output, print results to the screen
with open(myinput, 'r') as f:
    for line in f:
        
        # get cigar string
        cigar = str(line.split('\t')[6])
        
        # total length of matched bases
        baseM = sum([int(i.split('=')[0]) for i in re.findall('\d+=', cigar)])
        
        # baseVar all unmatched bases
        baseVar = [int(i.split('X')[0]) for i in re.findall('\d+X', cigar)]
        
        # baseVar can be divided to SNV and MNS
        # total length of SNV
        baseSNV = sum([var for var in baseVar if var == 1])
        # total length of multiple nucleotide substitution
        baseMNS = sum([var for var in baseVar if var > 1])
        
        # baseINS
        baseINS = [int(i.split('I')[0]) for i in re.findall('\d+I', cigar)]
        # baseSINS, small insertions with length >=1 and <50
        baseSINS = sum([var for var in baseINS if var >= 1 & var < 50])
        # baseLINS, large insertions with length >=50
        baseLINS = sum([var for var in baseINS if var >= 50])

        # baseDEL
        baseDEL = [int(i.split('D')[0]) for i in re.findall('\d+D', cigar)]
        # baseSDEL, small deletions with length >=1 and <50
        baseSDEL = sum([var for var in baseDEL if var >= 1 & var < 50])
        # baseLDEL, large insertions with length >=50
        baseLDEL = sum([var for var in baseDEL if var >= 50])

        # total length of soft clipping
        baseSC = sum([int(i.split('S')[0]) for i in re.findall('\d+S', cigar)])

        # the totla length of alignment based on the cigar string, inlcuding mathced, unmatched, all INDELs and softclipped 
        alignLen = baseM + baseSNV + baseMNS + baseSINS + baseLINS + baseSDEL + baseLDEL + baseSC

        # get the query chr, start, end and length
        query = str(line.split('\t')[3])
        qchr = query.split('_')[0]
        qstart = query.split('_')[1].split('-')[0]
        qend = query.split('_')[1].split('-')[1]
        qlen = int(qend) - int(qstart) + 1

        # the header of print out
        # rid, rs, re, qid, mapping_quality, strand, cigar, matched_bases, SNV, MNS, SINS, LINS, SDEL, LDEL, SC, alignlen
        print(line.rstrip() + '\t' + 
        	  str(qchr)	+ '\t' +
        	  str(qstart) + '\t' +
        	  str(qend) + '\t' +
        	  str(qlen) + '\t' +
        	  str(baseM) + '\t' + 
        	  str(baseSNV) + '\t' +
        	  str(baseMNS) + '\t' + 
        	  str(baseSINS) + '\t' +
        	  str(baseLINS) + '\t' +
        	  str(baseSDEL) + '\t' +
        	  str(baseLDEL) + '\t' +
        	  str(baseSC) + '\t' + 
        	  str(alignLen))
