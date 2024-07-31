# check phase, show density of SNPs in each 1 Mb window
# 14.08.2023

# setwd
setwd("P:/Yutangchen/DH647/DH647_nextpolish/phase_grass/whatshap_phase")

library(stringr)
library(ggplot2)

# import phase
phase_kyuss <- read.table("DH647_phase.txt", header = F, stringsAsFactors = F)
phase_h1 <- read.table('../../validate_phase/h1_chop/h1.all.snps', header = F, stringsAsFactors = F)
phase_h2 <- read.table('../../validate_phase/h2_chop/h2.all.snps', header = F, stringsAsFactors = F)

# import chr length
fai <- read.table("../asm/kyuss.nextdenovo.juicer.fasta.chr.fa.fai", header = F, stringsAsFactors = F)

# name phase_kyuss
colnames(phase_kyuss) <- c('chr', 'pos', 'allele')

# to not be confused with the alleles in phase_kyuss, delete that column
phase_kyuss <- phase_kyuss[, 1:2]

# select chr, pos, allel from h1 and h2
h1 <- phase_h1[, c(11, 1, 2, 3)]
colnames(h1) <- c('chr', 'pos', 'a1', 'a2')

h2 <- phase_h2[, c(11, 1, 2, 3)]
colnames(h2) <- c('chr', 'pos', 'a1', 'a2')

# merge kyuss, h1, h2
h1_merge <- NULL
h2_merge <- NULL
for(x in paste('chr', 1:7, sep = '')){
  
  xx <- phase_kyuss[phase_kyuss$chr == x, ]
  yy <- h1[h1$chr == x, ]
  zz <- h2[h2$chr == x, ]
  
  h1_merge <- type.convert(rbind.data.frame(h1_merge, merge(xx, yy, by = 'pos', all.x = T, all.y = F)), as.is = T)
  h2_merge <- type.convert(rbind.data.frame(h2_merge, merge(xx, zz, by = 'pos', all.x = T, all.y = F)), as.is = T)
  
}

#################################################################################################################
# function to have the last value, https://stackoverflow.com/questions/28419281/missing-last-sequence-in-seq-in-r 
# use seqlast instead of seq
seqlast <- function (from, to, by) 
{
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

h1_merge_n <- na.omit(h1_merge)
h1_merge_n <- h1
# calculate the SNP density in every 1 Mb
# make an empty data structure to store data
h1_check <- NULL

for(j in 1:7){
  
  # get the phase for each chr
  # chr_snp <- h1_merge_n[h1_merge_n$chr.x == paste('chr', j, sep = ''), ]
  chr_snp <- h1_merge_n[h1_merge_n$chr == paste('chr', j, sep = ''), ]
  
  # chop chr to 1Mb window
  mywindow <- seqlast(0, fai[j, 2], 1000000)
  
  # calculate switch in every window
  num_snp <- c()
  for(i in 1:(length(mywindow) - 1)){
    
    chr_window <- chr_snp[chr_snp$pos >= mywindow[i] & chr_snp$pos < mywindow[i + 1], ]
    num_snp <- c(num_snp, nrow(chr_window))
    
  }
  
  # get the middle position in each window
  myposition <- floor((mywindow[1:(length(mywindow)-1)] + mywindow[2:(length(mywindow))])/2)
  mystart <- mywindow[1:(length(mywindow)-1)] 
  myend <- mywindow[2:(length(mywindow))]
  
  # get a verctor of chr name with the same length as myswitch and myposition
  mychr <- rep(paste('chr', j, sep = ''), length(myposition)) # chr_snp$chr.x[1:length(myposition)]
  
  # get the snp density in every 1 Mb
  mysnpdensity <- num_snp
  
  # cbind mychr, myposition, myswitch to form a data frame
  mydf <- cbind.data.frame(mychr, myposition, mystart, myend, mysnpdensity)
  
  # rename the colums in the mydf
  colnames(mydf) <- c("chr", "pos", 'start', 'end', "snp_density")
  
  # convert format of mydf
  mydf <- type.convert(mydf, as.is = T)
  
  # store all the information in phase_check
  h1_check <- rbind(h1_check, mydf)
  
  # print j
  print(paste('Job', j, 'has been completed!', sep = ' '))
  
}

h1_check$hap <- rep('H1', nrow(h1_check))


h2_merge_n <- na.omit(h2_merge)
h2_merge_n <- h2
# calculate the SNP density in every 1 Mb
# make an empty data structure to store data
h2_check <- NULL

for(j in 1:7){
  
  # get the phase for each chr
  # chr_snp <- h2_merge_n[h2_merge_n$chr.x == paste('chr', j, sep = ''), ]
  chr_snp <- h2_merge_n[h2_merge_n$chr == paste('chr', j, sep = ''), ]
  
  # chop chr to 1Mb window
  mywindow <- seqlast(0, fai[j, 2], 1000000)
  
  # calculate switch in every window
  num_snp <- c()
  for(i in 1:(length(mywindow) - 1)){
    
    chr_window <- chr_snp[chr_snp$pos >= mywindow[i] & chr_snp$pos < mywindow[i + 1], ]
    num_snp <- c(num_snp, nrow(chr_window))
    
  }
  
  # get the middle position in each window
  myposition <- floor((mywindow[1:(length(mywindow)-1)] + mywindow[2:(length(mywindow))])/2)
  mystart <- mywindow[1:(length(mywindow)-1)] 
  myend <- mywindow[2:(length(mywindow))]
  
  # get a verctor of chr name with the same length as myswitch and myposition
  mychr <- rep(paste('chr', j, sep = ''), length(myposition)) # chr_snp$chr.x[1:length(myposition)]
  
  # get the snp density in every 1 Mb
  mysnpdensity <- num_snp
  
  # cbind mychr, myposition, myswitch to form a data frame
  mydf <- cbind.data.frame(mychr, myposition, mystart, myend, mysnpdensity)
  
  # rename the colums in the mydf
  colnames(mydf) <- c("chr", "pos", 'start', 'end', "snp_density")
  
  # convert format of mydf
  mydf <- type.convert(mydf, as.is = T)
  
  # store all the information in phase_check
  h2_check <- rbind(h2_check, mydf)
  
  # print j
  print(paste('Job', j, 'has been completed!', sep = ' '))
  
}

h2_check$hap <- rep('H2', nrow(h2_check))


# combine h1 and h2
combine_check <- rbind.data.frame(h1_check, h2_check)
combine_check$chr <- str_to_title(combine_check$chr)
combine_check$hap <- factor(combine_check$hap, levels = c('H2', 'H1'))
# this is SNPs merged with SNPs called from short reads
# write.table(combine_check, 'SNPs_hap_vs_kyuss_merged.txt', quote = F, sep = '\t', row.names = F, col.names = T)
# this is SNPs from NUCMER + show-snps
write.table(combine_check, 'SNPs_hap_vs_kyuss_nucmer.txt', quote = F, sep = '\t', row.names = F, col.names = T)

ggplot(data = combine_check, aes(x = pos, y = hap, fill = snp_density)) + 
  geom_tile() + 
  facet_wrap(~chr, ncol = 1) +
  scale_fill_gradient(low = 'white', high = 'red4') +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = 'Base pair position in Kyuss pseudo-chromosomes (Mb)', y = 'Number of SNPs found between DH647 haplotypes and Kyuss', fill = 'SNP density') +
  scale_x_continuous(labels = seq(0, max(fai$V2), 5*10^7)/10^6,
                     breaks = seq(0, max(fai$V2), 5*10^7)) +
  theme(legend.position = c(0.85, 0.95),
        legend.direction = 'horizontal',
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size = 20, vjust = 1),
        legend.text = element_text(size = 15),
        legend.text.position = 'bottom',  
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  # guides(fill = guide_legend(title.position = 'top', lable.position = "bottom")) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(axis.title.y.right = element_text(angle = 90)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 30))
# ggsave('DH647_phase_grass_check_phase_chr_alignment.png', width = 15, height = 12, units = 'in', dpi = 1500)
# ggsave('DH647_phase_grass_check_phase_chr_alignment_heatmap.pdf', width = 15, height = 12)
ggsave('DH647_phase_grass_check_phase_chr_alignment_heatmap_nucmer.pdf', width = 15, height = 12)
