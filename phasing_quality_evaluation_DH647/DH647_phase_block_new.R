# set working directory
setwd("P:/Yutangchen/DH647/DH647_nextpolish/phase_grass/whatshap_stats")

# import packages
library(stringr)
library(ggplot2)
library(ggh4x)

# import data
DH647_phase <- read.table('phased.list', header = F, stringsAsFactors = F, sep = '\t')
DH647_fai <- read.table('../asm/kyuss.nextdenovo.juicer.fasta.chr.fa.fai', header = F, stringsAsFactors = F)
DH647_snp <- read.table('../whatshap_phase/DH647_phase_block.txt', header = F, stringsAsFactors = F)

# name the df
colnames(DH647_phase) <- c('genotype', 'chr', 'block', 'start', 'end', 'variants')
colnames(DH647_snp) <- c('chr', 'pos', 'block')

# make chr names capital, namely, Chr1
DH647_phase$chr <- str_to_title(DH647_phase$chr)

# get the largest and small phase block for each genotype
DH647_large <- NULL
DH647_small <- NULL

for(x in unique(DH647_phase$genotype)){
  
  # get the genotype df
  xx <- DH647_phase[DH647_phase$genotype == x, ]
  
  # get the large and small block in each genotype
  # yy_large <- NULL
  # yy_small <- NULL
  for(y in paste('Chr', 1:7, sep = '')){
    
    yy <- xx[xx$chr == y, ]
    yy <- yy[order(yy$variants, decreasing = T), ]
    
    DH647_large <- type.convert(rbind.data.frame(DH647_large, yy[1, ]), as.is = T)
    DH647_small <- type.convert(rbind.data.frame(DH647_small, yy[-1, ]), as.is = T)
    
  }
  
  
}

# add ymin and ymax 
DH647_large$ystart <- rep(1, nrow(DH647_large))
DH647_large$yend <- rep(1.2, nrow(DH647_large))

# add y position for DH647_small
DH647_small$ypos <- rep(2, nrow(DH647_small))
write.table(DH647_small, 'DH647_small_blocks.txt', quote = F, sep = '\t', row.names = F, col.names = T)

#########################################################################################################################################################
# add length to DH647_large
DH647_large$length <- DH647_fai$V2

# select SNPs from the largest block
DH647_snp_large <- DH647_snp[DH647_snp$block %in% DH647_large$block, ]

# function to have the last value, https://stackoverflow.com/questions/28419281/missing-last-sequence-in-seq-in-r 
seqlast <- function (from, to, by) 
{
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

# calculate SNP density per 1 Mb for the largest block of each chr
block_snp <- NULL
for(j in 1:7){
  
  # get the snp for the largest block of each chr
  chr_snp <- DH647_snp_large[DH647_snp_large$chr == paste('chr', j, sep = ''), ]
  
  # chop chr to 1Mb window
  mywindow <- seqlast(min(chr_snp$pos), max(chr_snp$pos), 1000000)

  
  # calculate number of snp in every window
  num_snp <- c()
  for(i in 1:(length(mywindow) - 1)){
    
    chr_window <- chr_snp[chr_snp$pos >= mywindow[i] & chr_snp$pos < mywindow[i + 1], ]
    num_snp <- c(num_snp, nrow(chr_window))
    print(i)
    
  }
  
  # get the start and end position of each window
  mystart <- mywindow[1:(length(mywindow)-1)]
  myend <- mywindow[2:length(mywindow)]
  
  # get a verctor of chr name with the same length as mystart and myend
  mychr <- rep(paste('chr', j, sep = ''), length(mystart))
  
  # cbind mychr, myposition, myswitch to form a data frame
  mydf <- cbind.data.frame(mychr, mystart, myend, num_snp)
  
  # rename the colums in the mydf
  colnames(mydf) <- c("chr", "start", "end", "nsnp")
  
  # convert format of mydf
  mydf <- type.convert(mydf, as.is = T)
  
  # store all the information in block_snp
  block_snp <- rbind(block_snp, mydf)
  
  # print j
  print(paste('Job', j, 'has been completed!', sep = ' '))
  
}

# make chr capital in block_snp
block_snp$chr <- str_to_title(block_snp$chr)
write.table(block_snp, 'DH647_largest_blocks.txt', quote = F, sep = '\t', row.names = F, col.names = T)

# # combine block_snp with DH647_small
# DH647_small <- DH647_small[, c(2, 4, 5, 6)]
# colnames(DH647_small) <- colnames(block_snp)
# DH647_small$block <- rep('Small', nrow(DH647_small))
# 
# # add block to block_snp
# block_snp$block <- rep("Largest", nrow(block_snp))
# 
# block_combine <- rbind.data.frame(DH647_small, block_snp)

# visualize the block
ggplot(data = DH647_large, aes(x = length, y = yend)) +
  geom_rect(data = block_snp, aes(xmin = start, ymin = 1, xmax = end, ymax = 1.2, fill = nsnp), inherit.aes = F) + 
  geom_point(data = DH647_small, 
             aes(x = block, y = 1.7, size = log(variants, 5)), 
             alpha = 0.3, color = 'black', fill = 'blue', shape = 21, stroke = 1, 
             position = position_jitter(width = 1000000, height = 0.2)) +
  # facet_grid(cols = vars(chr), scales = 'free')+ 
  facet_grid(rows = vars(chr), scales = 'free') +
  ylim(c(0.8, 2)) +
  xlab('Base pair position in Kyuss psuedo-chromosomes (Mb)') +
  ylab('') +
  scale_fill_gradient(low = 'white', high = 'red', limits = c(0, 8000), 
                      guide = guide_legend(title = 'Number of SNPs in chromosome-level phase blocks',
                                           theme(legend.position = 'top',
                                                 legend.direction = 'horizontal',
                                                 legend.key.width = unit(1, 'cm'),
                                                 legend.title = element_text(size = 15, vjust = 0.9),
                                                 legend.text = element_text(size = 15),
                                                 legend.text.position = 'bottom',  
                                                 legend.key = element_rect(fill = NA, color='black')))) +
  scale_size('Number of SNPs in samll phase blocks', breaks = c(1, 2, 3), labels = c('5', '25', '>100'), range = c(0, 5)) +
  scale_x_continuous(labels = seq(0, max(DH647_large$length), 5*10^7)/10^6, breaks = seq(0, max(DH647_large$length), 5*10^7)) +
  theme(strip.text.y = element_text(face = "plain", size = 15),
        strip.background = element_blank(),
        legend.position="top", 
        legend.text = element_text(size=15),
        legend.title = element_text(size = 15),
        legend.key = element_rect(fill = NA, colour = NA),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray95"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_text(color="black", size=20, face="plain"),
        axis.title.y = element_text(color="black", size=15, face="plain"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
# ggsave('DH647_phase_with_SNP_density.pdf', device = 'pdf', width = 20, height = 3)
ggsave('DH647_phase_grass_check_phase_with_snp_density_1Mb_vertical.pdf', device = 'pdf', width = 15, height = 12)
