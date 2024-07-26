# set working directory
setwd("P:/Yutangchen/pan_asm_evaluation/phasegrass_manual")

# import packages
library(stringr)
library(ggplot2)
library(ggh4x)

# import data
# long read alignment depth
long_hifi <- read.table("depth/sikem_hifi.regions.hifi.ont.hic.bed", header = F, stringsAsFactors = F)
long_ONT <- read.table("depth/DH647_ont.regions.bed", header = F, stringsAsFactors = F)

myformat <- function(long_dep = df){
  
  # add col names
  colnames(long_dep) <- c('chr', 'start', 'end', 'depth')
  
  # split chr col to chr and haplotype
  long_dep$hap <- gsub('chr.*_', '', long_dep$chr)
  long_dep$chr <- gsub('_h.*', '', long_dep$chr)
  
  # only select chr lines
  long_dep <- long_dep[grep('chr', long_dep$chr), ]
  
  # make the chr and hap capital
  long_dep$chr <- str_to_title(long_dep$chr)
  long_dep$hap <- str_to_title(long_dep$hap)
  
  return(long_dep)
  
}

long_ONT <- myformat(long_dep = long_ONT)
long_hifi <- myformat(long_dep = long_hifi)

# visualize the depth hifi
ggplot(long_hifi, aes(x = ((start + end)/2)/10^6, y = depth)) +
  geom_area(color = 'black', fill = 'gray', alpha = 0.5) + 
  facet_grid(row = vars(chr), col = vars(hap), scales = 'free') +
  xlab('Base pair position in Sikem pseudo-chromosomes (Mb)') +
  ylab('PacBio HiFi alignment depth') + 
  ylim(0, 50) + 
  theme(strip.text.x = element_text(face = "plain", size = 25),
        strip.text.y = element_text(face = "plain", size = 25),
        strip.background = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(color="black", size=25, face="plain"),
        axis.title.y = element_text(color="black", size=25, face="plain"))
ggsave('sikem_hifi_alignment_depth_hifi_ont_hic.pdf', width = 15, height = 10)

# visualize the depth ont
ggplot(long_ONT, aes(x = ((start + end)/2)/10^6, y = depth)) +
  geom_area(color = 'black', fill = 'gray', alpha = 0.5) + 
  facet_grid(row = vars(chr), col = vars(hap), scales = 'free') +
  xlab('Base pair position in DH647 pseudo-chromosomes (Mb)') +
  ylab('ONT alignment depth') + 
  ylim(0, 70) + 
  theme(strip.text.x = element_text(face = "plain", size = 25),
        strip.text.y = element_text(face = "plain", size = 25),
        strip.background = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(color="black", size=25, face="plain"),
        axis.title.y = element_text(color="black", size=25, face="plain"))
ggsave('DH647_ont_alignment_depth.pdf', width = 15, height = 10)
