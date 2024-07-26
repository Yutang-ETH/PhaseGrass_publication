# set working directory
setwd("P:/Yutangchen/pan_asm_evaluation")

# import packages
library(stringr)
library(ggplot2)
library(ggh4x)

# import data
mybin <- read.table('Binned_ont_DH647_Sikem.txt', header = T, stringsAsFactors = F, sep = '\t')

# replace Alignment with whatshap, K-mer with SortGrass
mybin$method[mybin$method == 'Alignment'] <- 'WhatsHap'
mybin$method[mybin$method == 'K-mer'] <- 'SortGrass'

mybin$method <- factor(mybin$method, levels = c("WhatsHap", 'SortGrass'))

# calculate the ratio of seq for each bin
seqr <- c()
for(i in 1:nrow(mybin)){
  
  seqr <- c(seqr, mybin$sequence[i]/sum(mybin[mybin$genotype == mybin$genotype[i] & mybin$method == mybin$method[i], 3]))
  
}
mybin$ratio <- seqr

# plot DH647 only
mybin_DH647 <- mybin[mybin$genotype == 'DH647', ]
mybin_Sikem <- mybin[mybin$genotype == 'Sikem', ]

ggplot(mybin_Sikem, aes(x = genotype, y = ratio*100, fill = haplotype)) + 
  geom_bar(stat = 'identity', color="black", width = 0.5, position = position_dodge()) +
  facet_nested(cols = vars(genotype, method), scales = 'free', space = 'free',
               nest_line = element_line(linetype = 1),
               strip = strip_nested(clip = 'off', text_x = elem_list_text(face = c("plain"), size = c(20, 15)), by_layer_x = TRUE)) +
  scale_fill_manual(values = c('steelblue2', 'maroon2', 'gray'), labels = c('H1', 'H2', 'None')) +
  scale_y_continuous(labels = c("0","10", "20", "30","40", "50"), breaks = c(0,10,20,30,40,50)) + 
  xlab('') +
  ylab('% of binned ONT reads') + 
  theme(
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size=15),
    legend.key = element_rect(fill = 'white'),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.5, "lines"), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5, linetype = "solid"))

# ggsave('DH647_Sikem_ONT_read_binning.pdf', device = 'pdf', width = 10, height = 5)
ggsave('DH647_ONT_read_binning.pdf', device = 'pdf', width = 5, height = 5)
ggsave('Sikem_ONT_read_binning.pdf', device = 'pdf', width = 5, height = 5)
