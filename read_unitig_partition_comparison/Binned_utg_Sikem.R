# set working directory
setwd("P:/Yutangchen/pan_asm_evaluation")

# import packages
library(stringr)
library(ggplot2)
library(ggh4x)

# import data
mybin <- read.table('Binned_utg_Sikem.txt', header = T, stringsAsFactors = F, sep = '\t')

mybin$method[mybin$method == 'K-mer'] <- 'SortGrass'
mybin$method[mybin$method == 'Hifiasm'] <- 'Hifiasm (UL + Hi-C)'
mybin$data[mybin$data == '0.55'] <- '-s 0.55'
mybin$data[mybin$data == '0.20'] <- '-s 0.20'
mybin$data[mybin$data == 'PacBio HiFi'] <- 'PacBio HiFi + Hi-C'
mybin$data[mybin$data == 'ONT'] <- 'ONT + Hi-C'

mybin$data <- factor(mybin$data, levels = c('-s 0.55', '-s 0.20', 'PacBio HiFi + Hi-C', 'ONT + Hi-C'))

ggplot(mybin, aes(x = genotype, y = sequence, fill = haplotype)) + 
  geom_bar(stat = 'identity', color="black", width = 0.5, position = position_dodge()) +
  facet_nested(cols = vars(genotype, method, data), scales = 'free', space = 'free',
               nest_line = element_line(linetype = 1),
               strip = strip_nested(clip = 'off', text_x = elem_list_text(face = c("plain"), size = c(30, 20, 20)), by_layer_x = TRUE)) +
  scale_fill_manual(values = c('steelblue2', 'maroon2', 'gray', 'gray90'), labels = c('H1', 'H2', 'None', 'Unitig')) +
  scale_y_continuous(labels = c("0","0.5", "1.0", "1.5","2.0", "2.5", '3.0', '3.5', '4.0', '4.5', '5.0'), breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)) + 
  xlab('') +
  ylab('Total length (Gb)') + 
  theme(
    strip.background = element_blank(),
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    legend.key = element_rect(fill = 'white'),
    axis.title.y = element_text(size = 25),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.5, "lines"), 
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5, linetype = "solid"))

ggsave('Sikem_utg_binning.pdf', device = 'pdf', width = 10, height = 8)

