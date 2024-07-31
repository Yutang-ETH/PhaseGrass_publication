# 05.04.2024

# compare phase obtained from LepMAP3 to phase obtained from PhaseGrass

# setwd
setwd("P:/Yutangchen/Sikem_hifi/validate_phase/genetic_linkage_map")

# import packages
library(stringr)
library(ggplot2)

# import the phased genetic linkage map
mymap1 <- read.table("lepmap_ont_phased_hifi_h1/scaffold_position_map_phased_with_ref_allele.txt", header = F, stringsAsFactors = F)
mymap2 <- read.table("lepmap_ont_phased_hifi_h2/scaffold_position_map_phased_with_ref_allele.txt", header = F, stringsAsFactors = F)

# name mymap
colnames(mymap1) <- c("chr", 'pos', 'FcM', 'McM', 'h1', 'h2', 'ref')
colnames(mymap2) <- c("chr", 'pos', 'FcM', 'McM', 'h1', 'h2', 'ref')

# add a color column, h1 ref
mycolor <- c()
for(i in unique(mymap1$chr)){
  
  xx <- mymap1[mymap1$chr == i, ]
  xx$color1 <- as.character(xx$h1 == xx$ref)
  xx$color2 <- as.character(xx$h2 == xx$ref)
  xx$color <- paste(xx$color1, xx$color2, sep="_")
  # there will be for comibinations in xx$color in theory
  # TRUE_FALSE, FALSE_TRUE, TRUE_TRUE, FALSE_FALSE
  if(length(which(xx$color == 'TRUE_FALSE')) >= length(which(xx$color == 'FALSE_TRUE'))){
    xx$color[xx$color == 'TRUE_FALSE'] <- 0
    xx$color[xx$color == 'FALSE_TRUE'] <- 1
    xx$color[xx$color == 'TRUE_TRUE'] <- 2
    xx$color[xx$color == 'FALSE_FALSE'] <- 2
  }else{
    xx$color[xx$color == 'TRUE_FALSE'] <- 1
    xx$color[xx$color == 'FALSE_TRUE'] <- 0
    xx$color[xx$color == 'TRUE_TRUE'] <- 2
    xx$color[xx$color == 'FALSE_FALSE'] <- 2
  }
  
  mycolor <- c(mycolor, xx$color)
  
}

mymap1$color <- mycolor

#h2 ref
mycolor <- c()
for(i in unique(mymap2$chr)){
  
  xx <- mymap2[mymap2$chr == i, ]
  xx$color1 <- as.character(xx$h1 == xx$ref)
  xx$color2 <- as.character(xx$h2 == xx$ref)
  xx$color <- paste(xx$color1, xx$color2, sep="_")
  # there will be for comibinations in xx$color in theory
  # TRUE_FALSE, FALSE_TRUE, TRUE_TRUE, FALSE_FALSE
  if(length(which(xx$color == 'TRUE_FALSE')) >= length(which(xx$color == 'FALSE_TRUE'))){
    xx$color[xx$color == 'TRUE_FALSE'] <- 0
    xx$color[xx$color == 'FALSE_TRUE'] <- 1
    xx$color[xx$color == 'TRUE_TRUE'] <- 2
    xx$color[xx$color == 'FALSE_FALSE'] <- 2
  }else{
    xx$color[xx$color == 'TRUE_FALSE'] <- 1
    xx$color[xx$color == 'FALSE_TRUE'] <- 0
    xx$color[xx$color == 'TRUE_TRUE'] <- 2
    xx$color[xx$color == 'FALSE_FALSE'] <- 2
  }
  
  mycolor <- c(mycolor, xx$color)
  
}

mymap2$color <- mycolor

# combine map 1 and 2
mymap_combine <- rbind.data.frame(mymap1, mymap2)

# get haplotype from chr column
mymap_combine$hap <- gsub('chr.*_', '', mymap_combine$chr)
mymap_combine$chr <- gsub('_h.*', '', mymap_combine$chr)

# capital chr and hap
mymap_combine$hap <- str_to_title(mymap_combine$hap)
mymap_combine$chr <- str_to_title(mymap_combine$chr)

mymap_combine$color <- as.factor(mymap_combine$color)

# get the number of 0 and 1 in each chr
# P0 <- c()
# P1 <- c()
# for(i in unique(mymap$chr)){
#   
#   xx <- mymap[mymap$chr == i, 8]
#   P0 <- c(P0, rep((length(xx[xx == 0])/length(xx))*100, length(xx)))
#   P1 <- c(P1, rep((length(xx[xx == 1])/length(xx))*100, length(xx)))
#   
# }
# 
# mymap$P0 <- P0
# mymap$P1 <- P1

# visualize phase between lepmap3 and phasegrass
ggplot(data = mymap_combine, aes(x = pos/10^6, y = FcM)) + 
  geom_point(data = mymap_combine, aes(x = pos/10^6, y = FcM, color = color)) +
  facet_wrap(hap~chr, nrow = 2, ncol = 7, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'right', legend.text = element_text(size=15), legend.title = element_text(size = 15)) + 
  labs(x = 'Base pair position in Sikem pseudo-chromosomes (Mb)', y = 'Position in genetic linkage maps (cM)', color = 'SNP') +
  scale_color_manual(values = c('0' =  'black', '1' = 'red', '2' = 'green'), labels = c('Same phase', 'Opposite phase', 'Genotyping error')) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text = element_text(size = 10)) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 20))
ggsave('LepMAP3_vs_PhaseGrass_phase.pdf', device = 'pdf', width = 20, height = 7.5)

# show the number of dots with the specific color
ggplot(data = mymap_combine, aes(x = color)) +
  ylim(0, 800) + 
  geom_bar(aes(fill = color), stat = 'count') +
  geom_text(stat = 'count', aes(label=..count..), vjust=-0.5, color="blue", size=3.5) +
  facet_wrap(hap~chr, nrow = 2, ncol = 7, scales = 'fixed') + 
  theme_bw() +
  theme(legend.position = 'right', legend.text = element_text(size=15), legend.title = element_text(size = 15)) + 
  labs(x = 'Comparison of phase between Lep-MAP3 and PhaseGrass', y = 'Number of SNPs', fill = 'SNP') +
  scale_fill_manual(values = c('0' =  'black', '1' = 'red', '2' = 'green'), labels = c('Same phase', 'Opposite phase', 'Genotyping error')) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text = element_text(size = 10), axis.text.x = element_blank()) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 20))
ggsave('LepMAP3_vs_PhaseGrass_number_of_SNPs.pdf', device = 'pdf', width = 20, height = 7.5)
