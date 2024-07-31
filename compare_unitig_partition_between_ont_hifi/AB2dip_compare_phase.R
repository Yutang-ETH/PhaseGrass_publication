# set working directory
setwd("P:/Yutangchen/Sikem_hifi/validate_phase/")

# import packages
library(stringr)
library(ggplot2)

# import paf
mypaf_hifi <- read.table("AB2dip_filtered.paf", header = F, stringsAsFactors = F)
colnames(mypaf_hifi) <- c('query', 'qlen', 'qs', 'qe', 'strand', 'ref', 'rlen', 'rs', 're', 'acov', 'acovg', 'qv')

mypaf_ont <- read.table('ontAB2dip_filtered.paf', header = F, stringsAsFactors = F)
colnames(mypaf_ont) <- c('query', 'qlen', 'qs', 'qe', 'strand', 'ref', 'rlen', 'rs', 're', 'acov', 'acovg', 'qv')

# import utg assignment information
utgA_hifi <- read.table('../sort_grass/binned_seq/A.fa.fai', header = F, stringsAsFactors = F)
utgB_hifi <- read.table('../sort_grass/binned_seq/B.fa.fai', header = F, stringsAsFactors = F)

utgA_ont <- read.table('../sort_grass_with_ont_phase//binned_seq/A.fasta.fai', header = F, stringsAsFactors = F)
utgB_ont <- read.table('../sort_grass_with_ont_phase//binned_seq/B.fasta.fai', header = F, stringsAsFactors = F)

# some utg were split to more than one alignment, probably due to manual curation 
# I found that in paf, the longest alignment (without gaps) of a utg is always the first one (the first row of all rows of the utg)
mypaf_hifi <- mypaf_hifi[!duplicated(mypaf_hifi$query), ]
mypaf_ont <- mypaf_ont[!duplicated(mypaf_ont$query), ]

# split reference to chr and hap
mypaf_hifi$hap <- gsub('chr.*_', '', mypaf_hifi$ref)
mypaf_hifi$ref <- gsub('_h.*', '', mypaf_hifi$ref)

mypaf_ont$hap <- gsub('chr.*_', '', mypaf_ont$ref)
mypaf_ont$ref <- gsub('_h.*', '', mypaf_ont$ref)

# add utg assignment group to mypaf
mypaf_hifi$bin <- rep(1, nrow(mypaf_hifi))
mypaf_hifi$bin[which(mypaf_hifi$query %in% utgA_hifi$V1)] <- 'A'
mypaf_hifi$bin[which(mypaf_hifi$query %in% utgB_hifi$V1)] <- 'B'

mypaf_ont$bin <- rep(1, nrow(mypaf_ont))
mypaf_ont$bin[which(mypaf_ont$query %in% utgA_ont$V1)] <- 'A'
mypaf_ont$bin[which(mypaf_ont$query %in% utgB_ont$V1)] <- 'B'

# make chr and hap capital
mypaf_hifi$ref <- str_to_title(mypaf_hifi$ref)
mypaf_hifi$hap <- str_to_title(mypaf_hifi$hap)

mypaf_ont$ref <- str_to_title(mypaf_ont$ref)
mypaf_ont$hap <- str_to_title(mypaf_ont$hap)

# add data type
mypaf_hifi$tech <- rep('HiFi', nrow(mypaf_hifi))
mypaf_ont$tech <- rep('ONT', nrow(mypaf_ont))

# combine two df
mypaf_combine <- rbind.data.frame(mypaf_hifi, mypaf_ont)

# convert hap to factor
mypaf_combine$hap <- factor(mypaf_combine$hap, levels = c("H1", "H2"))

# convert data type to factor
mypaf_combine$tech <- factor(mypaf_combine$tech, levels = c('HiFi', 'ONT'))

# now visualize hap vs bin per chr
ggplot(data = mypaf_combine, aes(x = hap, y = rs)) +
  geom_point(data = mypaf_combine, aes(color = tech), size = 0.5, shape = 19) +
  facet_grid(cols = vars(ref), rows = vars(bin, tech), scales = 'free') +
  ylab('Base pair position in Sikem pseudo-chromosomes (Mb)') +
  xlab('Sikem haplomes') +
  scale_y_continuous(labels = seq(0, max(mypaf_combine$rlen), 10*10^7)/10^6,
                     breaks = seq(0, max(mypaf_combine$rlen), 10*10^7)) +
  scale_color_manual(values=c("red", "blue"), breaks = c("HiFi", "ONT")) + 
  labs(color = 'Phasing technology') +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(strip.text.x = element_text(face = "plain", size = 15),
        strip.text.y = element_text(face = "plain", size = 15),
        strip.background = element_blank(),
        legend.position="top", 
        legend.text = element_text(size=15),
        legend.title = element_text(size = 15),
        legend.key = element_rect(fill = 'white'),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.5, linetype = 'solid', colour = "gray90"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(color="black", size=20, face="plain"),
        axis.title.y = element_text(color="black", size=20, face="plain"))
# ggsave('phase_block_hifi_vs_ont.png', device = 'png', width = 15, height = 3, dpi = 480)
ggsave('validate_phase_sikem_hifihic_vs_hifionthic_new.pdf', device = 'pdf', width = 10, height = 12)
