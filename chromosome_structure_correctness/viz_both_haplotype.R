# visualize h1 to kyuss and h2 to kyuss together

# set working directory
setwd("P:/Yutangchen/DH647/DH647_nextpolish/chr_align_chr")

library(ggplot2)
library(stringr)

# import paf
h1paf <- read.table("DH647_h1_vs_kyuss/q2r_filtered.paf", header = F, stringsAsFactors = F)
h2paf <- read.table("DH647_h2_vs_kyuss/q2r_filtered.paf", header = F, stringsAsFactors = F)

# rename mypaf
colnames(h1paf) <- c("qchr", 'qlen', 'qs', 'qe', 'strand', 'rchr', 'rlen', 'rs', 're', 'Mmatch', 'Mall', 'MQ')
colnames(h2paf) <- c("qchr", 'qlen', 'qs', 'qe', 'strand', 'rchr', 'rlen', 'rs', 're', 'Mmatch', 'Mall', 'MQ')

# remove lines with unplaced
h1paf <- h1paf[grepl('chr', h1paf$qchr) & grepl('chr', h1paf$rchr), ]
h2paf <- h2paf[grepl('chr', h2paf$qchr) & grepl('chr', h2paf$rchr), ]

# split qchr to new qchr without _. and qid without chr._
h1paf$qid <- as.numeric(gsub('-.*', '', (gsub('chr._', '', h1paf$qchr))))
h1paf$qchr <- gsub('_.*', '', h1paf$qchr)
h2paf$qid <- as.numeric(gsub('-.*', '', (gsub('chr._', '', h2paf$qchr))))
h2paf$qchr <- gsub('_.*', '', h2paf$qchr)

# only selecte alignment between the same chrs
h1paf <- h1paf[h1paf$qchr == h1paf$rchr, ]
h2paf <- h2paf[h2paf$qchr == h2paf$rchr, ]

# filtering based on alignment length, and MQ >= 50
h1paf <- h1paf[h1paf$Mmatch >= 2000 & h1paf$MQ >= 60, ]
h2paf <- h2paf[h2paf$Mmatch >= 2000 & h2paf$MQ >= 60, ]

# add haplotype to table
h1paf$hap <- rep('H1', nrow(h1paf))
h2paf$hap <- rep('H2', nrow(h2paf))

# make chr capital
h1paf$qchr <- str_to_title(h1paf$qchr)
h1paf$rchr <- str_to_title(h1paf$rchr)
h2paf$qchr <- str_to_title(h2paf$qchr)
h2paf$rchr <- str_to_title(h2paf$rchr)

# combin h1 and h2 paf
mypaf <- rbind.data.frame(h1paf, h2paf)

# visualize the paf, single by single
ggplot(data = mypaf, aes(x = rs/1000000, y = qid/1000000)) +
  geom_point(size = 0.25, shape = 19, color = 'black') +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        panel.grid = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(), 
        strip.text = element_text(size = 15)) +
  labs(x = "Base pair position in Kyuss pseudo-chromosomes (Mb)", y = "Base pair position in DH647 pseudo-chromosomes (Mb)", color = "Mapping\nquality") +
  facet_wrap(hap ~ qchr,
             nrow = 2,
             ncol = 7,
             scales = 'free')

ggsave('DH647_h1_h2_vs_kyuss_ont.pdf', width = 20, height = 7.5)
