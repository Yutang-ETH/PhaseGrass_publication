# check phase, show density of SNPs in each 1 Mb window
# 14.08.2023

# setwd
setwd("P:/Yutangchen/DH647/DH647_nextpolish/phase_grass/whatshap_phase")

library(stringr)
library(ggplot2)

# import phase
phase <- read.table("DH647_phase.txt", header = F, stringsAsFactors = F)

# import chr length
fai <- read.table("../asm/kyuss.nextdenovo.juicer.fasta.chr.fa.fai", header = F, stringsAsFactors = F)


# output the switch count for all chr

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

# make an empty data structure to store data
phase_check <- NULL

for(j in 1:7){
  
  # get the phase for each chr
  chr_phase <- phase[phase$V1 == paste('chr', j, sep = ''), ]
  
  # chop chr to 1Mb window
  mywindow <- seqlast(0, fai[j, 2], 1000000)
  
  # calculate switch in every window
  myswitch <- c()
  num_snp <- c()
  for(i in 1:(length(mywindow) - 1)){
    
    chr_window <- chr_phase[chr_phase$V2 >= mywindow[i] & chr_phase$V2 < mywindow[i + 1], ]
    myswitch <- c(myswitch, round(nrow(chr_window[chr_window$V3 == 1, ])/nrow(chr_window), 2))
    num_snp <- c(num_snp, nrow(chr_window))
    print(i)
    
  }
  
  # get the middle position in each window
  myposition <- floor((mywindow[1:(length(mywindow)-1)] + mywindow[2:(length(mywindow))])/2)
  mystart <- mywindow[1:(length(mywindow)-1)] 
  myend <- mywindow[2:(length(mywindow))]
  
  # get a verctor of chr name with the same length as myswitch and myposition
  mychr <- chr_phase$V1[1:length(myposition)]
  
  # get the snp density in every 1 Mb
  mysnpdensity <- num_snp/1000000
  
  # cbind mychr, myposition, myswitch to form a data frame
  mydf <- cbind.data.frame(mychr, myposition, mystart, myend, myswitch, num_snp, mysnpdensity)
  
  # rename the colums in the mydf
  colnames(mydf) <- c("chr", "pos", "start", "end", "switch_count", "snp_number", "snp_density")
  
  # convert format of mydf
  mydf <- type.convert(mydf, as.is = T)
  
  # store all the information in phase_check
  phase_check <- rbind(phase_check, mydf)
  
  # print j
  print(paste('Job', j, 'has been completed!', sep = ' '))
  
}

# now visualize the check_phase table using dot plot
# make the chr capital
phase_check$chr <- str_to_title(phase_check$chr)
write.table(phase_check, 'phase_check.txt', quote = F, sep = '\t', col.names = T, row.names = F)

ggplot(phase_check, aes(x = pos, y = switch_count, color = switch_count)) +
  geom_point(size = 1) + 
  geom_line(aes(x = pos, y = snp_density*100), color = 'black') + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Number of heterozygous SNPs", breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1)*10000)) +
  scale_colour_gradient2(low = "blue", mid = 'gray', high = "red", midpoint = 0.5) +
  facet_wrap(~chr, ncol = 1) +
  theme_bw() +
  theme(legend.position = 'none') + 
  # theme(legend.position = c(0.85, 0.95), legend.direction = 'horizontal', legend.key.width = unit(1, 'cm'), legend.title = element_blank(), legend.key = element_rect(fill = NA)) +
  labs(x = 'Base pair position in Kyuss pesudo-chromosomes (Mb)', y = 'Ratio of non-reference alleles in the reference haplotype ') +
  scale_x_continuous(labels = seq(0, max(fai$V2), 5*10^7)/10^6, breaks = seq(0, max(fai$V2), 5*10^7)) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_line(linewidth = 0.2)) +
  theme(axis.title.y.right = element_text(angle = 90)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 30))

#ggsave('DH647_phase_grass_check_phase_with_snp_density_1Mb.png', width = 15, height = 12, units = 'in', dpi = 1500)
ggsave('DH647_phase_grass_check_phase_with_snp_density_1Mb.pdf', width = 15, height = 12)
