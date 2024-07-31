# 14.04.2024

# PhaseGrass ONT vs PhaseGrass hifi
setwd("P:/Yutangchen/Sikem_hifi/plot_phase")

library(stringr)

# import data
phase_ont <- read.table('sikem_phase_ont.txt', header = F, stringsAsFactors = F)
phase_hifi <- read.table("sikem_phase_hifi.txt", header = F, stringsAsFactors = F)

# import chr length
fai <- read.table("sikem.hifiasm.juicer.fasta.chr.fa.fai", header = F, stringsAsFactors = F)

# add col name
colnames(phase_ont) <- c("chr", 'pos', 'h1', 'h2')
colnames(phase_hifi) <- c('chr', 'pos', 'h1', 'h2')

# add snp tag
phase_ont$snp <- paste(phase_ont$chr, phase_ont$pos, sep = '_')
phase_hifi$snp <- paste(phase_hifi$chr, phase_hifi$pos, sep = '_')

# select common snp
common_snp <- intersect(phase_ont$snp, phase_hifi$snp)
phase_ont <- phase_ont[phase_ont$snp %in% common_snp, ]
phase_hifi <- phase_hifi[phase_hifi$snp %in% common_snp, ]

phase_hifi$snp <- phase_ont$h1
colnames(phase_hifi)[5] <- 'ont'

# compare both hifi h1 and h2 to ONT 
# the thing is h1 in chr1 might be the same to the ont hap, 
# but in chr2, h1 might be the opposite phase to the ont hap
# so here, I want to find out in each chr, which haplotype is same to the ont hap
# the number of phase switches should be low, therefore, based on this I could tell whether h1 or h2 is the same to the ref haplotype
hap_hifi <- c()
for(i in unique(phase_hifi$chr)){
  
  xx <- phase_hifi[phase_hifi$chr == i, ]
  
  # there will be two combinations in xx$combination in theory
  # TRUE_FALSE, FALSE_TRUE
  xx$combination <- paste(as.character(xx$h1 == xx$ont), as.character(xx$h2 == xx$ont), sep="_")
  
  # if there are more TRUE_FALSE, then xx$h1 should be the one same to the ont phase
  # otherwise, xx$h2 should be the one same to the ont phase
  if(length(which(xx$combination == 'TRUE_FALSE')) >= length(which(xx$combination == 'FALSE_TRUE'))){
    print(length(which(xx$combination == 'TRUE_FALSE')))
    print(length(which(xx$combination == 'FALSE_TRUE')))
    hap_hifi <- c(hap_hifi, xx$h1)
  }else{
    print(length(which(xx$combination == 'TRUE_FALSE')))
    print(length(which(xx$combination == 'FALSE_TRUE')))
    hap_hifi <- c(hap_hifi, xx$h2)
  }
  
}

phase_hifi$hifi <- hap_hifi
phase_hifi$opposite <- as.character(phase_hifi$ont == phase_hifi$hifi)

######################################################################################################################

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


# now I need to check the number of haplotype switches between ont and hifi haplotype per 1 Mb per chr
# make an empty data structure to store data
phase_check <- NULL

for(j in 1:7){
  
  # get the phase for each chr
  chr_phase <- phase_hifi[phase_hifi$chr == paste('chr', j, sep = ''), ]
  
  # chop chr to 1Mb window
  mywindow <- seqlast(0, fai[j, 2], 1000000)
  
  # calculate switch in every window
  myswitch <- c()
  hifivsref <- c()
  ontvsref <- c()
  for(i in 1:(length(mywindow) - 1)){
    
    chr_window <- chr_phase[chr_phase$pos >= mywindow[i] & chr_phase$pos < mywindow[i + 1], ]
    myswitch <- c(myswitch, round(nrow(chr_window[chr_window$opposite == 'FALSE', ])/nrow(chr_window), 2))
    hifivsref <- c(hifivsref, round(nrow(chr_window[chr_window$hifi == 1, ])/nrow(chr_window), 2))
    ontvsref <- c(ontvsref, round(nrow(chr_window[chr_window$ont == 1, ])/nrow(chr_window), 2))
    print(i)
    
  }
  
  # get the middle position in each window
  mystart <- mywindow[1:(length(mywindow)-1)] 
  myend <- mywindow[2:(length(mywindow))]
  
  # get a verctor of chr name with the same length as myswitch and myposition
  mychr <- chr_phase$chr[1:length(myswitch)]
  
  # cbind mychr, mystart, myend, myswitch to form a data frame
  mydf <- cbind.data.frame(mychr, mystart, myend, myswitch, hifivsref, ontvsref)
  
  # rename the colums in the mydf
  colnames(mydf) <- c("chr", "start", "end", "hifivsont", 'hifivsref', 'ontvsref')
  
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
write.table(phase_check, 'phase_check_hifi_vs_ont.txt', quote = F, sep = '\t', col.names = T, row.names = F)
