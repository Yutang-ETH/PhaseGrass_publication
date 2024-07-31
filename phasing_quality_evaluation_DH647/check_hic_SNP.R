# 13.04.2024
# Hi-C phased SNP density or distribution
setwd("P:/Yutangchen/DH647/DH647_nextpolish/phase_grass/hapcut2/vcf_concatenate")

library(stringr)

# import phase
phase <- read.table("all.hapcut.phased.largest.pos.txt", header = F, stringsAsFactors = F)

# import chr length
fai <- read.table("../../asm/kyuss.nextdenovo.juicer.fasta.chr.fa.fai", header = F, stringsAsFactors = F)


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
  
  num_snp <- c()
  for(i in 1:(length(mywindow) - 1)){
    
    chr_window <- chr_phase[chr_phase$V2 >= mywindow[i] & chr_phase$V2 < mywindow[i + 1], ]
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
  mysnpdensity <- num_snp
  
  # cbind mychr, myposition, myswitch to form a data frame
  mydf <- cbind.data.frame(mychr, myposition, mystart, myend, num_snp)
  
  # rename the colums in the mydf
  colnames(mydf) <- c("chr", "pos", "start", "end", "snp_number")
  
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
write.table(phase_check, 'HiC_phased_SNP.txt', quote = F, sep = '\t', col.names = T, row.names = F)
