# make pedigree file
args <- commandArgs(trailingOnly=TRUE)

sample_name = args[1]
p1 = args[2] # either 1 or 2
p2 = args[3] # either 1 or 2, the number has to be different to p1

# import fiels
mysample <- read.table(sample_name, header = F, stringsAsFactors = F)

# how many samples I have
sample_number <- length(mysample$V1) # 282

# make the pedigree file for lep-map3
myline1 <- c("CHR", "POS", rep("SR", sample_number))
myline2 <- c("CHR", "POS", mysample$V1)
myline3 <- c("CHR", "POS", "0", "0", rep(mysample$V1[as.numeric(p1)], (sample_number - 2)))
myline4 <- c("CHR", "POS", "0", "0", rep(mysample$V1[as.numeric(p2)], (sample_number - 2)))
myline5 <- c("CHR", "POS", p1, p2, rep("0", (sample_number-2))) 
myline6 <- c("CHR", "POS", rep("0", sample_number))

mypedi <- rbind.data.frame(myline1, myline2, myline3, myline4, myline5, myline6)

write.table(mypedi, "pedigree.txt", col.names = F, row.names = F, quote = F, sep = "\t")