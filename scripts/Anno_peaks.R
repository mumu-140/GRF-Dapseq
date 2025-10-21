#!/usr/bin/Rscript

library(ChIPseeker)
library(ggplot2)
library(GenomicFeatures)

Args <- commandArgs()
gff_f=Args[6]
idr_f=Args[7]
# 
p2=Args[8]
spst=Args[9]
dost=Args[10]

# input
in_idr_f=paste0(p2,"/",idr_f)

out_idr=paste0(p2,"/",idr_f,"2")

out_peak=paste0(p2,"/",idr_f,"_peak.txt")

txdb <- makeTxDbFromGFF(gff_f)
peak_idr <- read.table(in_idr_f,sep = "\t", header = F)
colnames(peak_idr) <- c("peak_chr","chrStart","chrEnd","name","score","strand","signalValue","p-value","q-value","summit","localIDR","globalIDR","rep1_chrStart","rep1_chrEnd","rep1_signalValue","rep1_summit","rep2_chrStart","rep2_chrEnd","rep2_signalValue","rep2_summit")
write.table(peak_idr,out_idr, sep = "\t", quote = FALSE, row.names = FALSE)
peak_idr <- readPeakFile(out_idr)2

# define promoter region (-2000 +2000)
peakAnno1 <- annotatePeak(peak_idr, tssRegion = c(-2000, 2000), TxDb = txdb) 

write.table(peakAnno1, file = out_peak, sep = '\t', quote = FALSE, col.names = T, row.names = FALSE)
