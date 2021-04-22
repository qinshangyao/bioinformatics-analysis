###  import narrowpeak files
setwd("~/login123/ATAC-seq/merge_peaks/narrowPeak")
list.files()
peaks <- dir(".", pattern = "*.narrowPeak",full.names = TRUE)
library(ChIPQC)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE) ### each narrowpeaks had different number of ranges
### import mdata
library(tidyverse)
anno <- read_delim("/ifs1/User/shangyao/login123/ATAC-seq/atac_sample.txt",delim = "\t") %>% select(c(1,9))
anno
names(myPeaks) <- anno$sample
Group <- factor(anno$sample %>% str_remove("_rep.*$"))
##### 
myGRangesList<-GRangesList(myPeaks) 
reduced <- reduce(unlist(myGRangesList)) ### get the 
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
consensusToCount <- reducedConsensus 
seqlevels(consensusToCount)
### filter the ranges in all 12 samples
blklist <- rtracklayer::import.bed("ENCFF226BDM.bed.gz")
blklist
seqlevels(blklist)
####chrM had been removed by samtools view..| grep -v chrM
###consensusToCount <- consensusToCount[!consensusToCount %over% blklist & !seqnames(consensusToCount) %in% "chrM"]
consensusToCount <- consensusToCount[!consensusToCount %over% blklist] 
#### find the overlaps between replicates
library(limma)
as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(starts_with("Foxj1_tdT_5dpi")) %>% 
  vennDiagram(main = "Overlap for Liver open regions")
as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(starts_with("Foxj1_Olig2_tdT_5dpi")) %>% 
  vennDiagram(main = "Overlap for Liver open regions")
as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(starts_with("Foxj1_Olig2_tdT_1dpi")) %>% 
  vennDiagram(main = "Overlap for Liver open regions")
as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(starts_with("Foxj1_Olig2_tdT_U")) %>% 
  vennDiagram(main = "Overlap for Liver open regions")
as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(starts_with("Foxj1_tdT_U")) %>% 
  vennDiagram(main = "Overlap for Liver open regions")
as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(starts_with("Foxj1_tdT_1")) %>% 
  vennDiagram(main = "Overlap for Liver open regions")

library(tidyr)
Group2 <- anno$sample %>% str_remove("_rep.*$")
myPlot <- as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(-consensusIDs) %>% 
  as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  mutate(Group = Group) %>% ggplot(aes(x = PC1, y = PC2,  colour = Group)) + geom_point(size = 5)
myPlot

library(Rsubread)
occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% 
  rowSums
table(occurrences) %>% rev %>% cumsum
consensusToCount <- consensusToCount[occurrences >= 2, ]
consensusToCount


#####load bam file 
bamsToCount <- dir("/ifs1/User/shangyao/login123/ATAC-seq/bowite2_result/last_bam", full.names = TRUE, pattern = "*.\\.bam$")
# indexBam(bamsToCount)
regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount), 
                                            start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount), 
                             Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))
fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, isPairedEnd = TRUE, 
                           countMultiMappingReads = FALSE, maxFragLength = 100)
myCounts <- fcResults$counts
colnames(myCounts) <- anno$sample
##save(myCounts, file = "countsFromATAC.RData")
########################################################
library(DESeq2)
metaData <- data.frame(Group, row.names = colnames(myCounts),olig2 = c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE),injury = c("5dpi","1dpi","5dpi","5dpi","Uninj","5dpi","Uninj","1dpi","1dpi","Uninj","Uninj","1dpi"))
metaData$olig2 <- factor(metaData$olig2,levels = c(FALSE,TRUE))
metaData$injury <- factor(metaData$injury,levels = c("Uninj","1dpi","5dpi"))
metaData

atacDDS <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount)
############################### Now we get the DESeqDataSet !!!!!!!
atacDDS <- DESeq(atacDDS)### enhance the DESeqDataSet!!
atac_Rlog <- rlog(atacDDS) ### normalization
counts(atacDDS) %>% head()
rowData(atacDDS) %>% head()
colData(atacDDS)
##pca plot!! fig2D
plotPCA(atac_Rlog, intgroup = "Group", ntop=500)
pcaData <- plotPCA(atac_Rlog, intgroup = "Group", returnData = TRUE,ntop = 500)
pcaData$olig2 <- metaData$olig2
pcaData$injury <- metaData$injury
attr(pcaData, "percentVar")
percentVar <- round(100 * attr(pcaData, "percentVar"))
#####
ggplot(pcaData, aes(x = -PC1, y = PC2, color = olig2, shape = injury)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + scale_color_manual(values = c("grey","red")) + 
  ggtitle("Global chromatin accessibility") + theme_bw()


##################
install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
install("tracktables")
library(tracktables)

atacDDS$Group
Uninj_Olig2vsTdt <- results(atacDDS, c("Group", "Foxj1_Olig2_tdT_Uninj", "Foxj1_tdT_Uninj"), format = "GRanges")
dpi1_Olig2vsTdt <- results(atacDDS, c("Group", "Foxj1_Olig2_tdT_1dpi", "Foxj1_tdT_1dpi"), format = "GRanges")
dpi5_Olig2vsTdt <- results(atacDDS, c("Group", "Foxj1_Olig2_tdT_5dpi", "Foxj1_tdT_5dpi"), format = "GRanges")

Uninj_Olig2vsTdt<- Uninj_Olig2vsTdt[order(Uninj_Olig2vsTdt$pvalue)]
dpi1_Olig2vsTdt <- dpi1_Olig2vsTdt[order(dpi1_Olig2vsTdt$pvalue)]
dpi5_Olig2vsTdt <- dpi5_Olig2vsTdt[order(dpi5_Olig2vsTdt$pvalue)]

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
toOverLap <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, 500,500)

dpi5_Olig2vsTdt <- dpi5_Olig2vsTdt[(!is.na(dpi5_Olig2vsTdt$padj) & dpi5_Olig2vsTdt$padj < 0.05) & dpi5_Olig2vsTdt %over% toOverLap, ]
makebedtable(dpi5_Olig2vsTdt, "dpi5_Olig2vsTdt.html", ".")

export.bed(dpi5_Olig2vsTdt, "dpi5_Olig2vsTdt.bed")

######
library(clusterProfiler)
install("ChIPseeker")
library(ChIPseeker)
anno_dpi5_Olig2vsTdt  <- annotatePeak(dpi5_Olig2vsTdt, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)

go1 <- enrichGO(as.data.frame(as.GRanges(anno_dpi5_Olig2vsTdt )[as.GRanges(anno_dpi5_Olig2vsTdt )$log2FoldChange > 
                                                                     0])$geneId, OrgDb = "org.Mm.eg.db", ont = "BP", maxGSSize = 5000)
go2 <- enrichGO(as.data.frame(as.GRanges(anno_dpi5_Olig2vsTdt )[as.GRanges(anno_dpi5_Olig2vsTdt )$log2FoldChange < 
                                                                     0])$geneId, OrgDb = "org.Mm.eg.db", ont = "BP", maxGSSize = 5000)

library(data.table)
??datatable
head(go1, 10) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% data.table(elementId = "goEle1")
head(go2, ) %>% dplyr::select(ID, Description, pvalue, p.adjust) %>% data.table(elementId = "goEle2")


###########Cutting sites from ATAC-seq data
library(MotifDb)
library(Biostrings)
install('BSgenome.Hsapiens.UCSC.hg19')
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)

SOX10 <- query(MotifDb, c("SOX10"))
sox10 <- as.list(SOX10)
rm(SOX10)
myRes <- matchPWM(sox10[[2]], BSgenome.Mmusculus.UCSC.mm10)
toCompare <- GRanges("chr20", ranges(myRes))
####export bed file
mybed <- data.frame(Uninj_Olig2vsTdt[1:2000,] %>% seqnames() ,Uninj_Olig2vsTdt[1:2000,] %>% ranges())
mybed <- mybed[1:3]
names(mybed) <- c("chr","start","end")
mybed


######QC
library("MotifDb")
install("ATACseqQC")
library(ATACseqQC)
bamfile1 <- "~/login123/ATAC-seq/bowite2_result/last_bam/SRR11510933.last.bam"
bamfile2 <- "~/login123/ATAC-seq/bowite2_result/last_bam/SRR11510930.last.bam"
bamfile.labels1 <- gsub(".last.bam", "", basename(bamfile1))
bamfile.labels2 <- gsub(".last.bam", "", basename(bamfile2))

fragSize1 <- fragSizeDist(bamfile1, bamfile.labels1)
fragSize2 <- fragSizeDist(bamfile1, bamfile.labels2)