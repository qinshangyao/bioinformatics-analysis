setwd("/ifs1/User/shangyao/ncbi/public/sra/seqData/RNA-seq")
list.files()
counts_done <- read.table("analysis_counts.txt",sep="\t",stringsAsFactors = F,header = T)
##head(counts_done)
##tmp <- counts_done$Foxj1_tdT_Uninj_rep1[1:10]
##rm(tmp)
##sub(",","",tmp)
library(tidyverse)
dropcom <- function(x){
  x <- as.numeric(sub(",","",x))
  return(x)
}
counts_done <- map_df(counts_done[-1],dropcom) %>% mutate(gene = counts_done$gene) %>% select("gene",1:12)###load the data has been processed!! 
##########################################################################################################
metadata <- read.csv("metadata.txt")[c("Run","Sample.Name","Assay.Type","Bases","Bytes")]
metadata <- metadata %>% dplyr::filter(Assay.Type == "RNA-Seq")
metadata
GSM <- group_info[[1]] %>% str_detect("GSM")
group_id <- group_info[[1]][!GSM]
GSM_id <- group_info[[1]][GSM]
GSM_id
group_id
metadata <- metadata[1:12,]
metadata$group_info <- group_id
metadata#####construct the metadata
metadata$address <- paste0(metadata$Run,"_quant")
metadata$files <- file.path("/ifs1/User/shangyao/ncbi/public/sra/seqData/RNA-seq/quants", metadata$address,"quant.sf")
file.exists(metadata$files)
metadata$Bases / metadata$Bytes
#################################################################
library("tximeta")
coldata <- metadata
coldata$names <- coldata$Run
rownames(coldata) <- coldata$Run
se <- tximeta(coldata = coldata)##import the count data

##################################################################
dim(se)
head(rownames(se))
library(SummarizedExperiment)
head(assay(se,"counts"))
gse <- summarizeToGene(se)
gse
colData(gse)
rowData(gse)
dim(gse)
assayNames(gse)
head(assay(gse,"abundance"), 3)
colSums(assay(gse))
rowRanges(gse)
seqinfo(rowRanges(gse))
round( colSums(assay(gse)) / 1e6, 1 )
colData(gse)$group_info 
cellcondition <- c("Uninj","Inj_1dpi","Inj_1dpi","Uninj","Inj_1dpi","Inj_5dpi","Inj_1dpi","Inj_5dpi","Inj_5dpi","Uninj","Inj_5dpi","Uninj")
mouse <- c(rep("olig2-tdt",3),rep("ctrl",5),rep("olig2-tdt",3),"ctrl")
colData(gse)$cellcondition <- factor(cellcondition,levels = c("Uninj","Inj_1dpi","Inj_5dpi"))
colData(gse)$mouse <- factor(mouse,levels = c("ctrl","olig2-tdt"))
colData(gse)
######################################
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ mouse + cellcondition)
dds
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
dds[dds$                                                                                                   
#########
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)
################
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

####sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDists
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$mouse, vsd$cellcondition, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( vsd$mouse, vsd$cellcondition, sep = " - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
plotPCA(vsd, intgroup = c("mouse", "cellcondition"))
pcaData <- plotPCA(vsd, intgroup = c("mouse", "cellcondition"), returnData = TRUE)
pcaData

##########
gene_fucus <- c("Olig2","Sox10","Nes","Ptprz1","Sox8","Gdnf","Mki67","Pdgfra","Smarca4","Cspg4","Gpr17")
Ens_id <- rowData(dds)$gene_id[rowData(dds)$symbol %in% gene_fucus]
dds
dds_focus <- dds[rownames(dds) %in% Ens_id,]
assay(dds_focus)
colData(dds_focus)
rowData(dds_focus)
plot2 <- as.data.frame(counts(dds_focus))
plot2$gene <- rowData(dds_focus)$symbol
paste(mouse,cellcondition,sep = "_") %>% rle()

#############
dds <- DESeq(dds)
dds
res <- results(dds)
results(dds,contrast = c("cellcondition","Inj_1dpi","Uninj"))

mcols(res, use.names = TRUE)
summary(res)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = Ens_id[1], intgroup=c("cellcondition"))
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = Ens_id[1], intgroup = c("cellcondition","mouse"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = cellcondition, y = count, color = mouse)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = cellcondition, y = count, group = mouse)) +
    geom_smooth(aes(color = mouse),se = F,method = "loess",level=0.99)
dds[rowData(dds)$symbol %in% gene_fucus ,] %>% counts(normalized = T)
dds[rowData(dds)$symbol %in% "Top2a" ,] %>% counts(normalized = T)
counts_done[counts_done$gene == "Top2a",]
#####
library("airway")
dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))
csvfile <- file.path(dir, "sample_table.csv")
coldata_2 <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata_2
coldata_2 <- coldata_2[1:2,]
coldata_2$names <- coldata_2$Run
coldata_2$files <- file.path(dir, "quants", coldata_2$names, "quant.sf.gz")
file.exists(coldata_2$files)
se2 <- tximeta(coldata_2)
gse2 <- summarizeToGene(se2)
gse
se
