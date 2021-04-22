library(tidyverse)
###import my count data
filenames <- list.files()
count_list <- vector("list",length = length(filenames))
for (i in seq_along(filenames)) {
  count_list[[i]] <- read_delim(filenames[i],"\t",col_names = F)
}
count_list 
##bind_cols(count_list[[1]],count_list[[2]]) %>% head()
count_data <- tibble(gene = count_list[[1]][[1]])
for (i in seq_along(count_list)) {
  count_data <<- bind_cols(count_data,count_list[[i]][2])
}
head(count_data)
#### import the author's data and metadata
count_done <- read_delim("count_done.csv",col_names = T,delim = ";")
metadata <- read_csv("Run.txt",col_names = T)
metadata <- metadata[,c(1,2,4,7,8,15,18,20)]

rna_metadata <- metadata %>% filter(`Assay Type` == "RNA-Seq" & method  == '1k_cells_RNAseq_Smartseq2') %>% select(c(1,3,4,6))
names(rna_metadata)[4] <- "GEO"
library(stringr)
"GSM4460176\tRNAseq_Foxj1_Olig2_tdT_Uninj_rep1\nGSM4460177\tRNAseq_Foxj1_Olig2_tdT_1dpi_rep1\nGSM4460178\tRNAseq_Foxj1_Olig2_tdT_1dpi_rep2\nGSM4460179\tRNAseq_Foxj1_tdT_Uninj_rep1\nGSM4460180\tRNAseq_Foxj1_tdT_1dpi_rep1\nGSM4460181\tRNAseq_Foxj1_tdT_5dpi_rep1\nGSM4460182\tRNAseq_Foxj1_tdT_1dpi_rep2\nGSM4460183\tRNAseq_Foxj1_tdT_5dpi_rep2\nGSM4460184\tRNAseq_Foxj1_Olig2_tdT_5dpi_rep1\nGSM4460185\tRNAseq_Foxj1_Olig2_tdT_Uninj_rep2\nGSM4460186\tRNAseq_Foxj1_Olig2_tdT_5dpi_rep2\nGSM4460187\tRNAseq_Foxj1_tdT_Uninj_rep2" %>% 
  str_split("\t") %>% unlist() %>% str_split("\n") %>% unlist() -> tmp
tmp[str_starts(tmp,"GSM")]
tmp[str_starts(tmp,"RNA")] %>% str_remove("RNAseq_") %>% str_remove("_rep(1|2)")
rna_metadata$sample <- tmp[str_starts(tmp,"RNA")] %>% str_remove("RNAseq_") 
rna_metadata$sample %>% str_detect("Olig2") -> rna_metadata$Olig2_positive
rna_metadata$injury <- c("Uninj","1dpi","1dpi","Uninj","1dpi","5dpi","1dpi","5dpi","5dpi","Uninj","5dpi","Uninj")
rna_metadata$injury <- factor(rna_metadata$injury,levels = c("Uninj","1dpi","5dpi"))
rna_metadata  <- rna_metadata %>% arrange(Olig2_positive,injury)
count_done %>% colnames()
colnames(count_data) <- c("gene",rna_metadata$sample)
count_data <- count_data %>% select(gene,rna_metadata$sample)

####annotation
library(Mus.musculus)
Mus.musculus %>% keytypes()
Mus.musculus %>% columns()
tmp <- count_data$gene %>% str_remove("\\..*")
length(tmp)
tmp %>% unique() %>% length()
names(count_data)[1] <- "ENSEMBL"
count_data$ENSEMBL <- tmp

genes <- select(Mus.musculus,keys = tmp,keytype = "ENSEMBL",columns = c("ENTREZID","SYMBOL")) %>%unique()
count_data$SYMBOL <- genes$SYMBOL[match(tmp,genes$ENSEMBL)]
count_data$ENTREZID <- genes$ENTREZID[match(tmp,genes$ENSEMBL)]
######## construct a se
library(SummarizedExperiment)
unique(count_data$ENSEMBL) %>% length()
counts <- as.matrix(count_data[2:13])
rownames(counts) <- count_data$ENSEMBL

colData <- DataFrame(rna_metadata)
rownames(colData) <- colnames(counts)
se <- SummarizedExperiment(assays = list(counts = counts),colData = colData)
se
metadata(se) <- count_data[c(1,14,15)]
colData(se)$Olig2_positive <- as.factor(colData(se)$Olig2_positive)

########################################
levels(se$Olig2_positive)
levels(se$injury)
round( colSums(assay(se)) / 1e6, 1 ) 
round( colSums(assay(se)) / 1e6, 1 ) /mean(round( colSums(assay(se)) / 1e6, 1 ) )## may related sizefactor
##########################################construct the DESeqDataSet object
library(DESeq2)
dds <- DESeqDataSet(se,design = ~ Olig2_positive + injury)
##pre-filter
nrow(dds)
keep <- rowSums(counts(dds) >= 1) >= 1
dds <- dds[keep,]
dds
###Normalization
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)
colData(dds)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
colData(rld)
#### calculate the Sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDists
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- vsd$sample
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
########################### PCA plot
plotPCA(vsd, intgroup = c("injury", "Olig2_positive"))
pcaData <- plotPCA(vsd,  c("injury", "Olig2_positive"), returnData = TRUE)
pcaData
#################PCA plot using Generalized PCA
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$dex <- dds$dex
gpca.dat$Olig2 <- dds$Olig2_positive
gpca.dat$injury <- dds$injury
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Olig2, shape = injury)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


###    Differential expression analysis

dds <- DESeq(dds) ## dds was enhanced
dds <- dds[str_starts(rownames(dds),"ENS"),]
res <- results(dds)
res
res <- results(dds, contrast=c("injury","5dpi","Uninj"))
res
##res is a DataFrame,it carries metadata 
mcols(res)
summary(res)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
resLFC1 <- results(dds, lfcThreshold=1,alpha = 0.05)
table(resLFC1$padj < 0.1)
table(resLFC1$padj < 0.1 & abs(resLFC1$log2FoldChange) >1)
##other compared
results(dds, contrast = c("Olig2_positive","TRUE","FALSE"))
####
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
metadata(res)

########## diagnostic plot
##M-A
resultsNames(dds)
library(apeglm)
res <- lfcShrink(dds, coef="injury_5dpi_vs_Uninj", type="apeglm")
plotMA(res, ylim = c(-18, 18))

with(res[genes$ENSEMBL[match("Pdgfra",genes$SYMBOL)], ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, genes$ENSEMBL[match("Pdgfra",genes$SYMBOL)], pos=2, col="dodgerblue")
})
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

#####Gene clustering
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat  <- assay(vsd)[topVarGenes, ]
mat  <- mat - rowMeans(mat)
rownames(mat) <- genes$SYMBOL[match(rownames(mat),genes$ENSEMBL)]
anno <- as.data.frame(colData(vsd)[, c("Olig2_positive","injury")])
pheatmap(mat, annotation_col = anno)


res$symbol <- mapIds(org.Mm.eg.db,
                     keys=rownames(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=rownames(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[1:1000, ]
write.csv(resOrderedDF, file = "results_dds.csv")

#########################Plotting results

dds2 <- dds
genes$SYMBOL[match(rownames(dds2),genes$ENSEMBL)] %>% unique %>% length()


BiocManager::install("tximeta") ### use the genes quantity
library(tximeta)
#tximeta::summarizeToGene(se)
rownames(se) %>% length()
head(genes)
genes %>% subset(genes$SYMBOL == "Olig2",1)

###############
#####################
##################### Plot the genes 
plotCounts(dds, gene = "ENSMUSG00000039830", intgroup=c("injury"))
###this is equal to the below
counts(dds,normalized = T) %>% as.data.frame() %>% filter(rownames(.) ==  "ENSMUSG00000039830")
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = "ENSMUSG00000039830", intgroup = c("injury","Olig2_positive"),returnData = TRUE)
ggplot(geneCounts, aes(x = injury, y = count, color = Olig2_positive)) +
  geom_beeswarm(cex = 3)
######***********************************************
geneCounts <- plotCounts(dds, gene = "ENSMUSG00000029231", intgroup = c("injury","Olig2_positive"),returnData = TRUE)
geneCounts$coor <- c(0,0,1,1,5,5)
geneCounts %>% group_by(coor,Olig2_positive) %>%
 summarise(mean_count = mean(count),sd = sd(count)/sqrt(2)) %>% 
 ggplot(aes(x = coor, y = mean_count, color = Olig2_positive, group = Olig2_positive)) +
 geom_ribbon(aes(ymin = mean_count -sd,ymax = mean_count + sd,fill = Olig2_positive),alpha= 0.3,linetype = 0)  + 
 geom_line(size = 1.2) + scale_x_continuous(breaks = c(0,1,5),labels = c("Uninj","1dpi","5dpi")) + 
  scale_color_manual(values = RColorBrewer::brewer.pal(8,"Set2")[c(8,2)]) + scale_fill_manual(values = RColorBrewer::brewer.pal(8,"Set2")[c(8,2)]) + 
  theme_bw() +  ylab("Normalized expression") +   xlab("Pdgfra")
###*****************************************************
genes$ENSEMBL[match("Pdgfra",genes$SYMBOL)]
selected_genes <- c("Nes","Ptprz1","Sox8","Gdnf","Mki67","Pdgfra","Smarca4","Cspg4","Gpr17")
genes$ENSEMBL[match(selected_genes,genes$SYMBOL)]
library(gridExtra)
plot_rnaseq <- function(gene){
  geneCounts <- plotCounts(dds, gene = genes$ENSEMBL[match(gene,genes$SYMBOL)], intgroup = c("injury","Olig2_positive"),returnData = TRUE)
  geneCounts$coor <- c(0,0,1,1,5,5)
  p <- geneCounts %>% group_by(coor,Olig2_positive) %>%
    summarise(mean_count = mean(count),sd = sd(count)/sqrt(2)) %>% 
    ggplot(aes(x = coor, y = mean_count, color = Olig2_positive, group = Olig2_positive)) +
    geom_ribbon(aes(ymin = mean_count -sd,ymax = mean_count + sd,fill = Olig2_positive),alpha= 0.3,linetype = 0)  + 
    geom_line(size = 1.2) + scale_x_continuous(breaks = c(0,1,5),labels = c("Uninj","1dpi","5dpi")) + 
    scale_color_manual(values = RColorBrewer::brewer.pal(8,"Set2")[c(8,2)]) + scale_fill_manual(values = RColorBrewer::brewer.pal(8,"Set2")[c(8,2)]) + 
    theme_bw() +  ylab("Normalized expression") +   xlab(gene)
  return(p)
}
grid.arrange(plot_rnaseq("Olig2"),plot_rnaseq("Sox10"))
grid.arrange(plot_rnaseq(selected_genes[1]),plot_rnaseq(selected_genes[2]),plot_rnaseq(selected_genes[3]),plot_rnaseq(selected_genes[4]),plot_rnaseq(selected_genes[5]),plot_rnaseq(selected_genes[6]),plot_rnaseq(selected_genes[7]),plot_rnaseq(selected_genes[8]),plot_rnaseq(selected_genes[9]),ncol = 3,nrow = 3)
####################
#####################
#########################





#################################annatation plot
colData(dds)

dds_5dpi_uninj <- dds2[,dds2$Olig2_positive == FALSE & dds2$injury %in% c("Uninj","5dpi") ] 
colData(dds_5dpi_uninj)
design(dds_5dpi_uninj) <- ~ injury
dds_5dpi_uninj <- DESeqDataSet(dds_5dpi_uninj,design = ~injury)
dds_5dpi_uninj <- DESeq(dds_5dpi_uninj)

dds_olig2_tdt <- dds2[,dds2$injury %in% c("5dpi") ] 
colData(dds_olig2_tdt)
design(dds_olig2_tdt) <- ~ Olig2_positive
dds_olig2_tdt <- DESeqDataSet(dds_olig2_tdt,design = ~Olig2_positive)
dds_olig2_tdt <- DESeq(dds_olig2_tdt)



resLFC1_5dpi_uninj <- results(dds_5dpi_uninj, lfcThreshold=1,alpha = 0.05)
table(resLFC1_5dpi_uninj$pvalue <= 0.05)
deg1 <- rownames(dds)[resLFC1_5dpi_uninj$pvalue <= 0.05]
resultsNames(dds_5dpi_uninj)

resLFC1_olig2_tdt<- results(dds_olig2_tdt , lfcThreshold=1,alpha = 0.05)
table(resLFC1_olig2_tdt$pvalue <= 0.05)
deg2 <- rownames(dds_olig2_tdt)[resLFC1_olig2_tdt$pvalue <= 0.05]
resultsNames(dds_olig2_tdt)

deg3 <- intersect(deg1,deg2)

library(apeglm)
res_5dpi <- lfcShrink(dds, coef="injury_5dpi_vs_Uninj", type="apeglm")
plotMA(res_5dpi, ylim = c(-18, 18))

vsd_1 <- vst(dds2,blind = F)
library("genefilter")
mat  <- assay(vsd_1)[match(deg1,rownames(assay(vsd_1))), ]
mat  <- mat - rowMeans(mat)
mat <- mat %>% as.data.frame() %>% select(rownames(anno))
anno <- as.data.frame(colData(vsd_1)[, c("injury","Olig2_positive")])
anno <- anno[order(anno$injury),]
p1 <- pheatmap(mat, annotation_col = anno,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames = F,color = colorRampPalette(c("navy","white","firebrick3"))(300),show_colnames = F)

mat  <- assay(vsd_1)[match(deg2,rownames(assay(vsd_1))), ]
mat  <- mat - rowMeans(mat)
mat <- mat %>% as.data.frame() %>% select(rownames(anno))
anno <- as.data.frame(colData(vsd_1)[, c("injury","Olig2_positive")])
anno <- anno[order(anno$injury),]
p2 <- pheatmap(mat, annotation_col = anno,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames = F,color = colorRampPalette(c("navy","white","firebrick3"))(300),show_colnames = F)

mat  <- assay(vsd_1)[match(deg3,rownames(assay(vsd_1))), ]
mat  <- mat - rowMeans(mat)
mat <- mat %>% as.data.frame() %>% select(rownames(anno))
anno <- as.data.frame(colData(vsd_1)[, c("injury","Olig2_positive")])
anno <- anno[order(anno$injury),]
p3 <- pheatmap(mat, annotation_col = anno,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames = F,color = colorRampPalette(c("navy","white","firebrick3"))(300),show_colnames = T)
p3

# This is a useful construction when users just want to compare
# specific groups which are combinations of variables.

dds$group <- factor(paste0(dds$genotype, dds$condition))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)

# the condition effect for genotypeIII
results(dds, contrast=c("group", "IIIB", "IIIA"))
