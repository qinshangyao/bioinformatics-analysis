setwd("/ifs1/User/shangyao/ncbi/public/sra/mapping_result/mapping_2")
coldata <- read.table(file = "/ifs1/User/shangyao/ncbi/public/sra/metadata.txt",stringsAsFactors = F)

head(coldata)
names(coldata) <- c("run","sample")
#coldata should contain "names" and "files"
coldata$names <- coldata$run
coldata$files <- file.path("/ifs1/User/shangyao/ncbi/public/sra/mapping_result/mapping_2",coldata$names,"quant.sf")
file.exists(coldata$files)

####load quant.sf file
library("tximeta")
se <- tximeta(coldata)#couldn't find matching transcriptome, returning non-ranged SummarizedExperiment

#### Linked transcriptomes
indexDir <- file.path("/ifs1/User/shangyao/ncbi/public/sra/danRer_genome/danRer_salmon_index")
annofile <- file.path("/ifs1/User/shangyao/ncbi/public/sra/danRer_genome/salmon_linked")
list.files(annofile)
fasta <- file.path(annofile,"transcripts.fa")
gtf <- file.path(annofile,"danRer11.ensGene.gtf")
tmp <- "/ifs1/User/shangyao/.cache/tximeta"
jsonFile <- file.path(tmp, paste0(basename(indexDir), ".json"))
makeLinkedTxome(indexDir=indexDir,
                source="UCSC", organism="Zebrafish",
                release="11", genome="GRCz11",
                fasta=fasta, gtf=gtf,
                write = T,
                jsonFile = jsonFile)
####load the data again
se <- tximeta(coldata)##still not finding matching transcriptome!! i don't know why!!!!

### transform to a sce!!!
library(SingleCellExperiment)
sce <- as(se, "SingleCellExperiment")

library(stringr)
#head(rownames(se)) %>% str_remove("\\.[0-9]*")
rownames(sce) <- rownames(sce) %>% str_remove("\\.[0-9]*")
dim(assay(sce))
colnames(colData(sce))
colnames(rowData(sce))
counts(sce)[1:6,1:6]

####Adding more assays
sce <- scater::logNormCounts(sce)
logcounts(sce)[1:4,1:4]
###adding more colDatasce <- scater::addPerCellQC(sce)
sce <- scater::addPerCellQC(sce)
colData(sce)
#sce$more_stuff <- runif(ncol(sce)) ## we can add some var in our need.
#colnames(colData(sce))

### add more rowData
sce <- scater::addPerFeatureQC(sce)
rowData(sce)

library(org.Dr.eg.db)
keys <- rownames(sce)
keytypes(org.Dr.eg.db)
anno <- select(org.Dr.eg.db,keys = keys,keytype = "ENSEMBLTRANS",column = c("REFSEQ","ENTREZID","SYMBOL"),multiVals ="first")
head(anno)

tx2gene <- anno[c(1,3)] %>% unique()
tx2gene <- tx2gene[match(rownames(sce),tx2gene$ENSEMBLTRANS),]
head(tx2gene)

library(TxDb.Drerio.UCSC.danRer11.refGene)
txdb <- TxDb.Drerio.UCSC.danRer11.refGene
seqlevels(txdb)
#seqlevels(txdb) <- "chr15"
#seqlevels(txdb) <- seqlevels0(txdb)

keytypes(txdb)
select(txdb, keys = keys, columns="GENEID", keytype="TXID")
GR <- transcripts(txdb)
GR$geneid <- anno$ENTREZID[match(GR$tx_name,anno$REFSEQ)]
GR$symbol <- anno$SYMBOL[match(GR$tx_name,anno$REFSEQ)]
head(GR)

## construct a gene_level tpm matrix
sce
tpm <- assay(sce,"abundance")
dim(tpm)
tpm <- as.data.frame(tpm)
tpm$geneid <- tx2gene$ENTREZID
head(tpm[,1:4])
library(tidyverse)
gene_tpm <- tpm %>% gather(key = "cell",value = "TPM",-geneid) %>% group_by(geneid,cell) %>% summarise(TPM_gene = sum(TPM,na.rm = T)) 
gene_tpm <- gene_tpm %>%　spread(key = "cell",value = "TPM_gene")
gene_matrix <- as.matrix(gene_tpm[-1])
row.names(gene_matrix) <- gene_tpm$geneid
colnames(gene_matrix)
gene_matrix[1:4,1:4]
dim(gene_matrix)
gene_matrix <- gene_matrix[! is.na(rownames(gene_matrix)),]
dim(gene_matrix)

row.names(gene_matrix)
anno_gene <- select(org.Dr.eg.db,row.names(gene_matrix),keytype = "ENTREZID",column = c("ENTREZID","ENSEMBL","ENSEMBLTRANS","SYMBOL"),multiVals ="first")
anno_gene <- anno_gene[match(row.names(gene_matrix),anno_gene$ENTREZID),]

gce <-  SingleCellExperiment(assays = list(counts = gene_matrix),colData = coldata[1:3])
gce
rowData(gce) <- anno_gene
rowRanges(gce)

gce <- scran::computeSumFactors(gce)
gce <- scater::addPerCellQC(gce)
colData(gce)
gce <- scater::addPerFeatureQC(gce)
rowData(gce)

library(scater)
library(scran)

GR <- GR[!is.na(GR$geneid),]
####QC
rowrange <- rowRanges(gce)
for (geneid in rownames(gce)) {
  rowrange[[geneid]] <- GR[GR$geneid ==geneid,]
}



####seurat object
library(Seurat) 
library(dplyr) 
library(patchwork)
rm(rowrange,qcstats,rd)

cell384 <- CreateSeuratObject(counts = counts(gce), project = "cell384", min.cells = 15, min.features = 500)
cell384
#cell384[["percent.mt"]] <- PercentageFeatureSet(cell384, pattern = "^mt-") 
head(cell384@meta.data, 5) 
VlnPlot(cell384, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
plot1 <- FeatureScatter(cell384, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cell384, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
cell384 <- NormalizeData(cell384, normalization.method = "LogNormalize", scale.factor = 10000)
cell384  <- FindVariableFeatures(cell384 , selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cell384), 10)
plot1 <- VariableFeaturePlot(cell384)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(cell384)
cell384 <- ScaleData(cell384, features = all.genes)
cell384 <- RunPCA(cell384, features = VariableFeatures(object =cell384))
print(cell384[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cell384, dims = 1:2, reduction = "pca")
DimPlot(cell384, reduction = "pca")
cell384 <- JackStraw(cell384, num.replicate = 100)
cell384 <- ScoreJackStraw(cell384, dims = 1:20)

JackStrawPlot(cell384, dims = 1:15)
ElbowPlot(cell384)

cell384 <- FindNeighbors(cell384, dims = 1:20)
cell384 <- FindClusters(cell384, resolution = 1)

cell384 <- RunTSNE(cell384,dims = 1:20)
DimPlot(cell384, reduction = "tsne")

cell384 <- RunUMAP(cell384, dims = 1:20)
DimPlot(cell384, reduction = "umap")

cluster1.markers <- FindMarkers(cell384, ident.1 = 1, min.pct = 0.25)
cluster1.markers$gene <- anno_gene$SYMBOL[match(rownames(cluster1.markers),anno_gene$ENTREZID)]
head(cluster1.markers, n = 5)

cell384.markers <- FindAllMarkers(cell384, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- cell384.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
markers$symbol <- anno_gene$SYMBOL[match(markers$gene,anno_gene$ENTREZID)]
markers

mygene <- c("sox10","olig2","nkx2.2a","ppp1r14bb","mbpa","plp1a")
mygene <- anno_gene$ENTREZID[match(mygene,anno_gene$SYMBOL)]
mygene

new.cluster.ids <- c("OPC#1", "VLMC1", "VLMC2", "mOL", "OPC#3", "OPC#4", "OPC#2")
names(new.cluster.ids) <- levels(cell384)
cell384 <- RenameIdents(cell384, new.cluster.ids)


DimPlot(cell384, reduction = "umap", label = TRUE, pt.size = 0.9) + NoLegend()
FeaturePlot(cell384, features = mygene)
VlnPlot(cell384,features = mygene)
# Reorder identity classes
levels(x = cell384)
levels(x = cell384) <- c("OPC#1","OPC#2", "OPC#3", "OPC#4","mOL","VLMC1", "VLMC2")
levels(x = cell384)
DimPlot(cell384, reduction = "umap", label = TRUE, pt.size = 0.9) + NoLegend()
FeaturePlot(cell384, features = mygene)
VlnPlot(cell384,features = mygene)
VlnPlot(cell384,features = gene2)
DimPlot(cell384, reduction = "umap", label = TRUE, pt.size = 0.9) + NoLegend()
table(cell384$RNA_snn_res.1)

#############################
s.genes <- cc.genes$s.genes
s.genes <- s.genes %>% stringr::str_to_lower(s.genes) 
s.genes <- anno_gene$ENTREZID[match(s.genes,anno_gene$SYMBOL)]
s.genes
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- g2m.genes %>% stringr::str_to_lower(g2m.genes) 
g2m.genes <- anno_gene$ENTREZID[match(g2m.genes,anno_gene$SYMBOL)]

cell354 <- CellCycleScoring(cell384, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(cell354[[]])
RidgePlot(cell354, features = gene2[c(3,4)], ncol = 2)
cell354<- RunPCA(cell354, features = c(s.genes, g2m.genes))
DimPlot(cell354)


cell_no_cc <- ScaleData(cell354, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(cell354))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
cell_no_cc<- RunPCA(cell_no_cc, features = VariableFeatures(cell_no_cc), nfeatures.print = 10)
# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
cell_no_cc <- RunPCA(cell_no_cc, features = c(s.genes, g2m.genes))
DimPlot(cell_no_cc)
############################################
opc <- cell384[,cell384$RNA_snn_res.1 != 1]
opc <- opc[,opc$RNA_snn_res.2 != 2]

opc <- cell384[,Idents(cell384) != "VLMC2"]
opc <- opc[,Idents(opc) != "VLMC1"]
opc

DimPlot(opc, reduction = "umap")
FeaturePlot(opc, features = mygene)
FeaturePlot(opc, features =gene2)

gene2 <- c("cspg4","ptprz1a","pcna","mki67","gpr17","myrf")
gene2 <- anno_gene$ENTREZID[match(gene2,anno_gene$SYMBOL)]
gene2

FeaturePlot(cell384, features = gene2)
VlnPlot(cell384,features = gene2)
VlnPlot(cell384, features = gene2, slot = "counts", log = TRUE)

VlnPlot(opc,features = mygene)
VlnPlot(opc,features = gene2)

opc.markers <- FindAllMarkers(opc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
opc.markers <- opc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
opc.markers$symbol <- anno_gene$SYMBOL[match(opc.markers$gene,anno_gene$ENTREZID)]
opc.markers

DoHeatmap(opc, features = opc.markers$gene)

# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)
#cluster1.markers <- FindMarkers(opc, ident.1 = 1, min.pct = 0.25)
saveRDS(cell384, file = "./cell384_final.rds")
saveRDS(opc,file = "./opc_final.rds")

glist <- c("nes","cd31","cd248")
glist <- c("100150939","569386")
VlnPlot(cell384,features ="PC_1",pt.size = 0)
FeaturePlot(cell384, features = "PC_2")


########################## Monocle3
library(monocle3)
cds <- list(gene_matrix=gene_matrix[anno_gene$ENTREZID,rownames(cell354@meta.data)],
                         cell_metadata = cell354@meta.data,
                         gene_metadata = anno_gene)
saveRDS(cds,file = "cds_3files.rds")
