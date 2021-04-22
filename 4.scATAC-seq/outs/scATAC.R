################scATAC-seq
setwd("~/login123/scATAC-seq/EP/outs")
list.files()
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)

counts <- Read10X_h5(filename = "./filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = 'fragments.tsv.gz',
  min.cells = 0,
  min.features = 200
)

ep <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

ep
ep[["peaks"]]
granges(ep)

# extract gene annotations from EnsDb
annotations_2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
# add the gene information to the object
Annotation(ep) <- annotations
rm(annotations_2)


# compute nucleosome signal score per cell
ep <- NucleosomeSignal(object = ep)

# compute TSS enrichment score per cell
ep <- TSSEnrichment(object = ep, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
ep$pct_reads_in_peaks <- ep$peak_region_fragments / ep$passed_filters * 100
ep$blacklist_ratio <- ep$blacklist_region_fragments / ep$peak_region_fragments
ep$high.tss <- ifelse(ep$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(ep, group.by = 'high.tss') + NoLegend()

ep$nucleosome_group <- ifelse(ep$nucleosome_signal >1, 'NS > 1', 'NS < 1')
ep$nucleosome_group %>%　table()
VlnPlot(
  object = ep,
  features = c('blacklist_ratio', 'nucleosome_signal','peak_region_fragments','pct_reads_in_peaks','TSS.enrichment'),
  pt.size =0,
  ncol = 5
)



ep <- subset(
  x = ep,
  subset = peak_region_fragments > 2000 &
    peak_region_fragments < 50000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 10 &
    TSS.enrichment > 2
)
ep

ep <- RunTFIDF(ep)
ep <- FindTopFeatures(ep, min.cutoff = 'q0')
ep <- RunSVD(ep)
DepthCor(ep)
ElbowPlot(ep,reduction = "lsi")

ep <- RunUMAP(object = ep, reduction = 'lsi', dims = 2:12)
ep <- FindNeighbors(object = ep, reduction = 'lsi', dims = 2:12)
ep <- FindClusters(object = ep,verbose = FALSE, algorithm = 3,resolution = 0.5)
DimPlot(object = ep, label = TRUE) + NoLegend()

#marker5 <- FindMarkers(ep, ident.1 = 5, min.pct = 0.25)
#head(marker5)
#markers <- FindAllMarkers(ep, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
############
gene.activities <- GeneActivity(ep)
ep[['RNA']] <- CreateAssayObject(counts = gene.activities)
ep <- NormalizeData(
  object = ep,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ep$nCount_RNA)
)

DefaultAssay(ep) <- 'RNA'

FeaturePlot(
  object = ep,
  features = c('Acot5',"Cd24a","Capsl", "Gfap","Slc7a10", "Gpr17", "Egfl7", "Foxd1"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4,
  label = F
)

library(dplyr)
markers <- FindAllMarkers(ep, only.pos = TRUE, min.pct = 0.8, logfc.threshold = 0.5)
markers <-  markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)


  #####Integrating with scRNA-seq data
#To help interpret the scATAC-seq data, we can classify cells based on an scRNA-seq experiment from the same biological system
#We aim to identify shared correlation patterns in the gene activity matrix and scRNA-seq dataset to identify matched biological
#states across the two modalities. This procedure returns a classification score for each cell for each scRNA-seq-defined cluster label
# Load the pre-processed scRNA-seq data for eps
ep_rna <- readRDS("~/ep_opc/glia2_scRNA.rds")

table(Idents(ep_rna))
#ep_rna <- ep_rna[,ep_rna$class != "PeripheralGlia"]
#ep_rna <- ep_rna[,ep_rna$class != "Excluded"]
table(ep_rna$class)
ep_rna$Idents <- Idents(ep_rna)

DimPlot(ep_rna, reduction = "umap", label = TRUE, pt.size = 0.1)+ NoLegend()
FeaturePlot(ep_rna,features = c('Acot5',"Cd24a","Capsl", "Gfap","Slc7a10", "Gpr17", "Egfl7", "Foxd1"),pt.size = 0.1,max.cutoff = 'q95',ncol = 4,label = F)

all_genes <- ep_rna@assays$RNA %>% rownames()
gene012 <- all_genes[match(markers$gene[markers$cluster %in% c(0,1,2)],all_genes)] %>% na.omit()
object <- AddModuleScore(object =ep_rna, features = list(gene012),name = "gene012",)
object
FeaturePlot(object = object, features = "gene0121")

gene3 <- all_genes[match(markers$gene[301:400],all_genes)] %>% na.omit()
object <- AddModuleScore(object =ep_rna, features = gene3[1:10],name = "gene3")
FeaturePlot(object = object, features = "gene31")

gene4<- all_genes[match(markers$gene[401:500],all_genes)] %>% na.omit()
object <- AddModuleScore(object =ep_rna, features = gene3[1:10],name = "gene4")
FeaturePlot(object = object, features = "gene41")

gene5<- all_genes[match(markers$gene[501:600],all_genes)] %>% na.omit()
object <- AddModuleScore(object =ep_rna, features = gene3[1:10],name = "gene5")
FeaturePlot(object = object, features = "gene51")

gene6<- all_genes[match(markers$gene[601:700],all_genes)] %>% na.omit()
object <- AddModuleScore(object =ep_rna, features = gene3,name = "gene6")
FeaturePlot(object = object, features = "gene61")

gene7<- all_genes[match(markers$gene[701:800],all_genes)] %>% na.omit()
object <- AddModuleScore(object =ep_rna, features = gene3,name = "gene7")
FeaturePlot(object = object, features = "gene71")

cluster1.markers <- FindMarkers(ep, ident.1 = c(0,1,2), min.pct = 0.25)
cluster1.markers <- cluster1.markers %>% top_n(n = 100, wt = avg_logFC)
gene1<- all_genes[match(rownames(cluster1.markers),all_genes)] %>% na.omit()
object <- AddModuleScore(object =ep_rna, features = gene1,name = "gene1")
FeaturePlot(object = object, features = "gene11")

plot1 <- DimPlot(
  object = ep_rna,
  group.by = 'Idents',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = ep,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

######
DefaultAssay(ep) <- "peaks"
TSSPlot(ep,group.by = ) + NoLegend()
cluster.ids <- c("EP1","EP2","EP3","ACNT1", "ACNT2", "OPC", "VEC", "PER")
levels(ep)
names(cluster.ids) <- levels(ep)
ep <- RenameIdents(ep, cluster.ids)

VlnPlot(
  object = ep,
  features = c('blacklist_ratio', 'nucleosome_signal','peak_region_fragments','pct_reads_in_peaks','TSS.enrichment'),
  
  pt.size =0,
  ncol = 5
)
####################
ep$nucleosome_group <- ifelse(ep$nucleosome_signal > 1, 'NS > 1', 'NS < 1')
FragmentHistogram(object = ep,
  region = "chr1-1-100000000",
  group.by = NULL)
####Find differentially accessible peaks between clusters
DefaultAssay(ep) <- 'peaks'

da_peaks <- FindAllMarkers(
  object = ep,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
head(da_peaks)

da_peaks %>%filter(p_val < 0.01 & abs(avg_logFC) > 0.25) %>% dim()
da_peaks %>% filter(cluster == "OPC") %>% dim()
top100 <- da_peaks  %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC) 
ep_scaled <- ScaleData(object = ep, features = rownames(ep))
DoHeatmap(ep_scaled, features = top100$gene,) + NoLegend()

#rownames(ep@assays$peaks)
#match(top10$gene[1:4],rownames(ep@assays$peaks))
#######################
####Find  Differentially accessible genes between clusters
DefaultAssay(ep) <- 'RNA'
markers <- FindAllMarkers(ep, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
top100_gene <- markers  %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
ep_scaled <- ScaleData(object = ep, features = rownames(ep))
DoHeatmap(ep_scaled, features = top100_gene$gene) + NoLegend()

genelist1 <- "Dnajc6 Celsr1 Zfp704 Slc7a10 Atp1a2 Sox10 Flt1 Apod"
library(stringr)
genelist1 %>% str_split(" ") %>% unlist() -> genelist1
DefaultAssay(ep) <- "RNA"
ExpressionPlot(
  ep,
  genelist1,
  slot = "data"
)
DefaultAssay(ep) <- "peaks"
CoveragePlot(
  object = ep,
  region = genelist1
)
##################################################################
# set plotting order
levels(ep)
DefaultAssay(ep) <- "peaks"
CoveragePlot(
  object = ep,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)

coverage_gene <- c("Ascl1","Sox8","Olig1","Olig2")
coverage_gene <- annotations[annotations$gene_name %in% coverage_gene,]
coverage_gene[coverage_gene$type == "cds"]

DefaultAssay(ep) <- "peaks"
CoveragePlot(
  object = ep,
  region = c("chr10-87492393-87493088","chr17-25567333-25570214","chr16-91269877-91270659","chr16-91190000-91240000"),
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1
)

VlnPlot(ep_rna,features = c("Ascl1","Sox8","Olig1","Olig2"),ncol = 1,pt.size = 0) + NoLegend()

genelist2 <- markers %>% 
  group_by(cluster) %>% top_n(n = 1, wt = avg_logFC) 
FeaturePlot(ep_rna,features = genelist2$gene)



#####sox10 enhancer
sox10_region <- "chr15-79192518-79194107"
?FeaturePlot()
DefaultAssay(ep) <- 'peaks'
ep
FeaturePlot(
  object = ep,
  features = sox10_region,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 4,
  label = F
)

##################motif analysis

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)

DimPlot(ep, label = TRUE, pt.size = 0.1) + NoLegend()
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)
# add motif information
ep <- AddMotifs(
  object = ep,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

da_peaks_2 <- FindAllMarkers(
  object = ep,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
da_peaks_2[da_peaks_2$p_val < 0.005, ] %>% group_by(cluster) %>% filter(p_val<0.005)%>% top_n(n = 100, wt = avg_logFC) -> top.da.peak2
top.da.peak2$gene
top.da.peak <- rownames(da_peaks_2[da_peaks_2$p_val < 0.005, ])
# test enrichment
DefaultAssay(ep) <- 'peaks'
enriched.motifs <- FindMotifs(
  object = ep,
  features = top.da.peak2$gene
)

MotifPlot(
  object = ep,
  motifs = head(rownames(enriched.motifs))
)

ep <- RunChromVAR(
  object = ep,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(ep) <- 'chromvar'

# look at the activity of sox10
FeaturePlot(
  object = ep,
  features = "MA0442.2",
  min.cutoff = 'q90',
  max.cutoff = 'q99',
  pt.size = 0.1
)

#####fig2D
top_100_motif <- da_peaks_2 %>%  group_by(cluster) %>% filter(p_val<0.005)%>% top_n(n = 100, wt = avg_logFC)
DefaultAssay(ep) <- 'peaks'
ep_scaled <- ScaleData(object = ep, features = rownames(ep))
DoHeatmap(ep_scaled,features = top_100_motif$gene)

DefaultAssay(ep) <- 'peaks'
ep_motif <- (top.da.peak2 %>% filter(cluster %in% c("EP1","EP2","EP3")))$gene
enriched.motifs_ep <- FindMotifs(
  object = ep,
  features = ep_motif
)
ac_motif <- (top.da.peak2 %>% filter(cluster %in% c("ACNT1","ACNT2")))$gene
enriched.motifs_ac <- FindMotifs(
  object = ep,
  features = ac_motif
)

opc_motif <- (top.da.peak2 %>% filter(cluster %in% c("OPC")))$gene
enriched.motifs_opc <- FindMotifs(
  object = ep,
  features = opc_motif
)

other_motif <- (top.da.peak2 %>% filter(cluster %in% c("PER","VEC")))$gene
enriched.motifs_other <- FindMotifs(
  object = ep,
  features = other_motif
)




DefaultAssay(ep) <- 'chromvar'
ep_scaled <- ScaleData(object = ep, features = rownames(ep))
myfeature <- c(rownames(enriched.motifs_ep[1:100,]),rownames(enriched.motifs_ac[1:100,]),rownames(enriched.motifs_opc[1:100,]),rownames(enriched.motifs_other[1:100,]))
ep_scaled <- ScaleData(object = ep, features =myfeature)
DoHeatmap(ep_scaled,features =myfeature)
DoHeatmap(ep_scaled,features =rownames(enriched.motifs_opc[1:100,]))
enriched.motifs_opc
"OLIG1 OLIG2 OLIG3 SOX8 ASCL1 SOX10 SOX2 NKX2-2 NKX2-3 SOX9" %>% str_split(" ") %>% unlist() -> opc_related_gene
rownames(enriched.motifs_opc)[match(opc_related_gene,enriched.motifs_opc$motif.name)] -> opc_related_gene
ep_scaled <- ScaleData(object = ep, features =opc_related_gene,)
DoHeatmap(ep_scaled,features =opc_related_gene)

FeaturePlot(
  object = ep,
  features = "MA0678.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

###############
library(MotifDb)
load("countsFromATAC.RData")
setwd("~/login123/scATAC-seq/EP/outs/binding_site")
library(rtracklayer)
olig2_bed <- import.bed("./GSM1040156_Olig2_ChIP-seq_OPC.bed")
sox10_bed <- import.bed("./GSM1869163_Sox10_d3_peaks.bed")
olig2_bed
sox10_bed

da_peaks_2 %>% dim()
data.frame(gene = (da_peaks_2$gene),stringsAsFactors = F)$gene %>% str_extract("chr(X|Y)|chr[0-9]*") -> chr
data.frame(gene = (da_peaks_2$gene),stringsAsFactors = F)$gene %>% str_remove("chr(X-|[0-9]*-|Y-)") ->irs
rm(start,starts)
gr <- GRanges(seqnames = Rle(chr),IRanges(irs),strand = "*")
gr
olig2_bed

olig2_peaks <- gr[gr %over% olig2_bed,]
sox10_peaks <- gr[gr %over% sox10_bed,]

setwd("~/ep_opc/")
export.bed(olig2_peaks,"olig2_binding_site.bed")
