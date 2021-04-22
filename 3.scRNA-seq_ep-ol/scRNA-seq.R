library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
##########################load the data, creat seurat objects################
dir = "./login123/scRNA-seq/olig2_inj/cellranger_counts"
samples <- paste0("SRR115109",36:41)
gse <- c("Foxj1_tdT_Uninj","Foxj1_tdT_2wpi","Foxj1_Olig2_tdT_Uninj","Foxj1_Olig2_tdT_2wpi","Foxj1_tdT_4wpi","Foxj1_Olig2_tdT_4wpi")
seqdata <- vector(mode = "list",length = length(samples))

for (id in seq_along(seqdata)){
  seqdata[[id]] <- Read10X(data.dir = file.path(dir,samples[id],"outs"))
}
seqdata

seurat_oj <- vector(mode = "list")
for (id in 1:6){
  index <- gse[[id]]
  seurat_oj[[index]] <- CreateSeuratObject(counts = seqdata[[id]], project = index, min.cells = 3, min.features = 200)
}
seurat_oj 
###########################################merge objects##########################

sc.all <- merge(seurat_oj$Foxj1_tdT_Uninj, y = c(seurat_oj$Foxj1_tdT_2wpi,seurat_oj$Foxj1_Olig2_tdT_Uninj,seurat_oj$Foxj1_Olig2_tdT_2wpi,seurat_oj$Foxj1_tdT_4wpi,seurat_oj$Foxj1_Olig2_tdT_4wpi), add.cell.ids = gse, project = "sc")
sc.all

head(colnames(sc.all))
tail(colnames(sc.all))
library(stringr)
colnames(sc.all) %>% str_trunc(20,side = "right") %>% unique()
table(sc.all$orig.ident)
sc.all
############################################QC
sc.all[["percent.mt"]] <- PercentageFeatureSet(sc.all, pattern = "^mt-")
head(sc.all@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(sc.all, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3,pt.size = 0)
sc.all <- subset(sc.all, subset = nFeature_RNA > 250 & nFeature_RNA < 5000 & percent.mt < 20)
#########normlization
sc.all <- NormalizeData(sc.all, normalization.method = "LogNormalize", scale.factor = 10000)
sc.all <- FindVariableFeatures(sc.all, selection.method = "vst", nfeatures = 8000)##########here author select 8050
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc.all), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc.all)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
########################scale the data
all.genes <- rownames(sc.all)
sc.all <- ScaleData(sc.all, features = all.genes)
###################Perform linear dimensional reduction
sc.all <- RunPCA(sc.all, features = VariableFeatures(object = sc.all))
# Examine and visualize PCA results a few different ways
print(sc.all[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sc.all, dims = 1:2, reduction = "pca")
DimPlot(sc.all, reduction = "pca")

######################################filter by tdtomato
match("tdtomato",all.genes)
length(all.genes)
sc.all@assays$RNA[12197,1:100] ==0 
########################################################
####plot pcs
DimHeatmap(sc.all, dims = 1, cells = 500, balanced = TRUE)
############## determine the PCs number
sc.all <- JackStraw(sc.all, num.replicate = 100)
sc.all <- ScoreJackStraw(sc.all, dims = 1:20)
JackStrawPlot(sc.all, dims = 1:15)
ElbowPlot(sc.all)


#################cluster
sc.all <- FindNeighbors(sc.all, dims = 1:10)
sc.all <- FindClusters(sc.all, resolution = 0.2)
head(Idents(sc.all), 5)
sc.all <- RunUMAP(sc.all, dims = 1:10)
DimPlot(sc.all, reduction = "umap",label = T)

############fiter cells
mg_gene <- c("Lyz1", "P2ry12", "Aif1", "Cx3cr1", "C1qa", "Fcrls")
FeaturePlot(sc.all,features = mg_gene,label = T)
fb_gene <- c("Col1a1", "Lum", "Dcn", "Fap")
FeaturePlot(sc.all,features = fb_gene,label = T)
FeaturePlot(sc.all,features = "tdtomato",label = T)
########################cluter 4 5 7 should be removed !!!
table(sc.all$RNA_snn_res.0.2 )
#                        0    1    2    3    4    5    6    7    8 
#                       1323  870  562  307  276  192  171   33   28 
sc <- sc.all[,! sc.all$RNA_snn_res.0.2 %in% c(5,7)]
sc
sc.all

#################cluster
sc <- FindNeighbors(sc, dims = 1:9)
sc <- FindClusters(sc, resolution = 0.25)
head(Idents(sc.all), 5)
sc <- RunUMAP(sc, dims = 1:9)
sc <- RunTSNE(sc,dims = 1:9)
DimPlot(sc, reduction = "umap",label = T)
DimPlot(sc, reduction = "tsne",label = T)
label_marker <- c("tdtomato","Foxj1","Sox9","Sox10")
FeaturePlot(sc,features = label_marker,label = T)


###############figI
figi <-c("Sox9","Cd27","Stmn2","Vim","Id4","Sox10","Pdgfra","Myrf","Mog")
VlnPlot(sc,features = figi,pt.size = 0)


#############
###4:213cells 7:28cells
sc6 <- sc[,!sc$RNA_snn_res.0.25 %in% c(4,7)]
table(sc6$RNA_snn_res.0.25)
############## rename
sc6 <- RunUMAP(sc6, dims = 1:9)


cluster.ids <- c("EP","AE1","AE2","AE3", "epOL1", "epOL2")
levels(sc6)
names(cluster.ids) <- levels(sc6)
sc6 <- RenameIdents(sc6, cluster.ids)

DimPlot(sc6, reduction = "umap",label = T)
figi <-c("Sox9","Cd27","Stmn2","Vim","Id4","Sox10","Pdgfra","Myrf","Mog")
VlnPlot(sc6,features = figi,pt.size = 0)
######################figj
"Ptgds
Mobp
Mag
Plp1
Ptprz1
Fyn
Cspg5
Gfap
Lgals3
Cryab
Nnat
Vim
Stmn2
Cd27
Rarres2" %>% str_split("\n") %>% unlist() %>% rev() ->gl1
gl1
DoHeatmap(sc6,features = gl1)
allmarkers <- FindAllMarkers(sc6,only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
topmarker <- allmarkers %>% group_by(cluster) %>% top_n(60,wt = avg_logFC)
DoHeatmap(sc6,features = topmarker$gene)
DimPlot(sc6[,sc6$class %in% c("Foxj1_Olig2_tdT_2wpi", "Foxj1_Olig2_tdT_4wpi","Foxj1_Olig2_tdT_Uninj")],reduction = "tsne")
#########################figK
colnames(sc6) %>% grep("gse")
rm(wpi2,wpi4,uninj)
rm(g1,g2,g3,g4,g5,g6)
sc6$class %>% unique()

sc6@meta.data$class <- colnames(sc6) %>% str_remove("[_ATCG]{17}") %>%　str_remove("-1")
d1 <- DimPlot(sc6, reduction = "umap",label = F, group.by = 'class',cols = c('grey', 'grey','red','grey','grey', 'black'))
d2 <- DimPlot(sc6, reduction = "umap",label = F, group.by = 'class',cols = c('red', 'grey','grey','black','grey', 'grey'))
d3 <- DimPlot(sc6, reduction = "umap",label = F, group.by = 'class',cols = c('grey', 'red','grey','grey','black', 'grey'))

d1 + d2 + d3 






















############################################################################ diffusion map analysis & Trajectory Analysis
install.packages("princurve")
library(princurve)
library(SingleCellExperiment) 
install("clusterExperiment")
library(clusterExperiment)
install("gam")
library(gam)
library(corrplot)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(knitr)
# Set folder location for saving output files. This is also the same location as input data.
mydir <- "~/ep_opc/trajectory/"
# Objects to save.
Rda.destiny.path <- paste0(mydir, "trajectory_destiny.Rda")
Rda.slingshot.path <- paste0(mydir, "trajectory_slingshot.Rda")
set.seed(1)  # set a seed for your random number generator to get reproducible results 
opts_chunk$set(fig.align = "center")
##############convert seurat obj to sc!!
sdata <- as.SingleCellExperiment(sc6)
structure(sdata)
sdata$class %>% table()
sdata <- sdata[,sdata$class %in% c("Foxj1_Olig2_tdT_2wpi", "Foxj1_Olig2_tdT_4wpi","Foxj1_Olig2_tdT_Uninj")]
############ the cell type counts
table(sdata$ident)
table(sdata$class)
# Run PCA . Use the runPCA function from the SingleCellExperiment package.
#sdata <- runPCA(sdata, ncomponents = 50)
pca <- reducedDim(sdata, "PCA")
head(pca)
dim(pca)
sdata$PC1 <- pca[, 1]
sdata$PC2 <- pca[, 2]
head(colData(sdata))

# Plot PC biplot 
# colData(sdata) accesses the cell metadata DataFrame object for sdata.
# Look at Figure 1A of the paper as a comparison to your PC biplot.
ggplot(as.data.frame(colData(sdata)), aes(x = PC1, y = PC2, color = ident)) + geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")

# PCA is a simple approach and can be good to compare to more complex algorithms 
# designed to capture differentiation processes. As a simple measure of pseudotime 
# we can use the coordinates of PC1.
sdata$pseudotime_PC1 <- rank(sdata$PC1)  # rank cells by their PC1 score

ggplot(as.data.frame(colData(sdata)), aes(x = pseudotime_PC1, y = ident,colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")
ggsave(paste0(mydir, "/pseudotime_PC1.png"))
# Try separating the cell types using other PCs. How does the separation look?
sdata$PC5 <- pca[, 5]
sdata$PC6 <- pca[, 6]
ggplot(as.data.frame(colData(sdata)), aes(x = PC5, y = PC6, color = ident)) + geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("PC5") + ylab("PC6") + ggtitle("PC biplot")

#14.4 Diffusion map pseudotime
#Let us see how a more advance trajectory inference method, diffusion maps and diffusion pseudotime, 
#performs at placing cells along the expected differentiation trajectory.
#Diffusion maps were introduced by Ronald Coifman and Stephane Lafon[http://www.sciencedirect.com/science/article/pii/S1063520306000546],
#and the underlying idea is to assume that the data are samples from a diffusion process.
#The method infers the low-dimensional manifold by estimating the eigenvalues and eigenvectors for the diffusion operator related to the data.
#Angerer et al have applied the diffusion maps concept to the analysis of single-cell RNA-seq data to create an R package called destiny.
#We will use two forms of pseudotime: the first diffusion component and the diffusion pseudotime.

#  Prepare a counts matrix with labeled rows and columns. 
#tmp <- logcounts(sdata)  # access log-transformed counts matrix
#assays(sdata)
#logdata <- assay(sdata,"logcounts")
###colnames(logdata) <- cellLabels


# Make a diffusion map.
library(destiny)
dm <- DiffusionMap(sdata)
print(dm)
show(dm)
# Optional: Try different sigma values when making diffusion map.
# dm <- DiffusionMap(t(logdata), sigma = "local")  # use local option to set sigma
# sigmas <- find_sigmas(t(logdata), verbose = FALSE)  # find optimal sigma
# dm <- DiffusionMap(t(logdata), sigma = optimal_sigma(sigmas))  

# Plot diffusion component 1 vs diffusion component 2 (DC1 vs DC2). 
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Timepoint = cellLabels)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() + scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

################
deng <- logcounts(sdata)  # access log-transformed counts matrix
cellLabels <- colnames(sdata)
colnames(deng) <- cellLabels
dim(deng)
deng <- deng[,1:1055] %>% as.matrix() 


##############################*****************************************************************************
# Make a diffusion map.
# Diusion maps are based on a distance metric (difusion distance) which is conceptually relevant to how dierentiating cells follow noisy difusion-like dynamics, moving from a pluripotent state
#towards more differentiated states.The R package destiny implements the formulation of diusion maps presented inHaghverdi et al.(2015)
#which is especially suited for analyzing single-cell gene expression data from time-course experiments. It
#implicitly arranges cells along their developmental path, with bifurcations where dierentiation events occur.
sdata
assay(sdata,"logcounts")
dm <- DiffusionMap(sdata[rownames(sdata) %in% VariableFeatures(sc.all),])
plot(dm)
plot(dm, 1:3,col_by = 'ident',
     legend_main = 'Cell stage')
plot(dm, pch=20,
     col_by = 'ident',
     legend_main = 'Cell stage')

tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  Timepoint = sdata$ident)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point(size= 3) + scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()
############################################################################################################
scdata2<-sc6[,sc6$class %in% c("Foxj1_Olig2_tdT_2wpi", "Foxj1_Olig2_tdT_4wpi","Foxj1_Olig2_tdT_Uninj")]
scdata2_dm <- CreateDimReducObject(
  embeddings = dm@eigenvectors,
  key = "DM_",
  assay = "RNA"
)
scdata2[["dm"]]  <-scdata2_dm
p1 <- DimPlot(scdata2, reduction = "dm", pt.size = 1)
p2 <- DimPlot(scdata2, reduction = "dm",group.by = 'class',cols = c("grey","grey","red"),pt.size = 0.5)
p3 <- DimPlot(scdata2, reduction = "dm",group.by = 'class',cols = c("red","grey","grey"),pt.size = 0.5)
p4 <- DimPlot(scdata2, reduction = "dm",group.by = 'class',cols = c("grey","red","grey"),pt.size = 0.5)

p1 + p2 + p3 + p4


# Next, let us use the first diffusion component (DC1) as a measure of pseudotime.
# How does the separation by cell stage look?
sdata$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])    # rank cells by their dpt
ggplot(as.data.frame(colData(sdata)), 
       aes(x = pseudotime_diffusionmap, 
           y = ident, colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Diffusion component 1 (DC1)") + ylab("Timepoint") +
  ggtitle("Cells ordered by DC1")

# how many PCs to use in downstream applications like clustering.
plot(eigenvalues(dm), ylim = 0:1, pch = 20, xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')

sdata2 <- slingshot(sdata, reducedDim = 'PCA')
tmp <- data.frame(DC1 = reducedDims(sdata2)$PCA[, 1],
                  DC2 = reducedDims(sdata2)$PCA[, 2],
                  DC3 = reducedDims(sdata2)$PCA[, 3],
                  Timepoint = sdata2$ident)


ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
       geom_point(size= 3) + scale_color_tableau() + 
       xlab("Diffusion component 1") + 
       ylab("Diffusion component 2") +
       theme_classic()
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sdata2)$PCA, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)
