library(monocle3)
cds <- readRDS("cds_3files.rds")
# cds$gene_matrix    cds$cell_metadata  cds$gene_metadata
# Make the CDS object
cds$gene_metadata$gene_short_name <- cds$gene_metadata$SYMBOL
cds <- new_cell_data_set(cds$gene_matrix,cell_metadata = cds$cell_metadata,gene_metadata = cds$gene_matrix)
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
plot_cells(cds)
#plot_cells(cds, color_cells_by="cao_cell_type")
cds <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE")
cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)


######################################## loom object
install.packages("hdf5r", configure.args="--with-hdf5=/ifs1/User/shangyao/miniconda3/bin/h5cc")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
library(loomR)
ls("package:loomR")
lfile <- connect(filename = "./downlaod/l1_spinalcord.loom", mode = "r+")
lfile

lfile[["matrix"]]
lfile[["col_attrs"]]
lfile[["layers"]]
lfile$row.attrs$Gene
lfile[["matrix"]][1:5, 1:5]

full.matrix <- lfile$matrix[, ]
dim(x = full.matrix)
lfile[["row_attrs"]]

lfile[["col_attrs/CellID"]]
lfile[["row_attrs/Gene"]]

gene.names <- lfile[["row_attrs/Gene"]][]
head(x = gene.names)
###gene.names not unique!!
index <- match(unique(gene.names),gene.names)

uniq_gene.names <- gene.names[index]
uniq_matrix <- lfile[["matrix"]][,index]
dim(uniq_matrix)


cell.id <- lfile[["col_attrs/CellID"]][]
index2 <- match(unique(cell.id),cell.id)
uniq_cell.id <- cell.id[index2]



class <- lfile[["col_attrs/Class"]][]
uniq_class <- class[index2]



lfile[["row_attrs/Gene"]]$dims
lfile[["row_attrs/Gene"]]$dims == lfile[["matrix"]]$dims[2]
lfile[["row_attrs/Gene"]]$dims == lfile$shape[1]
#data.subset <- lfile[["matrix"]][1:5, ]
#dim(x = data.subset)
uniq_matrix <- uniq_matrix[,index2]

data.Top2a <- lfile[["matrix"]][, lfile$row.attrs$Gene[] == "Top2a"]
sum(x = data.Top2a)

library(Seurat)
library(SingleCellExperiment)
#CreateSeuratObject
uniq_matrix <- t(x=uniq_matrix)
dim(uniq_matrix)
colnames(uniq_matrix) <- uniq_cell.id
rownames(uniq_matrix) <- uniq_gene.names
sc <- CreateSeuratObject(counts = uniq_matrix, project = "sc", min.cells = 15, min.features = 500,meta.data =data.frame(class = uniq_class,row.names = colnames(uniq_matrix)) )
########################
sc <- sc[,startsWith(colnames(sc),prefix = "10X53")]
sc[["sc.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA","sc.mt"), ncol = 3)
FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & sc.mt < 5)

sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc  <- FindVariableFeatures(sc , selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sc), 10)
top10
plot1 <- VariableFeaturePlot(sc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
sc <- RunPCA(sc, features = VariableFeatures(object =sc))
print(sc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sc, dims = 1:2, reduction = "pca")
DimPlot(sc, reduction = "pca")
sc <- JackStraw(sc, num.replicate = 100)
sc <- ScoreJackStraw(sc, dims = 1:10)
ElbowPlot(sc)

####remove neuron
glia <- sc[,sc@meta.data$class != "Neurons"]
glia <- glia[,glia@meta.data$class != "Blood"]

glia <- NormalizeData(glia, normalization.method = "LogNormalize", scale.factor = 10000)
glia  <- FindVariableFeatures(glia , selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(glia), 10)
top10
plot1 <- VariableFeaturePlot(glia)

glia <- ScaleData(glia, features = all.genes)
glia <- RunPCA(glia, features = VariableFeatures(object =glia))

glia <- JackStraw(sc, num.replicate = 100)
glia <- ScoreJackStraw(glia, dims = 1:10)
ElbowPlot(glia)
print(glia[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(glia, dims = 1:2, reduction = "pca")
DimPlot(glia, reduction = "pca")


glia
glia <- FindNeighbors(glia, dims = 1:10)
glia <- FindClusters(glia, resolution = 0.2)

glia <- RunUMAP(glia, dims = 1:10)
DimPlot(glia, reduction = "umap",label = T)


### find all markers

markers <- FindAllMarkers(glia, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
m1 <-  markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
m1
cluster1.markers <- FindMarkers(glia, ident.1 = 1, min.pct = 0.25)
top_n(cluster1.markers,n = 100, wt = avg_logFC)
cluster6.markers <- FindMarkers(glia, ident.1 = 6, min.pct = 0.25)
top_n(cluster6.markers,n = 100, wt = avg_logFC)
cluster13.markers <- FindMarkers(glia, ident.1 = 13, min.pct = 0.25)
top_n(cluster13.markers,n = 100, wt = avg_logFC)

FeaturePlot(glia,features = c("Acot5","Cd24a","Capsl","Gfap","Slc7a10","Gpr17","Egfl7","Foxd1"))
FeaturePlot(glia,features = "Hba-a2")


ep.marker <- c("1700012B09Rik", "Cfap126", "Fam183b","Dynlrb2", "1500015O10Rik")
acnt1.marker <- c("Slc7a10", "Agt", "Slc6a11","Fgfr3", "Cldn10")
atnt2.marker <- c("Igsf1", "Agt", "Slc6a11","Itih3")
opc.marker  <-c("Pdgfra","C1ql1","Sapcd2","Emid1","Neu4")
cop.marker <- c("Bmp4", "Neu4", "Gpr17","Pak4", "Lims2")
nfol.marker <- c("Rras2", "Tmem2","Il23a", "Tmem163")
mfol.marker <- c("2210011C24Rik", "Opalin", "Birc2","Tmem141", "Fam214a")
mol.marker <- c("Opalin", "Mal", "Mog", "Ninj2", "Klk6", "Anln","Cdkn1c", "Gjb1", "Prr5l", "Ugt8a", "Nkx2-9")
vec.marker <- c("Ly6c1", "Cldn5", "Ly6a","Pltp", "Pglyrp1")
per.marker <- c("Higd1b", "Ndufa4l2", "Rgs5","Kcnj8", "Flt1")



###############################
gene_list <- list(ep=ep.marker,acnt1 = acnt1.marker,atnt2 = atnt2.marker,opc=opc.marker,cop=cop.marker,nfol=nfol.marker,mfol = mfol.marker,mol=mol.marker,vec= vec.marker,per = per.marker)
object <- AddModuleScore(object = glia, features = gene_list$ep,name = "ep")
FeaturePlot(object = object, features = "ep1")
object <- AddModuleScore(object = glia, features = gene_list$acnt1,name = "acnta")
FeaturePlot(object = object, features = "acnta1")
object <- AddModuleScore(object = glia, features = gene_list$atnt2,name = "acntb")
FeaturePlot(object = object, features = "acntb1")
object <- AddModuleScore(object = glia, features = opc.marker,name = "opc")
FeaturePlot(object = object, features = "opc1")
object <- AddModuleScore(object = glia, features = cop.marker,name = "cop")
FeaturePlot(object = object, features = "cop1")
object <- AddModuleScore(object = glia, features = nfol.marker,name = "nfol")
FeaturePlot(object = object, features = "nfol1")
object <- AddModuleScore(object = glia, features = mfol.marker,name = "mfol")
FeaturePlot(object = object, features = "mfol1")
object <- AddModuleScore(object = glia, features = mol.marker ,name = "mol")
FeaturePlot(object = object, features = "mol1")
object <- AddModuleScore(object = glia, features = vec.marker ,name = "vec")
FeaturePlot(object = object, features = "vec1")
object <- AddModuleScore(object = glia, features = per.marker ,name = "per")
FeaturePlot(object = object, features = "per1")

#################################################
cluster.ids <- c("EP", "ACNT1", "ACNT2", "OPC", "COP/NFOL", "NFOL/MFOL", "MOL","VEC","PER","MG")
new.cluster.ids<-c("MOL","MOL","VEC","NFOL/MFOL","MOL","COP/NFOL","MG","ACNT1","ACNT2","PER", "MOL","EP","OPC","13")
names(new.cluster.ids) <- levels(glia)
glia <- RenameIdents(glia, new.cluster.ids)


DimPlot(glia, reduction = "umap", label = TRUE, pt.size = 0.1)+ NoLegend()
glia2 <- glia[,glia$RNA_snn_res.0.2 != "12"]
glia2 <- FindNeighbors(glia2, dims = 1:10)

glia2 <- FindClusters(glia2, resolution = 0.3)
glia2 <- RunUMAP(glia2, dims = 1:10)
DimPlot(glia2, reduction = "umap",label = F)
DimPlot(glia2, reduction = "umap", label = TRUE, pt.size = 0.1)+ NoLegend()

# Reorder identity classes
levels(x = glia2)
levels(x = glia2) <- cluster.ids
levels(x = glia2)
DimPlot(glia2, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend()


##################################


####################
FeaturePlot(glia2,features = c("Acot5","Cd24a","Capsl","Gfap","Slc7a10","Gpr17","Egfl7","Foxd1"))
markers <- FindAllMarkers(glia2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
m1 <-  markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
m1
clusterAC.markers <- FindMarkers(glia, ident.1 = "ACNT2",ident.2 = "ACNT1", min.pct = 0.25)
top_n(clusterAC.markers,n=20,wt=avg_logFC)

gene_plot <- c("Fam183b",'Slc6a11',"Gfap","Pdgfra","Gpr17","Enpp6","Birc2","Mog","Pltp",'Higd1b',"Tmem119")
FeaturePlot(glia2,features = gene_plot,label = T,label.size=3,repel = T) 
library(ggplot2)
ggsave("cluster11.tiff",dpi = 300,units = "cm",width = 60,height = 30)

FeaturePlot(glia2,features = "Sox10",label = F,repel = T,pt.size = 0.5) 





