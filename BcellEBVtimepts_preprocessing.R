#---set up packages
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(biomaRt)
library(celldex)
library(org.Hs.eg.db)
library(SingleR)
#---enter all count files and make into Seurat objects 
setwd('/home/em/Downloads/B_cell_EBV_tmpts_exp/SRR31408058_BcellEBV_d0_cts/outs/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data <- CreateSeuratObject(counts=matrix,
                           min.cells=1,
                           min.features=1,
                           project = 'day0')
setwd('/home/em/Downloads/B_cell_EBV_tmpts_exp/SRR31408054_BcellEBV_d2_cts/outs/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data2 <- CreateSeuratObject(counts=matrix,
                           min.cells=1,
                           min.features=1,
                           project = 'day2')
setwd('/home/em/Downloads/B_cell_EBV_tmpts_exp/SRR31408050_BcellEBV_d5_cts/outs/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data3 <- CreateSeuratObject(counts=matrix,
                            min.cells=1,
                            min.features=1,
                            project = 'day5')
setwd('/home/em/Downloads/B_cell_EBV_tmpts_exp/SRR31408046_BcellEBV_d8_cts/outs/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data4 <- CreateSeuratObject(counts=matrix,
                            min.cells=1,
                            min.features=1,
                            project = 'day8')
#---merge all seurat objects into one
data<-merge(data,
            y=c(data2,data3,data4))
dim(data)
table(Idents(data))
#---subset cells with 200 < features < 9000 and with 15% or less mitochondrial genes 
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt <15)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
#--normalize object
data <- NormalizeData(data)
#---find variable genes in object
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top1000 <- head(VariableFeatures(data), 1000)
head(top1000)
gc()
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
head(all.genes)
#---scale object
data <- ScaleData(data, features = all.genes)
#---run PCA
data <- RunPCA(data, features = VariableFeatures(object = data))
gc()
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
#---run UMAP
data <- RunUMAP(data, dims = 1:15)
data <-JoinLayers(data)
head(data@meta.data)
#--auto annotate cell subsets
surveyReferences()
searchReferences('human')
ref <- fetchReference("blueprint_encode", "2024-02-26")
ref
#--- convert seurat object into single cell experiment object 
rna.sce <- as.SingleCellExperiment(DietSeurat(data))
rna.sce
#---get general celltype annotations
ref.main <- SingleR(test = rna.sce,
                    assay.type.test = 1,
                    ref = ref,
                    labels = ref$label.main)
gc()
#---get fine celltype annotations
ref.fine <- SingleR(test = rna.sce,
                    assay.type.test = 1,
                    ref = ref,
                    labels = ref$label.fine)
table(ref.main$pruned.labels)
table(ref.fine$pruned.labels)
head(data@meta.data)
#---add general celltypes to metadata
data@meta.data$gen_celltype <- ref.main$pruned.labels
#---add fine celltypes to metadata
data@meta.data$fine_celltype <- ref.fine$pruned.labels
head(data@meta.data)
#---set object identities to fine celltype
data <- SetIdent(data, value = "fine_celltype")
table(Idents(data))
#---subset object by b-cell types
data<-subset(x = data, idents = c('naive B-cells','Class-switched memory B-cells','Memory B-cells','Plasma cells'),
             invert = FALSE)
table(Idents(data))
#---visualize final data
p1<-DimPlot(data,
            reduction = 'umap',
            label=TRUE,
            label.size =4,
            label.box = FALSE,
            raster=FALSE,
            pt.size = 0.5,
            seed=1,
            cols.highlight = c('grey'),
            dims = c(1,2),
            repel = TRUE,
            group.by = 'orig.ident')
p2<-DimPlot(data,
            reduction = 'umap',
            label=TRUE,
            label.size =4,
            label.box = FALSE,
            raster=FALSE,
            pt.size = 0.5,
            seed=1,
            cols.highlight = c('grey'),
            dims = c(1,2),
            repel = TRUE,
            group.by = 'fine_celltype')
p1+p2
#---write Seurat object to rda file 
SaveSeuratRds(data,
              file='/home/em/Downloads/B_cell_EBV_tmpts_exp/GSE282400_BcellEBV_viralpanel.rda')
