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
#---enter in B cell EBV day 8 post infection count file
setwd('<path_to>/B_cell_EBV_tmpts_exp/SRR31408046_BcellEBV_d8_viral_cts/outs/filtered_feature_bc_matrix')
barcodes_path <- 'barcodes.tsv.gz'
features_path <- 'features.tsv.gz'
matrix_path <- 'matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data <- CreateSeuratObject(counts=matrix,
                           min.cells=1,
                           min.features=1,
                           project = 'viral_seq_tester')
data <- NormalizeData(data)
all.genes <- rownames(data)
head(all.genes)
data <- ScaleData(data, features = all.genes)
#---get list of viral sequences
genes<-read.csv('<path_to>/B_cell_EBV_tmpts_exp/genes/viral_sequences.csv')
genes<-genes$x
genes
#---isolate viral sequences in dataset
genes<-FetchData(data,vars = c('ident',genes),layer = 'counts')
genes<-colnames(genes)[-1]
genes
DoHeatmap(data,
          features = c(genes))
#---Scan viral sequences in chunks of 10 for BLAST validations
#---for candidate viral sequences, use the samtools toolkit to view 
#---sequences by RefSeq id numbers
VlnPlot(data,
        features = c(genes[1:10]),
        layer='counts',
        pt.size = 2)
genes[1:10]
VlnPlot(data,
        features = c(genes[11:20]),
        layer='counts',
        pt.size = 2)
genes[11:20]
VlnPlot(data,
        features = c(genes[21:30]),
        layer='counts',
        pt.size = 2)
genes[21:30]
VlnPlot(data,
        features = c(genes[31:40]),
        layer='counts',
        pt.size = 2)
genes[31:40]
VlnPlot(data,
        features = c(genes[41:50]),
        layer='counts',
        pt.size = 2)
genes[41:50]
VlnPlot(data,
        features = c(genes[51:60]),
        layer='counts',
        pt.size = 2)
genes[51:60]
VlnPlot(data,
        features = c(genes[61:70]),
        layer='counts',
        pt.size = 2)
genes[61:70]
VlnPlot(data,
        features = c(genes[71:80]),
        layer='counts',
        pt.size = 2)
genes[71:80]
VlnPlot(data,
        features = c(genes[81:90]),
        layer='counts',
        pt.size = 2)
genes[81:90]
VlnPlot(data,
        features = c(genes[91:100]),
        layer='counts',
        pt.size = 2)
genes[91:100]
VlnPlot(data,
        features = c(genes[100:109]),
        layer='counts',
        pt.size = 2)
genes[100:109]
#---viral BLAST validated results
VlnPlot(data,features = c('NC-022518.1-viral'),
        layer = 'count',
        pt.size=4)+ggtitle(('Human endogenous retrovirus K113 complete genome'))
VlnPlot(data,features = c('NC-001669.1-viral'),
        layer = 'count',
        pt.size=4)+ggtitle(('Simian virus 40 complete genome'))
VlnPlot(data,features = c('NC-049922.1-viral'),
        layer = 'count',
        pt.size=4)+ggtitle(('Stx converting phage vB_EcoS_ST2-8624, complete genome'))

