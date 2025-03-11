library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(biomaRt)
library(celldex)
library(org.Hs.eg.db)
library(SingleR)
library(plotly)
library(htmlwidgets)
library(ReactomePA)
library(DOSE)
library(clusterProfiler)
data<-readRDS(file='GBM_R1-4.VIR.rds')
data <- SetIdent(data, value = "orig.ident")
cat('Human endogenous retrovirus K113 complete genome')
data <- RenameIdents(data,
                     `gbm_peripheralR1` = 'peripheral')
data <- RenameIdents(data,
                     `gbm_peripheralR2` = 'peripheral')
data <- RenameIdents(data,
                     `gbm_coreR1` = 'core')
data <- RenameIdents(data,
                     `gbm_coreR2` = 'core')
table(Idents(data))
head(data@meta.data)
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-022518.1-viral`>0)
table(Idents(infected_cells))
#---Figure 3a
VlnPlot(infected_cells,
        features = c('NC-022518.1-viral'),
        layer = 'data',
        pt.size=3,
        cols = c('red','grey'))+ ggtitle('Human endogenous retrovirus K113 
        complete genome (NC-022518.1) logNormalized counts')+NoLegend()
FindMarkers(infected_cells,
            features=c('NC-022518.1-viral'),
            ident.1 = 'core',
            test.use='bimod')
#---Figure 3b
VlnPlot(infected_cells,
        features = c('NC-022518.1-viral'),
        layer='data',
        group.by = 'gen_celltype',
        pt.size = 3)+ ggtitle('Human endogenous retrovirus K113 
        complete genome (NC-022518.1) logNormalized counts')+NoLegend()
all.genes<-row.names(data)
DEG<-FindMarkers(infected_cells,
                 ident.1 = 'core',
                 ident.2 = 'peripheral',
                 features=c(all.genes),
                 logfc.threshold = 1,
                 min.pct=0.4,
                 only.pos=FALSE,
                 test.use='bimod')
head(DEG)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
geneid <- row.names(DEG)
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
reactome <- enrichPathway(genes$entrezgene_id,
                          organism = 'human')
cat('Reactome pathway DEG annotations')
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Human endogenous retrovirus 
K113+ cells comparing core to 
peripheral cells',
        label_format=50)
