library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(biomaRt)
library(celldex)
library(dplyr)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)
library(clusterProfiler)
library(SingleR)
library(plotly)
library(htmlwidgets)
data<-readRDS('GSE270109_glioma_RVk113.rds')
dim(data)
table(Idents(data))
dim(data)
genes<-'NC-022518.1'
data <- SetIdent(data, value = "orig.ident")
table(Idents(data))
DoHeatmap(data,
          features = c('PTPRC','GFAP','OLIG1','OLIG2','IFNG','TNF',
                       'CXCL8',genes),
          size = 3)
cat('Human endogenous retrovirus K113 complete genome: NC-022518.1')
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-022518.1`>0)
table(Idents(infected_cells))
#---Figure 2a
VlnPlot(infected_cells,
        features = c(genes),
        layer='data',
        group.by = 'orig.ident',
        pt.size = 3,
        cols = c('grey','red'))+ggtitle('Human endogenous retrovirus K113 complete 
genome (NC-022518.1) logNormalized counts')+NoLegend()
FindMarkers(infected_cells,
            features=c('NC-022518.1'),
            ident.1 = 'glioma_tumor',
           test.use='bimod')
#---Figure 2b
VlnPlot(infected_cells,
        features = c('NC-022518.1'),
        layer='data',
        group.by = 'gen_celltype',
        pt.size = 3,
        cols = c())+ggtitle('Human endogenous retrovirus K113 complete 
genome (NC-022518.1) logNormalized counts')+NoLegend()
all.genes<-row.names(data)
DEG<-FindMarkers(infected_cells,
                 ident.1 = 'glioma_tumor',
                 ident.2 = 'glioma_normal',
                 features=c(all.genes),
                 logfc.threshold = 1,
                 min.pct=0.4,
                 only.pos=TRUE)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
geneid <- row.names(DEG)[1:502]
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
#---Figure 2c
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        label_format=50)+ggtitle('GO: Biological Process HERVK113+ cells 
comparing glioma to normal controls')
