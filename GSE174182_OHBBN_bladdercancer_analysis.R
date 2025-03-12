library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(biomaRt)
library(celldex)
library(org.Mm.eg.db)
library(SingleR)
library(plotly)
library(htmlwidgets)
library(ReactomePA)
library(DOSE)
library(clusterProfiler)
data<-readRDS(file='GSE174182_BladderCancer_viralpanel.rds')
cat('GSE174182: Single-cell RNAseq and longitudinal proteomic analysis of a novel semi-spontaneous urothelial cancer model reveals tumor cell heterogeneity and pretumoral urine protein alterations')
cat('Public on Jun 09, 2021')
cat('Drug used: N-butyl-N-(4-hydroxybutyl) nitrosamine, a model compound that induces high-grade, invasive tumors in the urinary bladder.')
cat('Bladder cancer, one of the most prevalent malignancies worldwide, remains hard to classify due to a staggering molecular complexity. Despite a plethora of diagnostic tools and therapies, it is hard to outline the key steps leading up to the transition from high-risk non–muscle-invasive bladder cancer (NMIBC) to muscle-invasive bladder cancer (MIBC). Carcinogen-induced murine models can recapitulate urothelial carcinogenesis and natural anti-tumor immunity. Herein, we have developed and profiled a novel model of progressive NMIBC based on 10 weeks of OH-BBN exposure in hepatocyte growth factor/cyclin dependent kinase 4 (R24C) (Hgf-Cdk4R24C) mice. The profiling of the model was performed by histology grading, single cell transcriptomic and proteomic analysis, while the derivation of a tumorigenic cell line was validated and used to assess in vivo anti-tumor effects in response to immunotherapy. Established NMIBC was present in females at 10 weeks post OH-BBN exposure while neoplasia was not as advanced in male mice, however all mice progressed to MIBC. Single cell RNA sequencing analysis revealed an intratumoral heterogeneity also described in the human disease trajectory. Moreover, although immune activation biomarkers were elevated in urine during carcinogen exposure, anti-programmed cell death protein 1 (anti-PD1) monotherapy did not prevent tumor progression. Furthermore, anti-PD1 immunotherapy did not control the growth of subcutaneous tumors formed by the newly derived urothelial cancer cell line. However, treatment with CpG-oligodeoxynucleotides (ODN) significantly decreased tumor volume, but only in females. In conclusion, the molecular map of this novel preclinical model of bladder cancer provides an opportunity to further investigate pharmacological therapies ahead with regards to both targeted drugs and immunotherapies to improve the strategies of how we should tackle the heterogeneous tumor microenvironment in urothelial bladder cancer to improve responses rates in the clinic.')
cat('Experimental design: Single cell analysis was performed for whole mouse bladders after a) 10 weeks of OH-BBN exposure and after b) progression to MIBC, as well as of c) healthy Hgf-Cdk4R24C bladders. On the 10-week endpoint Hgf-Cdk4R24C tumors required pooling of samples due to low cell numbers per bladder (one pool of n= 4 females and one of n= 5 males).')
cat('Uppsala University Sweden')
cat('Kerzeli IK, Lord M, Doroszko M, Elgendy R et al. Single-cell RNAseq and longitudinal proteomic analysis of a novel semi-spontaneous urothelial cancer model reveals tumor cell heterogeneity and pretumoral urine protein alterations. PLoS One 2021;16(7):e0253178. PMID: 34232958')
data <- SetIdent(data, value = "orig.ident")
#-------------------------------------------------------------------------------
cat('Spleen focus-forming virus, complete genome')
cat('Wolff L, Ruscetti S. The spleen focus-forming virus (SFFV) envelope gene, when introduced into mice in the absence of other SFFV genes, induces acute erythroleukemia. J Virol. 1988 Jun;62(6):2158-63. doi: 10.1128/JVI.62.6.2158-2163.1988. PMID: 2835516; PMCID: PMC253318.')
table(Idents(data))
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-001500.1-viral`>0)
cat('SFFV infected cells increase after OH-BBN treatments') 
table(Idents(infected_cells))
#---Fig5a
VlnPlot(infected_cells,
        features = c('NC-001500.1-viral'),
        layer='data',
        group.by = 'orig.ident',
        pt.size = 3,
        cols = c('grey','red'))+ ggtitle('Spleen focus-forming virus (NC-001500.1)
logNormalized counts')+NoLegend()
FindMarkers(data,
            features=c('NC-001500.1-viral'),
            ident.1 = 'OH-BBN')
all.genes<-row.names(data)
#---differntially expressed genes using healthy as control
DEG<-FindMarkers(infected_cells,
                 ident.1 = 'OH-BBN',
                 ident.2 = 'healthy',
                 features=c(all.genes),
                 logfc.threshold = 0.5,
                 min.pct=0.4,
                 test.use='bimod',
                 only.pos=TRUE)
cat('Differentially expressed genes comparing treated to untreated')
head(DEG)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
geneid <- row.names(DEG)[1:1086]
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Mm.eg.db,
               ont = "BP")
#---Figure 5b
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        label_format=50)+ ggtitle('GO:Biological Process:Spleen focus-forming 
virus+ cells. Upregulated genes comparing 
OH-BBN treated and untreated. Count
indicate number of genes')
#-------------------------------------------------------------------------------
data<-readRDS(file='GSE174182_BladderCancer_viralpanel.rds')
data <- SetIdent(data, value = "orig.ident")
cat('Murine type C retrovirus, complete genome')
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-001702.1-viral`>0)
cat('Murine type C RV cells increase after OH-BBN treatments')
table(Idents(infected_cells))
#---Figure 5c
VlnPlot(infected_cells,
        features = c('NC-001702.1-viral'),
        layer = 'data',
        pt.size=3,
        cols = c('grey','red'))+ 
  ggtitle('Murine type C retrovirus, complete genome (NC_001702.1)
logNormalized counts')+NoLegend()
FindMarkers(data,
            features=c('NC-001702.1-viral'),
            ident.1 = 'OH-BBN',
            test.use='bimod')
all.genes<-row.names(data)
DEG<-FindMarkers(infected_cells,
                 ident.1 = 'OH-BBN',
                 ident.2 = 'healthy',
                 features=c(all.genes),
                 logfc.threshold = 0.7,
                 min.pct=0.5,
                 only.pos=TRUE,
                 test.use='bimod')
cat('Differentially expressed genes comparing OH-BBN treated and healthy cells')
head(DEG)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
geneid <- row.names(DEG)[1:900]
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
reactome <- enrichPathway(genes$entrezgene_id,
                          organism = 'mouse')
cat('Reactome pathway DEG annotations')
#---Fig5d
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Murine type C RV+ 
cells. Upregulated DEGs comparing 
OH-BBN treatments and untreated controls,
Counts indicate number of genes.',
        label_format=50)
#---------------------------------------------------------------------
data<-readRDS(file='GSE174182_BladderCancer_viralpanel.rds')
data <- SetIdent(data, value = "orig.ident")
cat('Mouse mammary tumor virus, complete genome')
table(Idents(data))
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-001503.1-viral`>0)
cat('Mouse mammary tumor virus cells increase after OH-BBN treatments')
table(Idents(infected_cells))
#---Fig5e
VlnPlot(infected_cells,
        features = c('NC-001503.1-viral'),
        layer = 'data',
        pt.size=3,
        cols = c('grey','red'))+ggtitle('Mouse mammary tumor virus complete 
genome (NC-001503.1) logNormalized counts')+NoLegend()
FindMarkers(data,
            features=c('NC-001503.1-viral'),
            ident.1 = 'OH-BBN',
            test.use='bimod')
all.genes<-row.names(data)
DEG<-FindMarkers(infected_cells,
                 ident.1 = 'OH-BBN',
                 ident.2 = 'healthy',
                 features=c(all.genes),
                 logfc.threshold = 0.5,
                 min.pct=0.4,
                 test.use='bimod',
                 only.pos=TRUE)
cat('Differentially expressed genes comparing OH-BBN treated and healthy')
head(DEG)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
geneid <- row.names(DEG)[1:292]
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Mm.eg.db,
               ont = "BP")
#---Figure 5f
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        label_format=50)+ ggtitle('GO:Biological Process:Mouse mammary 
tumor virus+ cells. Upregulated DEGs in OH-BBN 
treatment compared to untreated cells,
Count indicates number of genes.')
cat('Split into uninfected and infected cells')
data$viral.groups <- 'viral.pos'
data$viral.groups[WhichCells(data, expression= `NC-001503.1-viral` < 0.1)] <- 'viral.neg'
DimPlot(data, reduction = 'pca',split.by = 'viral.groups')
VlnPlot(data,
        features = c('NC-001503.1-viral'),
        group.by = 'viral.groups')+ggtitle('Mouse mammary tumor virus, complete genome:NC-001503.1')
head(data@meta.data)
data <- SetIdent(data, value = "viral.groups")
table(Idents(data))
DEG2<-FindMarkers(data,
                  ident.1 = 'viral.pos',
                  ident.2 = 'viral.neg',
                  logfc.threshold = 0.5,
                  min.pct = 0.5,
                  test.use = 'bimod',
                  only.pos=FALSE)
cat('Differentially expressed genes comparing infected and uninfected cells')
head(DEG2)
geneid <- row.names(DEG2)[1:643]
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Mm.eg.db,
               ont = "BP")
#---Figure 5g
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process:Mouse mammary
tumor virus+ compared to -negative cells.  
Count indicates number of genes',
        showCategory=12,
        label_format=50)


