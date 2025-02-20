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
data<-readRDS(file='<path_to>/B_cell_EBV_tmpts_methods/GSE282400_BcellEBV_viralpanel.rda')
cat('GSE282400: Sp140L Is a Herpesvirus Restriction Factor [scRNA-seq]')
cat('Herpesviruses, including the oncogenic Epstein-Barr Virus (EBV), must bypass host DNA sensing mechanisms to drive infection and pathogenesis. The first viral latency protein expressed, EBNA-LP, is essential for the transformation of naïve B cells, yet its role in evading host defenses remains unclear. Using single-cell RNA sequencing of EBNA-LP-Knockout (LPKO)-infected B cells, we reveal an antiviral response landscape implicating the ‘speckled proteins’ as key restriction factors countered by EBNA-LP. Specifically, loss of SP100 or the primate-specific SP140L reverses the restriction of LPKO, suppresses a subset of canonically interferon-stimulated genes, and restores viral gene transcription and cellular proliferation. Notably, we also identify Sp140L as a restriction target of the herpesvirus saimiri ORF3 protein, implying a role in immunity to other DNA viruses. This study reveals Sp140L as a restriction factor that we propose links sensing and transcriptional suppression of viral DNA to an IFN-independent innate immune response, likely relevant to all nuclear DNA viruses.')
cat('Public on Nov 20, 2024')
cat('Cable JM, Wongwiwat W, Grabowski JC, White RE et al. Sp140L Is a Novel Herpesvirus Restriction Factor. bioRxiv 2024 Dec 14. PMID: 39713285')
cat('Duke University. Molecular Genetics and Microbiology')
cat('Experimental design: Single cell RNAseq data from B cells in a timecourse of LPKO or WT EBV infection. Samples including uninfected (Day 0), and 2, 5, and 8 days post-infection.')

ebv1<-read.csv('<path_to>/B_cell_EBV_tmpts_methods/genes/ebv_genes.csv')
ebv1<-ebv1$x
ebv1<-FetchData(data,vars = c('ident',ebv1),layer = 'counts')
ebv1<-colnames(ebv1)[-1]
cat('All viral sequence reads detected')
ebv1

ebv2<-read.csv('<path_to>/B_cell_EBV_tmpts_methods/genes/ebv2_genes.csv')
ebv2<-ebv2$x
ebv2<-FetchData(data,vars = c('ident',ebv2),layer = 'counts')
ebv2<-colnames(ebv2)[-1]
cat('All viral sequence reads detected')
ebv2

hervk113<-read.csv('<path_to>/B_cell_EBV_tmpts_methods/genes/HERVK113_sequence.csv')
hervk113<-hervk113$x
hervk113

sv40<-read.csv('<path_to>/B_cell_EBV_tmpts_methods/genes/SV40_sequence.csv')
sv40<-sv40$x
sv40

stxPhagevB<-read.csv('<path_to>/B_cell_EBV_tmpts_methods/genes/stxPhage_sequence.csv')
stxPhagevB<-stxPhagevB$x
stxPhagevB
data <- SetIdent(data, value = "orig.ident")
#---Figure 1
DoHeatmap(data,
          features = c('PTPRC','CD19','MS4A1','CD22','CD274','PDCD1LG2',
                       'IL10','TNF','MYC','MKI67',ebv1,ebv2,hervk113,sv40,stxPhagevB),
          size = 3)
FindMarkers(data,
            features = c('PTPRC','CD19','MS4A1','CD22','CD274','PDCD1LG2',
                         'IL10','TNF','MYC','MKI67'),
            test.use='bimod',
            ident.1 = 'day8',
            ident.2 = 'day0')
#-------------------------------------------------------------------------------
cat('Human endogenous retrovirus K113 complete genome')
table(Idents(data))
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-022518.1-HERVK113`>0)
cat('HERVK113+ cells increase from day0 to day8')
table(Idents(infected_cells))
infected_subset <- SetIdent(infected_cells, value = "fine_celltype")
table(Idents(infected_subset))
infected_subset <- subset(infected_subset,
                          idents = c('naive B-cells','Class-switched memory B-cells',
                                     'Memory B-cells','Plasma cells'), invert = FALSE)
table(Idents(infected_subset))
#---Figure 4a
VlnPlot(infected_cells,
        features = c('NC-022518.1-HERVK113'),
        layer = 'data',
        pt.size=3,
        cols = c('red','orange','yellow','skyblue'))+ ggtitle('Human endogenous retrovirus K113 
        complete genome (NC-022518.1) logNormalized counts')+NoLegend()
FindMarkers(infected_cells,
            features=c('NC-022518.1-HERVK113'),
            ident.1 = 'day8',
            test.use='bimod')
all.genes<-row.names(data)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
cat('Split into uninfected and infected cells')
data$viral.groups <- 'viral.pos'
data$viral.groups[WhichCells(data, expression= `NC-022518.1-HERVK113` < 0.1)] <- 'viral.neg'
DimPlot(data, reduction = 'pca',split.by = 'viral.groups')
VlnPlot(data,
        features = c('NC-022518.1-HERVK113'),
        group.by = 'viral.groups')+ ggtitle('Human endogenous retrovirus K113 complete genome:NC-022518.1')
head(data@meta.data)
data <- SetIdent(data, value = "viral.groups")
table(Idents(data))
DEG2<-FindMarkers(data,
                  ident.1 = 'viral.pos',
                  ident.2 = 'viral.neg',
                  logfc.threshold = 0.4,
                  min.pct = 0.3,
                  test.use = 'bimod',
                  only.pos=FALSE)
cat('Differentially expressed genes comparing viral positive and negative subsets')
head(DEG2)
geneid <- row.names(DEG2)[1:776]
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
reactome <- enrichPathway(genes$entrezgene_id,
                          organism = 'human')
cat('Reactome pathway DEG annotations')
#---
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Human endogenous 
retrovirus K113+ positive comparison to 
negative- cells. Count axis indicates 
number of genes',
        label_format=50)
#-------------------------------------------------------------------------------
data<-readRDS(file='<path_to>/B_cell_EBV_tmpts_methods/GSE282400_BcellEBV_viralpanel.rda')
data <- SetIdent(data, value = "orig.ident")
cat('Simian virus 40 complete genome')
cat('SV40 is a potent DNA tumor virus that induces tumors in rodents and transforms many types of cells in culture, including those of human origin.
Reference: Janet S. Butel, John A. Lednicky, Cell and Molecular Biology of Simian Virus 40: Implications for Human Infections and Disease , JNCI: Journal of the National Cancer Institute, Volume 91, Issue 2, 20 January 1999, Pages 119–134, https://doi.org/10.1093/jnci/91.2.119')
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-001669.1-SV40`>0)
cat('Simian virus 40+ cells increase from day2 to day8')
table(Idents(infected_cells))
#--Fig4c
VlnPlot(infected_cells,
        features = c('NC-001669.1-SV40'),
        layer = 'data',
        pt.size=3,
        cols = c('yellow','orange','red'))+ ggtitle('Simian virus 40 complete genome (NC-001669.1)
logNormalized counts')+NoLegend()
FindMarkers(infected_cells,
            features=c('NC-001669.1-SV40'),
            ident.1 = 'day8',
            ident.2 = 'day2',
            test.use='bimod')
all.genes<-row.names(data)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
cat('Split into uninfected and infected cells')
data$viral.groups <- 'viral.pos'
data$viral.groups[WhichCells(data, expression= `NC-001669.1-SV40` < 0.1)] <- 'viral.neg'
DimPlot(data, reduction = 'pca',split.by = 'viral.groups')
VlnPlot(data,
        features = c('NC-001669.1-SV40'),
        group.by = 'viral.groups')+ ggtitle('Simian virus 40')
head(data@meta.data)
data <- SetIdent(data, value = "viral.groups")
table(Idents(data))
DEG2<-FindMarkers(data,
                  ident.1 = 'viral.pos',
                  ident.2 = 'viral.neg',
                  logfc.threshold = 0.4,
                  min.pct = 0.3,
                  test.use = 'bimod',
                  only.pos=FALSE)
cat('Differentially expressed genes comparing infected and uninfected cells')
head(DEG2)
geneid <- row.names(DEG2)
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
reactome <- enrichPathway(genes$entrezgene_id,
                          organism = 'human')
cat('Reactome pathway DEG annotations')
#---Fig4d
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Simian virus 40+ 
positive comparison to negative- cells.
Count axis indicates number of genes.',
        label_format=50)
#---------------------------------------------------------------------
data<-readRDS(file='<path_to>/B_cell_EBV_tmpts_methods/GSE282400_BcellEBV_viralpanel.rda')
data <- SetIdent(data, value = "orig.ident")
cat('Stx converting phage vB_EcoS_ST2-8624, complete genome')
cat('Viruses; Duplodnaviria; Heunggongvirae; Uroviricota;
Caudoviricetes; Sepvirinae; Traversvirus; Traversvirus ST28624.')
cat('Submitted by: Laboratory of Molecular Biology, Institute of Biochemistry and Biophysics, Polish Academy of Sciences, Wita Stwosza 59, Gdansk 80-308, Poland')
cat('Shiga toxins (Stx) comprise a family of potent cytotoxins that are involved in severe human disease. Stx are mainly produced by Escherichia coli isolated from human and nonhuman sources, and by Shigella dysenteriae type 1.
Reference: Schmidt H. Shiga-toxin-converting bacteriophages. Res Microbiol. 2001 Oct;152(8):687-95. doi: 10.1016/s0923-2508(01)01249-9. PMID: 11686382.')
table(Idents(data))
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-049922.1-stxPhage`>0)
cat('Stx converting phage vB+ cells increase from day8 to day2')
table(Idents(infected_cells))
#---Fig4e
VlnPlot(infected_cells,
        features = c('NC-049922.1-stxPhage'),
        layer = 'data',
        pt.size=3)+ggtitle('Stx converting phage vB_EcoS_ST2-8624, 
        complete genome (NC_049922.1). logNormalized counts')
FindMarkers(data,
            features=c('NC-049922.1-stxPhage'),
            ident.1 = 'day8',
            ident.2 = 'day2',
            test.use='bimod')
all.genes<-row.names(data)
cat('Split into uninfected and infected cells')
data$viral.groups <- 'viral.pos'
data$viral.groups[WhichCells(data, expression= `NC-049922.1-stxPhage` < 0.1)] <- 'viral.neg'
DimPlot(data, reduction = 'pca',split.by = 'viral.groups')
VlnPlot(data,
        features = c('NC-049922.1-stxPhage'),
        group.by = 'viral.groups')+ggtitle('Stx converting phage vB_EcoS_ST2-8624:NC-049922.1')
head(data@meta.data)
data <- SetIdent(data, value = "viral.groups")
table(Idents(data))
DEG2<-FindMarkers(data,
                  ident.1 = 'viral.pos',
                  ident.2 = 'viral.neg',
                  logfc.threshold = 0.5,
                  min.pct = 0.4,
                  test.use = 'bimod',
                  only.pos=FALSE)
cat('Differentially expressed genes comparing infected to uninfected cells')
head(DEG2)
geneid <- row.names(DEG2)
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
#---Fig4f
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process:Stx converting phage 
vB_EcoS_ST2-8624 Infected v non-infected genes',
        showCategory=6,
        label_format=50)
