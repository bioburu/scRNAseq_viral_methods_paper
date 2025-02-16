#--set up packages
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(biomaRt)
library(celldex)
library(org.Hs.eg.db)
library(SingleR)
library(ReactomePA)
library(DOSE)
library(clusterProfiler)
#---read in rda file
data<-readRDS(file='<path_to>/B_cell_EBV_tmpts_exp/GSE282400_BcellEBV_viralpanel.rda')
cat('GSE282400: Sp140L Is a Herpesvirus Restriction Factor [scRNA-seq]')
cat('Herpesviruses, including the oncogenic Epstein-Barr Virus (EBV), must bypass host DNA sensing mechanisms to drive infection and pathogenesis. The first viral latency protein expressed, EBNA-LP, is essential for the transformation of naïve B cells, yet its role in evading host defenses remains unclear. Using single-cell RNA sequencing of EBNA-LP-Knockout (LPKO)-infected B cells, we reveal an antiviral response landscape implicating the ‘speckled proteins’ as key restriction factors countered by EBNA-LP. Specifically, loss of SP100 or the primate-specific SP140L reverses the restriction of LPKO, suppresses a subset of canonically interferon-stimulated genes, and restores viral gene transcription and cellular proliferation. Notably, we also identify Sp140L as a restriction target of the herpesvirus saimiri ORF3 protein, implying a role in immunity to other DNA viruses. This study reveals Sp140L as a restriction factor that we propose links sensing and transcriptional suppression of viral DNA to an IFN-independent innate immune response, likely relevant to all nuclear DNA viruses.')
cat('Public on Nov 20, 2024')
cat('Cable JM, Wongwiwat W, Grabowski JC, White RE et al. Sp140L Is a Novel Herpesvirus Restriction Factor. bioRxiv 2024 Dec 14. PMID: 39713285')
cat('Duke University. Molecular Genetics and Microbiology')
cat('Experimental design: Single cell RNAseq data from B cells in a timecourse of LPKO or WT EBV infection. Samples including uninfected (Day 0), and 2, 5, and 8 days post-infection.')
#---get EBV type I CDS sequences
ebv<-read.csv('<path_to>/B_cell_EBV_tmpts_exp/genes/ebv_genes.csv')
ebv<-ebv$x
ebv<-FetchData(data,vars = c('ident',ebv),layer = 'counts')
ebv<-colnames(ebv)[-1]
cat('All ebv sequence reads detected')
ebv
data <- SetIdent(data, value = "orig.ident")
#---get EBV type II CDS sequences
ebv2<-read.csv('<path_to>/B_cell_EBV_tmpts_exp/genes/ebv2_genes.csv')
ebv2<-ebv2$x
ebv2<-FetchData(data,vars = c('ident',ebv2),layer = 'counts')
ebv2<-colnames(ebv2)[-1]
cat('All ebv2 sequence reads detected')
ebv2
#---get HERVK113 genome sequence
hervK113<-read.csv('<path_to>/B_cell_EBV_tmpts_exp/genes/HERVK113_sequence.csv')
hervK113<-hervK113$x
hervK113
#---get SV40 genome sequence
sv40<-read.csv('<path_to>/B_cell_EBV_tmpts_exp/genes/SV40_sequence.csv')
sv40<-sv40$x
sv40
#---get stxPhage genome sequence
stxPhage<-read.csv('<path_to>/B_cell_EBV_tmpts_exp/genes/stxPhage_sequence.csv')
stxPhage<-stxPhage$x
stxPhage
DoHeatmap(data,
          features = c('PTPRC','CD19','MS4A1','CD22','CD274','PDCD1LG2',
                       'IL10','TNF','MYC','MKI67',
                       ebv,ebv2,hervK113,sv40,stxPhage),
          size = 3)
#-------------------------------------------------------------------------------
cat('Human endogenous retrovirus K113 complete genome')
VlnPlot(data,
        features = c(hervK113),
        layer = 'data',
        pt.size=3)+ ggtitle('Human endogenous retrovirus K113 complete genome:NC-022518.1')
table(Idents(data))
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-022518.1-HERVK113` >0)
cat('HERVK113+ cell counts increase from day0 to day8')
table(Idents(infected_cells))
cat('HERVK113 reads in B-cell subsets')
VlnPlot(infected_cells,
        features = c(hervK113),
        layer='data',
        group.by = 'fine_celltype',
        pt.size = 3)+ ggtitle('Human endogenous retrovirus K113 complete 
                  genome:NC-022518.1')
cat('HERVK113 reads at timepoints')
VlnPlot(infected_cells,
        features = c(hervK113),
        layer='data',
        group.by = 'orig.ident',
        pt.size = 3)+ ggtitle('Human endogenous retrovirus K113 complete 
                  genome:NC-022518.1')
all.genes<-row.names(data)
cat('Differentially expressed genes comparing day8 to day0')
DEG<-FindMarkers(infected_cells,
                 ident.1 = 'day8',
                 ident.2 = 'day0',
                 features=c(all.genes),
                 logfc.threshold = 0.6,
                 min.pct=0.4,
                 only.pos=FALSE)
head(DEG)
cat('genes for DAVID Bioinformatics functional analysis')
cat(row.names(DEG))
DoHeatmap(infected_cells,
          features = c(row.names(DEG)),
          size = 4)+ ggtitle('Human endogenous retrovirus K113 infected cells only
                             ')
cat('get entrezgene ids for functional analysis')
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
geneid <- row.names(DEG)
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
cat('Reactome pathway DEG annotations')
reactome <- enrichPathway(genes$entrezgene_id,
                          organism = 'human')
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Human endogenous retrovirus 
K113 infected cells only. DEGs at day8 v day0',
        label_format=50)
edox <- setReadable(reactome, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox,
         foldChange=NULL,
         node_label='all',
         colorEdge=TRUE,
         cex_label_gene=0.1,
         circular=TRUE,
         cex_gene=0.1,
         showCategory = 4,
         node_label_size=0.1,
         cex_label_category=1)+ ggtitle(' Reactome pathways:Human endogenous retrovirus K113 infected cells
 DEGs at day8 v day0')
cat('Biological Process DEG annotations')
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        label_format=50)+ ggtitle(' GO:Biological Process:Human endogenous retrovirus 
 K113 infected cells only, DEGs at day8 v day0')
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
cat('Differentially expressed genes comparing viral positive and negative subsets')
DEG2<-FindMarkers(data,
                  ident.1 = 'viral.pos',
                  ident.2 = 'viral.neg',
                  logfc.threshold = 0.4,
                  min.pct = 0.3,
                  test.use = 'bimod',
                  only.pos=FALSE)
head(DEG2)
cat('Even out subsets')
data_subset<-subset(data,
                    downsample=71)
table(Idents(data_subset))
DoHeatmap(data_subset,
          features = c(row.names(DEG2)),
          size = 4)+ ggtitle('Human endogenous retrovirus K113 infected v non-infected cells
                             ')
cat('genes for DAVID Bioinformatics functional analysis')
cat(row.names(DEG2))
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
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Human endogenous retrovirus 
K113 Infected v non-infected cells',
        label_format=50)
cat('Loss-of-function mutations in the BRCA1 and BRCA2 genes increase the risk of cancer.
Reference: Zámborszky, J., Szikriszt, B., Gervai, J. et al. Loss of BRCA1 or BRCA2 markedly increases the rate of base substitution mutagenesis and has distinct effects on genomic deletions. Oncogene 36, 746–755 (2017). https://doi.org/10.1038/onc.2016.243')
edox <- setReadable(reactome, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox,
         foldChange=NULL,
         node_label='all',
         colorEdge=TRUE,
         cex_label_gene=0.5,
         circular=TRUE,
         cex_gene=0.1,
         showCategory = 5,
         node_label_size=1,
         cex_label_category=1)+ ggtitle(' Reactome pathways:Human endogenous retrovirus K113
 Infected v non-infected cells')
cat('Biological Process DEG annotations')
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        showCategory=10,
        label_format=50)+ ggtitle(' GO:Biological Process:Human endogenous retrovirus 
 K113 Infected v non-infected cells')
p1<-FeaturePlot(data_subset,
                features = c('NC-022518.1-HERVK113'),
                cols = c('lightgrey','red'),
                #keep.scale = 'all',
                pt.size = 0.5)+ggtitle(' Human endogenous retrovirus K113 complete genome:NC-022518.1')
p2<-DimPlot(data_subset,
            group.by = 'orig.ident',
            label = TRUE,
            label.size = 6,
            cols = c())+ggtitle('Timepoints')
p1+p2
cat('Epstein Barr Type I nuclear antigen leader protein is expressed in HERVK113+ cells')
VlnPlot(data_subset,
        features = c('EBNA-LP-ebv1'),
        layer = 'data',
        pt.size=3,
        group.by = 'viral.groups')+ ggtitle('EBNA-LP-ebv Type I is present HERVK113 infected cells')
cat('Epstein-Barr virus (EBV) nuclear antigen leader protein (EBNA-LP) plays a critical role in transformation of primary B lymphocytes to continuously proliferating lymphoblastoid cell lines (LCLs).
Reference: Kanamori M, Watanabe S, Honma R, Kuroda M, Imai S, Takada K, Yamamoto N, Nishiyama Y, Kawaguchi Y. Epstein-Barr virus nuclear antigen leader protein induces expression of thymus- and activation-regulated chemokine in B cells. J Virol. 2004 Apr;78(8):3984-93. doi: 10.1128/jvi.78.8.3984-3993.2004. PMID: 15047814; PMCID: PMC374277.')
#-------------------------------------------------------------------------------
data<-readRDS(file='<path_to>/B_cell_EBV_tmpts_exp/GSE282400_BcellEBV_viralpanel.rda')
data <- SetIdent(data, value = "orig.ident")
cat('Simian virus 40 complete genome')
cat('SV40 is a potent DNA tumor virus that induces tumors in rodents and transforms many types of cells in culture, including those of human origin.
Reference: Janet S. Butel, John A. Lednicky, Cell and Molecular Biology of Simian Virus 40: Implications for Human Infections and Disease , JNCI: Journal of the National Cancer Institute, Volume 91, Issue 2, 20 January 1999, Pages 119–134, https://doi.org/10.1093/jnci/91.2.119')
cat('SV40+ cells are less abundant in plasma cells')
VlnPlot(data,
        features = c(sv40),
        layer = 'data',
        pt.size=3,
        group.by = 'fine_celltype')+ ggtitle('Simian virus 40 complete genome:NC-001669.1')
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-001669.1-SV40`>0)
cat('Simian virus 40+ cells increase from day0 to day8')
table(Idents(infected_cells))
VlnPlot(data,
        features = c(sv40),
        layer = 'data',
        pt.size=3,
        group.by = 'orig.ident')+ ggtitle('Simian virus 40 complete genome:NC-001669.1')
cat('Our observation that SV40 may successfully infect and persist within human B cells and constitutes an important prerequisite for further assessment of whether the virus has a role in the development of B-cell lymphoproliferative disorders. No evidence of SV40 sequences was found in the marmoset B95.8 cell line from which infectious EBV virions were produced and in mock samples, thus ruling out laboratory contamination.')
cat('Dolcetti R, Martini F, Quaia M, Gloghini A, Vignocchi B, Cariati R, Martinelli M, Carbone A, Boiocchi M, Tognon M. Simian virus 40 sequences in human lymphoblastoid B-cell lines. J Virol. 2003 Jan;77(2):1595-7. doi: 10.1128/jvi.77.2.1595-1597.2003. PMID: 12502874; PMCID: PMC140833')
all.genes<-row.names(data)
cat('Differentially expressed genes comparing day8 to day2 in infected cells only')
DEG<-FindMarkers(infected_cells,
                 ident.1 = 'day8',
                 ident.2 = 'day2',
                 features=c(all.genes),
                 logfc.threshold = 0.5,
                 min.pct=0.4,
                 only.pos=FALSE)
head(DEG)
cat('genes for DAVID Bioinformatics functional annotations')
cat(row.names(DEG))
DoHeatmap(infected_cells,
          features = c(row.names(DEG)),
          size = 4)+ ggtitle('Simian virus 40 infected cells only
                             ')
cat('get entrezgene ids for functional analysis')
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
geneid <- row.names(DEG)
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
cat('Reactome pathway DEG annotations')
reactome <- enrichPathway(genes$entrezgene_id,
                          organism = 'human')
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Simian virus 40 infected
cells only. DEGs day8 v day2',
        label_format=50)
edox <- setReadable(reactome, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox,
         foldChange=NULL,
         node_label='all',
         colorEdge=TRUE,
         cex_label_gene=0.1,
         circular=TRUE,
         cex_gene=0.1,
         showCategory = 8,
         node_label_size=1,
         cex_label_category=1)+ ggtitle(' Reactome pathways:Simian virus 40 infected
 cells. DEGs day8 v day2')
cat('Biological Process DEG annotations')
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        label_format=50)+ ggtitle(' GO:Biological Process:Simian virus 40 infected cells
 DEGs day8 v day0')
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
cat('Differentially expressed genes comparing infected and uninfected cells')
DEG2<-FindMarkers(data,
                  ident.1 = 'viral.pos',
                  ident.2 = 'viral.neg',
                  logfc.threshold = 0.4,
                  min.pct = 0.3,
                  test.use = 'bimod',
                  only.pos=FALSE)
head(DEG2)
cat('Even out subsets')
data_subset<-subset(data,
                    downsample=352)
table(Idents(data_subset))
DoHeatmap(data_subset,
          features = c(row.names(DEG2)),
          size = 4)+ ggtitle('Simian virus 40 infected v non-infected cells')
cat('genes for DAVID Bioinformatics functional annotations')
cat(row.names(DEG2))
geneid <- row.names(DEG2)
head(geneid)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
cat('Reactome pathway DEG annotations')
reactome <- enrichPathway(genes$entrezgene_id,
                          organism = 'human')
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Simian virus 40 
Infected v non-infected genes',
        label_format=50)
edox <- setReadable(reactome, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox,
         foldChange=NULL,
         node_label='all',
         colorEdge=TRUE,
         cex_label_gene=0.5,
         circular=TRUE,
         cex_gene=0.1,
         showCategory = 3,
         node_label_size=1,
         cex_label_category=1)+ ggtitle(' Reactome pathways:Simian virus 40 
 Infected v non-infected genes')
cat('Biological Process DEG annotations')
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        showCategory=6,
        label_format=50)+ ggtitle(' GO:Biological Process:Simian virus 40
 Infected v non-infected genes')
p1<-FeaturePlot(data_subset,
                features = c('NC-001669.1-SV40'),
                cols = c('lightgrey','red'),
                #keep.scale = 'all',
                pt.size = 0.5)+ggtitle('Simian virus 40')
p2<-DimPlot(data_subset,
            group.by = 'orig.ident',
            label = TRUE,
            label.size = 5,
            cols = c())+ggtitle('Timepoints')
p1+p2
cat('Epstein Barr Type I nuclear antigen leader protein is expressed in SV40+ cells')
VlnPlot(data_subset,
        features = c('EBNA-LP-ebv1'),
        layer = 'data',
        pt.size=3,
        group.by = 'viral.groups')+ ggtitle('EBNA-LP-ebv Type I in SV40 infected cells')
#---------------------------------------------------------------------
data<-readRDS(file='<path_to>/B_cell_EBV_tmpts_exp/GSE282400_BcellEBV_viralpanel.rda')
data <- SetIdent(data, value = "orig.ident")
cat('Stx converting phage vB_EcoS_ST2-8624, complete genome')
cat('Viruses; Duplodnaviria; Heunggongvirae; Uroviricota;
Caudoviricetes; Sepvirinae; Traversvirus; Traversvirus ST28624.')
cat('Submitted by: Laboratory of Molecular Biology, Institute of Biochemistry and Biophysics, Polish Academy of Sciences, Wita Stwosza 59, Gdansk 80-308, Poland')
cat('Shiga toxins (Stx) comprise a family of potent cytotoxins that are involved in severe human disease. Stx are mainly produced by Escherichia coli isolated from human and nonhuman sources, and by Shigella dysenteriae type 1.
Reference: Schmidt H. Shiga-toxin-converting bacteriophages. Res Microbiol. 2001 Oct;152(8):687-95. doi: 10.1016/s0923-2508(01)01249-9. PMID: 11686382.')
table(Idents(data))
VlnPlot(data,
        features = c(stxPhage),
        layer = 'data',
        pt.size=3)+ggtitle('Stx converting phage vB_EcoS_ST2-8624, 
                  complete genome:NC-042057.1')
cat('Isolate only virally infected cells')
infected_cells<-subset(data, subset = `NC-049922.1-stxPhage`>0)
cat('Stx converting phage vB+ cells increase from day0 to day8')
table(Idents(infected_cells))
cat('Stx converting phage vB reads are present in all B-cell subsets')
VlnPlot(infected_cells,
        features = c(stxPhage),
        layer='data',
        group.by = 'fine_celltype',
        pt.size = 3)+ggtitle('Stx converting phage vB_EcoS_ST2-8624,
                  complete genome:NC-042057.1')
all.genes<-row.names(data)
cat('Differentially expressed genes comparing day8 to day2')
DEG<-FindMarkers(infected_cells,
                 ident.1 = 'day8',
                 ident.2 = 'day2',
                 features=c(all.genes),
                 logfc.threshold = 0.4,
                 min.pct=0.4,
                 #group.by = 'gen_celltype',
                 only.pos=FALSE)
head(DEG)
cat('genes for DAVID bioinformatics functional annotations')
cat(row.names(DEG))
DoHeatmap(infected_cells,
          features = c(row.names(DEG)),
          size = 4)+ ggtitle('Stx converting phage vB_EcoS_ST2-8624 infected cells only
                             ')
cat('Get entrez gene ids for functional annotations')
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
cat('Reactome pathway DEG annotations')
reactome <- enrichPathway(genes$entrezgene_id,
                          organism = 'human')
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Stx converting phage 
vB_EcoS_ST2-8624 infected cells only. DEGs comparing 
day8 v day2',
        label_format=50)
edox <- setReadable(reactome, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox,
         foldChange=NULL,
         node_label='all',
         colorEdge=TRUE,
         cex_label_gene=0.5,
         circular=TRUE,
         cex_gene=0.1,
         showCategory = 3,
         node_label_size=1,
         cex_label_category=1)+ ggtitle(' Reactome pathways:Stx converting phage vB_EcoS_ST2-8624 infected cells only
 DEGs day8 v day2')
cat('Biologicial process annotations')
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
dotplot(go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        label_format=50)+ ggtitle(' GO:Biological Process:Stx converting phage 
 vB_EcoS_ST2-8624 infected cells only. DEGs 
 day8 v day2')
cat('Split into uninfected and infected cells')
data$viral.groups <- 'viral.pos'
data$viral.groups[WhichCells(data, expression= `NC-049922.1-stxPhage` < 0.1)] <- 'viral.neg'
DimPlot(data, reduction = 'pca',split.by = 'viral.groups')
VlnPlot(data,
        features = c(stxPhage),
        group.by = 'viral.groups')+ggtitle('Stx converting phage vB_EcoS_ST2-8624:NC-049922.1')
head(data@meta.data)
data <- SetIdent(data, value = "viral.groups")
table(Idents(data))
cat('Differentially expressed genes comparing viral positive to negative cells')
DEG2<-FindMarkers(data,
                  ident.1 = 'viral.pos',
                  ident.2 = 'viral.neg',
                  logfc.threshold = 0.5,
                  min.pct = 0.4,
                  test.use = 'bimod',
                  only.pos=FALSE)
head(DEG2)
cat('Even out subsets')
data_subset<-subset(data,
                    downsample=115)
DoHeatmap(data_subset,
          features = c(row.names(DEG2)),
          size = 4)+ ggtitle('Stx converting phage vB_EcoS_ST2-8624 infected v non-infected cells
                             ')
cat('genes for DAVID Bioinformatics functional analysis')
cat(row.names(DEG2))
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
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways:Stx converting phage 
vB_EcoS_ST2-8624 Infected v non-infected genes',
        label_format=50)
table(Idents(data_subset))
edox <- setReadable(reactome, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox,
         foldChange=NULL,
         node_label='all',
         colorEdge=TRUE,
         cex_label_gene=0.5,
         circular=TRUE,
         cex_gene=0.1,
         showCategory = 3,
         node_label_size=1,
         cex_label_category=1)+ ggtitle(' Reactome pathways:Stx converting phage 
 vB_EcoS_ST2-8624 Infected v non-infected genes')
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "BP")
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
p1<-FeaturePlot(data_subset,
                features = c(stxPhage),
                cols = c('lightgrey','red'),
                #keep.scale = 'all',
                pt.size = 0.5)+ggtitle('Stx converting phage vB_EcoS_ST2-8624, complete genome')
p2<-DimPlot(data_subset,
            group.by = 'orig.ident',
            label = TRUE,
            label.size = 5,
            cols = c())+ggtitle('Timepoints')
p1+p2
cat('Epstein Barr Type I nuclear antigen leader protein is expressed in Stx converting phage vB+ cells')
VlnPlot(data_subset,
        features = c('EBNA-LP-ebv1'),
        layer = 'data',
        pt.size=3,
        group.by = 'viral.groups')+ ggtitle('EBNA-LP-ebv Type I in Stx converting phage 
vB_EcoS_ST2-8624 infected cells')
