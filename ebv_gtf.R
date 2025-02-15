library(tidyr)
library(tidyverse)
library(dplyr)
#---read in fasta file
file<-read.delim('/home/em/Downloads/B_cell_EBV_tmpts_exp/ebv.fa',
                 header = FALSE,sep = ' ')
#---retrieve sequence ids
selected_rows <- file[grep(">", file$V1), ]
contig_names <- gsub('>','',selected_rows$V1)
contig_names
#---check if sequence ids have duplicates
duplicated(contig_names)
#---retrieve CDS sequence names
gene_names<-gsub(".*[=]([^.]+)[]].*", "\\1",selected_rows$V2)
gene_names
#---check if CDS sequencs names have duplicates
duplicated(gene_names)
#----number duplicated gene names
df <- data.frame("GENES"=c(gene_names)) 
gene_names<-df |> 
  mutate(dupl = if_else(duplicated(GENES), 1, 0)) |> 
  group_by(GENES) |> 
  mutate(dupl = cumsum(dupl),
         GENES = paste(GENES, dupl, sep = "-")) |> 
  select(-dupl)
gene_names<-gene_names$GENES
gene_names
duplicated(gene_names)
gene_names<-gsub("-0","",as.character(gene_names))
#---attach identifier to coding sequence names
gene_names<-paste0(gene_names, "-ebv1")
#gene_names <- sub("_", "-", gene_names)
gene_names
#---write gene names to csv file
write.csv(gene_names,file='/home/em/Downloads/B_cell_EBV_tmpts_exp/ebv_genes.csv')
#---read in CDS sequence sizes
file2<-read.delim('/home/em/Downloads/B_cell_EBV_tmpts_exp/sizes.genome',
                  header = FALSE)
#---make gtf file
gene_sizes<-file2$V2
gene_sizes
gene_names_edit<-paste0('"',gene_names,'";')
gene_names_edit
genes_transcripts<-cbind('gene_id',gene_names_edit,'transcript_id',gene_names_edit,
                         'gene_name',gene_names_edit,'gene_biotype "protein_coding";')
genes_transcripts<-data.frame(genes_transcripts)
genes_transcripts_final <- paste(genes_transcripts$V1,
                                 genes_transcripts$gene_names_edit,
                                 genes_transcripts$V3,
                                 genes_transcripts$gene_names_edit.1,
                                 genes_transcripts$V5,
                                 genes_transcripts$gene_names_edit.2,
                                 genes_transcripts$V7)
final<-cbind(contig_names,
             'unknown',
             'exon',
             1,
             gene_sizes,
             '.',
             '+',
             '.',
             genes_transcripts_final)
#---view finalized gtf
final<-data.frame(final)
head(final)
write.table(data.frame(final),'/home/em/Downloads/B_cell_EBV_tmpts_exp/ebv.gtf',
            sep="\t",
            row.names=FALSE,
            quote = FALSE,
            col.names = FALSE)

