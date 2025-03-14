library(tidyr)
library(tidyverse)
library(dplyr)
#---read in fasta file
file<-read.delim('<path_to>/B_cell_EBV_tmpts_exp/NC_001669.1_SV40.fa',
                 header = FALSE,sep = ' ')
#---retrieve sequence ids
selected_rows <- file[grep(">", file$V1), ]
contig_names <- gsub('>','',selected_rows$V1)
contig_names
#---use RefSeq id as name
name<-contig_names
#---attach identifier to sequence name
name<-paste0(name, "-SV40")
name <- sub("_", "-", name)
name
#---write gene names to csv file
write.csv(name,file='<path_to>/B_cell_EBV_tmpts_exp/SV40_sequence.csv')
#---read in CDS sequence sizes
file2<-read.delim('<path_to>/B_cell_EBV_tmpts_exp/sizes.genome',
                  header = FALSE)
#---make gtf file
gene_sizes<-file2$V2
gene_sizes
gene_names_edit<-paste0('"',name,'";')
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
write.table(data.frame(final),'<path_to>/B_cell_EBV_tmpts_exp/NC_001669.1_SV40.gtf',
            sep="\t",
            row.names=FALSE,
            quote = FALSE,
            col.names = FALSE)
