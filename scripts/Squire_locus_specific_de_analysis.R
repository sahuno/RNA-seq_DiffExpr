#squire locus specific de analysis
library(data.table)
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(viridis)
library(magrittr)
library(pheatmap)
library(ggrepel)
library(ggfortify)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
library(cowplot)
library(RColorBrewer)
library(NbClust)
options("width"=200)

rna_path <- "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA"
counts_path <- paste0(rna_path, "/CT/squire_te_fwd.tsv")
# counts_df <- read_tsv(counts_path)
counts_df <- fread(counts_path)

dim(counts_df)
rownames(counts_df) <- counts_df$te.id
head(counts_df)
counts_df_2 <- counts_df %>% dplyr::select(!c("te.id", "te.name"))


workflow_dir <- "/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/"
metadata_path <- paste0("/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/metadata_triplicates_recoded.csv")
blind_transform <- TRUE #should the rlog transformation be blind
# drop_samples <- "R.S.2" #samples to drop from analysis
drop_samples <- c("R.S.2", "R.C.3") #samples to drop from analysis
# drop_samples <- NULL #samples to drop from analysis
ref_variable = "DMSO"
MIN_ReadsCounts = 50
smallestGroupSize <- 3

source(paste0(workflow_dir,"scripts/helper_functions.R"))
metadata_df <- read.csv(file = metadata_path, sep="," ,header = TRUE)

metadata_df <- metadata_df[!metadata_df$samples %in% drop_samples,]
# counts_df <- counts_df[,..metadata_df$samples]
counts_df_2 <- counts_df_2 %>% dplyr::select(!all_of(drop_samples))


# dds <- DESeqDataSetFromMatrix(countData = data.matrix(counts_df_2),
#                               colData = metadata_df,
#                               design = ~ condition)

# sapply(data.matrix(counts_df_2), is.na)
# counts_df_2_rowSums <- rowSums(counts_df_2)
# counts_df_2_rowSums_less1 <- counts_df_2_rowSums[counts_df_2_rowSums < 1.0]
# counts_df_2_rowSums_less1[counts_df_2_rowSums_less1 > 0.00]


# dim(counts_df_2)
#make a matrix out of squire counts
counts_df_te.id <- counts_df$`te.id`
counts_df_matrix <- data.matrix(counts_df[,!c(1,2)])
rownames(counts_df_matrix) <- counts_df_te.id
counts_df_matrix_noZeros <- counts_df_matrix[rowSums(counts_df_matrix) > 1,]
# dim(counts_df_noZeros)
counts_df_matrix_noZeros_int <- apply(counts_df_matrix_noZeros, 2, as.integer)

counts_df_matrix_noZeros_int_samplesRemov <- counts_df_matrix_noZeros_int[,-which(colnames(counts_df_matrix_noZeros_int) %in% drop_samples)]


dds <- DESeqDataSetFromMatrix(countData = data.matrix(counts_df_matrix_noZeros_int_samplesRemov),
                              colData = metadata_df,
                              design = ~ condition)
message("setting CTRL as reference group")
dds$condition <- relevel(dds$condition, ref = ref_variable)
keep <- rowSums(counts(dds) >= MIN_ReadsCounts) >= smallestGroupSize
# keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
dds <- dds[keep,] #filter dds
# head(counts(dds))
names(assays(dds))
names(mcols(dds))
dds <- DESeq(dds)


##########################
### get cordinates
TEcounts <- fread("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/mapping/R-0-1/squire/aln_TEcounts.txt.gz")
# make a bed file of TE cordinates from the TEcounts 
namesTEs <- c("TE_chr",   "TE_start",   "TE_stop" ,  "TE_name" ,   "TE_strand",  "milliDiv" ,  "tx_chr",      "tx_start",    "tx_stop"  ,   "TE_ID"   ,    "fpkm" ,    "tx_strand"  , "Sample"   ,   "alignedsize","uniq_counts", "tot_counts" , "tot_reads",   "score")
TEcounts_ColnamesRearranged <- TEcounts[, ..namesTEs]

setnames(TEcounts_ColnamesRearranged, c("TE_chr",   "TE_start",   "TE_stop" ,  "TE_name" ,   "TE_strand"), c("chr",   "start",   "end" ,  "TE_name" ,   "strand"))

fwrite(TEcounts_ColnamesRearranged, "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/R01_TEcounts_DNA_RNA.bed", sep="\t", quote = FALSE, row.names = FALSE)