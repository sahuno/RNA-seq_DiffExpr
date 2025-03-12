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
library(GenomicRanges)

options("width"=200)


###general approach is to run 
#1: `squire_te_fwd.tsv` this is the counts file you need
#2: normalized with coding genes
#remove coding genes
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

#check if there's is line 1
#filter out L1
countsL1_df <- counts_df %>% dplyr::filter(str_detect(te.id, "L1"))
countsL1_df_withoutHLA <- counts_df %>% dplyr::filter(str_detect(te.id, "L1")) #!str_detect(te.id, "HLA")

#check length of line1
HeadcountsL1_df <- head(countsL1_df, 10)[,c(1:5)]
# HeadcountsL1_df <- countsL1_df[,c(1:5)] #uncomment to run all

# HeadcountsL1_df
HeadcountsL1_df_sep <- HeadcountsL1_df %>% separate_wider_delim(te.name, "|" ,names = c("chr", "start", "end", "name","number","strand"))
HeadcountsL1_df_sep <- HeadcountsL1_df_sep %>% separate_wider_delim(name, ":" ,names = c("consensus","clade", "class"))
L1HS_filter <- HeadcountsL1_df_sep %>% dplyr::filter(str_detect(consensus, "L1Md")) #!str_detect(te.id, "HLA")
unique(L1HS_filter$consensus)

unique(HeadcountsL1_df_sep$consensus)
HeadcountsL1_gr_sep <- makeGRangesFromDataFrame(HeadcountsL1_df_sep, keep.extra.columns = TRUE)
width(HeadcountsL1_gr_sep)
# counts_dfp
# dim(counts_df_2)
#make a matrix out of squire counts
counts_df_te.id <- counts_df$`te.id`
counts_df_matrix <- data.matrix(counts_df[,!c(1,2)])
rownames(counts_df_matrix) <- counts_df_te.id
counts_df_matrix_noZeros <- counts_df_matrix[rowSums(counts_df_matrix) > 1,]
# dim(counts_df_noZeros)
counts_df_matrix_noZeros_int <- apply(counts_df_matrix_noZeros, 2, as.integer)
rownames(counts_df_matrix_noZeros_int) <- rownames(counts_df_matrix_noZeros) #put geneIDs back

#drop low variable samples
counts_df_matrix_noZeros_int_samplesRemov <- counts_df_matrix_noZeros_int[,-which(colnames(counts_df_matrix_noZeros_int) %in% drop_samples)]


keep_TEs <- rowSums(counts_df_matrix_noZeros_int_samplesRemov >= MIN_ReadsCounts) >= smallestGroupSize
# keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
counts_df_matrix_noZeros_int_samplesRemov <- counts_df_matrix_noZeros_int_samplesRemov[keep_TEs,] 




################################################################################################
#load and subset protein coding datasets from rna-seq results folder - sasha
################################################################################################
# counts_annot <- fread(paste0(rna_path,'/CT/counts_annot.tsv'))
counts_annot <- read_tsv(paste0(rna_path,'/CT/counts_annot.tsv'))
counts_annot_proteinCodingOnly <- counts_annot %>% dplyr::filter(gene.type == "protein_coding")
countsOnly_proteinCodingOnly_DF <- counts_annot_proteinCodingOnly[,11 : ncol(counts_annot_proteinCodingOnly)]
# unique(counts_annot$gene.type)
countsOnly_proteinCodingOnly_DF_removeOutlierSamples <- countsOnly_proteinCodingOnly_DF %>% dplyr::select(!all_of(drop_samples))
# dim(countsOnly_proteinCodingOnly_DF_removeOutlierSamples)

keep <- rowSums(countsOnly_proteinCodingOnly_DF_removeOutlierSamples >= MIN_ReadsCounts) >= smallestGroupSize
# keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
countsOnly_proteinCodingOnly_DF_removeOutlierSamples <- countsOnly_proteinCodingOnly_DF_removeOutlierSamples[keep,] #filter dds
countsOnly_proteinCodingOnly_DF_removeOutlierSamples <- as.data.frame(countsOnly_proteinCodingOnly_DF_removeOutlierSamples) #tibble doesn't allow rownames so first convert to data.frame
rownames(countsOnly_proteinCodingOnly_DF_removeOutlierSamples) <- counts_annot_proteinCodingOnly$gene.id[keep]

#cbind repeats and protein codinggenes
countsRepeats_codingGenes <- rbind(countsOnly_proteinCodingOnly_DF_removeOutlierSamples, counts_df_matrix_noZeros_int_samplesRemov)

##keep dmso & Aza only
metadata_df_AZA_DMSO <- metadata_df %>% dplyr::filter(condition %in% c("DMSO", "AZA"))

#make dseq2 objects
rownames(metadata_df) <- metadata_df$samples # this is important for the DESeqDataSetFromMatrix function
# metadata_df$condition <- factor(metadata_df$condition, ref)

countsRepeats_codingGenes_AZA_DMSO <- countsRepeats_codingGenes %>% dplyr::select(all_of(metadata_df_AZA_DMSO$samples))

metadata_df_AZA_DMSO$condition <- factor(metadata_df_AZA_DMSO$condition, levels = c("DMSO", "AZA"))



dds <- DESeqDataSetFromMatrix(countData = data.matrix(countsRepeats_codingGenes_AZA_DMSO),
                              colData = metadata_df_AZA_DMSO,
                              design = ~ condition)

message("setting CTRL as reference group")
# dds$condition <- relevel(dds$condition, ref = ref_variable)

#sanity checks
all(rownames(metadata_df) == colnames(countsRepeats_codingGenes))


# keep <- rowSums(counts(dds) >= MIN_ReadsCounts) >= smallestGroupSize
# # keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
# dds <- dds[keep,] #filter dds
# # head(counts(dds))
# names(assays(dds))
# names(mcols(dds))

#either i do this or run step by step
dds <- DESeq(dds)


#run DESeq step by step, normalizing with coding genes; 
dds_sf <- estimateSizeFactors(dds)
dds_sf_Disp <- estimateDispersions(dds_sf)
dds_sf_Disp_TEs <- dds_sf_Disp[!grepl("^ENSMUSG", rownames(dds_sf_Disp)), ] #now remove protein coding genes for the test
dds_sf_Disp_TEs <- nbinomWaldTest(dds_sf_Disp_TEs)



#do PCA to check sample for outliers
rlogTEs <- rlog(dds_sf_Disp_TEs, blind = blind_transform)
rlogTEs_withcoding <- rlog(dds, blind = blind_transform)

pltPCA_samples <- plotPCA(rlogTEs, intgroup=c("condition", "samples"))
pltPCA_ConditionsOnly <- plotPCA(rlogTEs, intgroup=c("condition"))
ggsave(pltPCA_samples, file = paste0("figures/PCA_by_condition_and_samples_locusSpecific_SquireRepeats_normalizedWithCoding.png"), width = 9, height = 7)
ggsave(pltPCA_ConditionsOnly, file = paste0("figures/PCA_by_conditionOnly_locusSpecific_SquireRepeats_normalizedWithCoding.png"), width = 9, height = 7)


pltPCA_samples2 <- plotPCA(rlogTEs_withcoding, intgroup=c("condition", "samples"))
pltPCA_ConditionsOnly2 <- plotPCA(rlogTEs_withcoding, intgroup=c("condition"))
ggsave(pltPCA_samples2, file = paste0("figures/PCA_by_condition_and_samples_locusSpecific_SquireRepeats_normalizedWithCoding_DESeq.png"), width = 9, height = 7)
ggsave(pltPCA_ConditionsOnly2, file = paste0("figures/PCA_by_conditionOnly_locusSpecific_SquireRepeats_normalizedWithCoding_DESeq.png"), width = 9, height = 7)


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rlogTEs)
rld_cor <- cor(rld_mat)
rld_df <- rld_mat %>% as.data.frame() %>% rownames_to_column("gene.id")







##########################
# ### get cordinates
# TEcounts <- fread("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/mapping/R-0-1/squire/aln_TEcounts.txt.gz")
# # make a bed file of TE cordinates from the TEcounts 
# namesTEs <- c("TE_chr",   "TE_start",   "TE_stop" ,  "TE_name" ,   "TE_strand",  "milliDiv" ,  "tx_chr",      "tx_start",    "tx_stop"  ,   "TE_ID"   ,    "fpkm" ,    "tx_strand"  , "Sample"   ,   "alignedsize","uniq_counts", "tot_counts" , "tot_reads",   "score")
# TEcounts_ColnamesRearranged <- TEcounts[, ..namesTEs]

# setnames(TEcounts_ColnamesRearranged, c("TE_chr",   "TE_start",   "TE_stop" ,  "TE_name" ,   "TE_strand"), c("chr",   "start",   "end" ,  "TE_name" ,   "strand"))

# fwrite(TEcounts_ColnamesRearranged, "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/R01_TEcounts_DNA_RNA.bed", sep="\t", quote = FALSE, row.names = FALSE)