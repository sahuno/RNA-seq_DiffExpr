#name: samuel ahuno
#load packages
library(DESeq2) #main desdeq analysis
library(tidyverse)
library(EnhancedVolcano) # used for making volcano plots
library(viridis)
library(magrittr)
library(pheatmap)
library(data.table)
library(ggrepel)
library(ggfortify)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
library(cowplot)
library(RColorBrewer)
library(NbClust)
options("width"=200)

########################################################################################
# set prroject  parameters
########################################################################################
source_dir = "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/" 
setwd("/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/sandbox/")
blind_transform <- TRUE #should the rlog transformation be blind
drop_samples <- NULL #samples to drop from analysis
ref_variable = "DMSO"
smallestGroupSize <- 3

message("workig dir is - ", getwd())
dir.create(paste0("figures/"), recursive = TRUE, showWarnings = TRUE)
dir.create(paste0("data/"), recursive = TRUE, showWarnings = TRUE)
dir.create(paste0("figures/geneWiseNormalizedCounts"), recursive = TRUE, showWarnings = TRUE)
########################################################################################
## #read in data, extract coding genes counts 
########################################################################################
#load datasets from rna-seq results folder - sasha
counts_annot <- fread(paste0(source_dir,'CT/counts_annot.tsv'))
# counts_coding_genes_df <- fread(paste0(source_dir,'CT/counts.tsv'))
countsOnlyDF <- counts_annot[,11 : ncol(counts_annot)]

### Step 2; get sample metadata
#note; we assume you already have sample metadata file saved some where
metadata_df <- read.csv(file = paste0(source_dir,'metadata_triplicates.csv'),sep="," ,header = TRUE)
metadata_df$condition <- factor(metadata_df$condition)

#create dseq object
ddsGenesAndTEs <- DESeqDataSetFromMatrix(countData = countsOnlyDF,
                              colData = metadata_df,
                              design = ~ condition)

#set the reference group
message("setting CTRL as reference group")
ddsGenesAndTEs$condition <- relevel(ddsGenesAndTEs$condition, ref = ref_variable)
mcols(ddsGenesAndTEs) <- DataFrame(mcols(dObjectGenesAndTEs), counts_annot[,1:11]) #add gene annotations to the dds object

dObjectGenesAndTEs <- DESeq(ddsGenesAndTEs) #this command estimates the size factors and dispersion estimates
# idxOfRepeats <- which(mcols(dObjectGenesAndTEs)$gene.id %in% counts_coding_genes_df$gene.id) 

# mcols(dObjectGenesAndTEs)[which(mcols(dObjectGenesAndTEs)$gene.id %in% counts_coding_genes_df$gene.id),]
# dObjectTEs <- dObjectGenesAndTEs[idxOfRepeats,]

# dObjectTEs <- estimateDispersionsGeneEst(dObjectTEs)
# dispersions(dObjectTEs) <- mcols(dObjectTEs)$dispGeneEst
# rlogTEs <- rlog(dObjectTEs, blind = blind_transform)

####isolate repeat elements. 
gs.rep <- annot %>% dplyr::filter(chr == "REP") %>% 
    dplyr::select(c(gene.id, gene.type)) %>% 
      mutate(tp = case_when(str_detect(gene.type,"^(Endogenous|ERV|LTR|Lomg_terminal_repeat|Endogenous_Retrovirus)") ~ "LTR", 
                            str_detect(gene.type,"^(Retroposon/SVA|MSAT)") ~ "SVA", 
                            str_detect(gene.type,"^(LINE|_LINE1|L1|L2|_LINE|RTE|CR1|R4|Penelope|Non-LTR_Retrotransposon)") ~ "LINE", 
                            str_detect(gene.type,"^tRNA") ~ "tRNA", 
                            str_detect(gene.type,"^rRNA") ~ "rRNA", 
                            str_detect(gene.type,"^SINE") ~ "SINE", 
                            str_detect(gene.type,"^(SAT|satellite|Satellite)") ~ "SAT", 
                            str_detect(gene.type,"^(DNA|hAT|_DNA|piggyBac|Mariner/Tc1|MuDR)") ~ "DNA", 
                            str_detect(gene.type,"^(X17_DNA)") ~ "DNA", 
                            str_detect(gene.type,"^(Simple_repeat|Low_complexity)") ~ "Simple_repeat", 
                            str_detect(gene.type,"^(scRNA|snRNA|srpRNA|RNA)") ~ "sRNA", 
                            str_detect(gene.type,"^Pseudogene") ~ "Pseudogene", 
                            TRUE ~ "Unclassified"))

coding_genes_df <- counts_annot %>% dplyr::filter(!(gene.id %in% gs.rep$gene.id)) #%>% dim()

idxOfRepeats <- which(mcols(dObjectGenesAndTEs)$gene.id %in% gs.rep$gene.id) 
mcols(dObjectGenesAndTEs)[which(mcols(dObjectGenesAndTEs)$gene.id %in% gs.rep$gene.id),]
dObjectGenesAndTEs[idxOfRepeats,]

dObjectTEs <- dObjectGenesAndTEs[idxOfRepeats,]

rlogTEs <- rlog(dObjectTEs, blind = blind_transform)


pltPCA_samples <- plotPCA(rlogTEs, intgroup=c("condition", "samples"))
pltPCA_ConditionsOnly <- plotPCA(rlogTEs, intgroup=c("condition"))
ggsave(pltPCA_samples, file = paste0("figures/PCA_by_condition_and_samples_Repeats.png"), width = 9, height = 7)
ggsave(pltPCA_ConditionsOnly, file = paste0("figures/PCA_by_conditionOnly_Repeats.png"), width = 9, height = 7)

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rlogTEs)
rld_cor <- cor(rld_mat)
rld_df <- rld_mat %>% as.data.frame() %>% rownames_to_column("gene.id")

#save a copy of the rlog matrix
# fwrite(rld_df, file=paste0("data/DESeq2_rlog_Transform_Blind",blind_transform,".tsv"), sep="\t")


# Plot heatmaps
pltHeatmap_samples <- pheatmap(rld_cor, annotation = metadata_rown_df[,c("condition"), drop = FALSE], scale = "row")
ggsave(pltHeatmap_samples, file = paste0("figures/heatmap_conditions_Repeats.png"), width = 12, height = 12)


# Compute the distance matrix of samples
sampleDists <- as.dist(1 - cor(rld_mat)) #what does, closet to 1 mean? similar samples
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, annotation = metadata_rown_df[,c("condition"), drop = FALSE] ,filename = "figures/sampleDistMatrix.pdf")
