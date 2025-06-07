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
library(ggrepel)
options("width"=200)



# Rscript /data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/scripts/coding_genes.r \
#   --source_dir "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/" \
#   --workflow_dir "/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/" \
#   --blind_transform TRUE \
#   --drop_samples "R.S.2,R.C.3" \
#   --ref_variable "DMSO" \
#   --min_read_counts 50 \
#   --smallest_group_size 3 \
#   --qc_metrics "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/triplicates_mouse/qc.tsv"


# /data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/


#   --sandbox_dir "/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/sandbox/" \


# Parse command-line options
library(optparse)
option_list <- list(
  make_option(c("-s", "--source_dir"), type="character",
              default="/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/",
              help="parent directory to RNA-seq count data"),
  make_option(c("-w", "--workflow_dir"), type="character",
              default="/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/",
              help="Path to workflow directory"),
  # make_option(c("--sandbox_dir"), type="character",
  #             default="/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/sandbox/",
              # help="Working sandbox directory"),
  make_option(c("--blind_transform"), type="logical", default=TRUE,
              help="Perform blind rlog transformation"),
  make_option(c("--drop_samples"), type="character", default="R.S.2,R.C.3",
              help="Comma-separated samples to drop"),
  make_option(c("--ref_variable"), type="character", default="DMSO",
              help="Reference condition variable"),
  make_option(c("--min_read_counts"), type="integer", default=50,
              help="Minimum read count threshold"),
  make_option(c("--smallest_group_size"), type="integer", default=3,
              help="Smallest group size for filtering"),
  make_option(c("--qc_metrics"), type="character",
              default="/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/triplicates_mouse/qc.tsv",
              help="Path to QC metrics TSV file"),
  make_option(c("--metadata_File"), type="character",
              default="/data1/greenbab/projects/triplicates_epigenetics_diyva/configs/metadata_triplicates_recoded_QSTAT_group.csv",
              help="Path metadata file")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Assign variables from options
source_dir <- opt$source_dir
workflow_dir <- opt$workflow_dir
message("setting working dir is - ", getwd())
setwd(getwd())

# Create output directories
dir.create("figures", recursive=TRUE, showWarnings=TRUE)
dir.create("data", recursive=TRUE, showWarnings=TRUE)
dir.create("figures/geneWiseNormalizedCounts", recursive=TRUE, showWarnings=TRUE)

blind_transform <- opt$blind_transform
drop_samples <- if (opt$drop_samples != "") strsplit(opt$drop_samples, ",")[[1]] else NULL
ref_variable <- opt$ref_variable
MIN_ReadsCounts <- opt$min_read_counts
smallestGroupSize <- opt$smallest_group_size

# Load helper functions
# source(file.path(workflow_dir, "scripts", "helper_functions.R"))

# Read data files
annot <- fread(file.path(source_dir, "annot.tsv"))
counts_annot <- fread(file.path(source_dir, "CT", "counts_annot.tsv"))
cts <- fread(file.path(source_dir, "CT", "counts.tsv"))
qc_metrics <- fread(opt$qc_metrics)

### source helper functions
source(paste0(opt$workflow_dir,"/scripts/de_helper_functions.R")) #split this script into a function

#/juno/work/greenbaum/projects/TRI_EPIGENETIC/RNASeq_DE_TriEpi
########################################################################################
# set prroject  parameters
########################################################################################

#get protein coding genes only
gs.coding_ids <- annot[gene.type=="protein_coding",.(gene.id)]
cts_coding_df <- cts[gs.coding_ids, on = "gene.id"] %>% as.data.frame()
cts_coding <- cts_coding_df %>% column_to_rownames(var = "gene.id") # #convert colnames to rowids
# head(cts_coding)

####isolate repeat elements. we don't want them in differential expression for now
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
########################################################################################




########################################################################################
### Step 2; get sample metadata
########################################################################################
#note; we assume you already have sample metadata file saved some where
metadata_df <- read.csv(file = paste0(opt$metadata_File),sep="," ,header = TRUE)
#metadata_df; needs 3 columns; samples-sample names matching rna-seq sample list, condition - treatment/ctrl, new_samples_name - new sample names

# TODO: check if metadata_df has the required columns
message(unique(metadata_df$condition_long))

#rename  samples = 
#Assign new sample names
cond_assign_new_sample_names = ncol(metadata_df) > 2 & any(str_detect(names(metadata_df), "new_samples_name"))

if(cond_assign_new_sample_names){
  #drop samples
if(!is.null(drop_samples)){
metadata_df <- metadata_df %>% dplyr::filter(!samples %in% drop_samples)
cts_coding <- cts_coding %>% dplyr::select(!all_of(drop_samples)) #%>% head()

## rename samples in matrix
message("Renaming samples in cts_coding matrix based on metadata_df$new_samples_name")
# 1. Check that every column in cts_coding has a matching "samples" entry:
old_cols <- colnames(cts_coding)
if (any(! old_cols %in% metadata_df$samples)) {
  warning("Some columns in cts_coding have no entry in metadata_df$samples:",
          paste(setdiff(old_cols, metadata_df$samples), collapse = ", "))
}
# 2. Create a vector of new names in the same order as cts_coding’s columns:
new_cols <- metadata_df$new_samples_name[ match(old_cols, metadata_df$samples) ]
# 3. Assign them back to cts_coding:
colnames(cts_coding) <- new_cols
cts_coding <- cts_coding[, !is.na(colnames(cts_coding))]
# cts_coding <- cts_coding %>% dplyr::select(all_of(metadata_df[,"new_samples_name"])) #%>% head()
}

# #assign new sample names
# #idx_samples <- match(metadata_df$samples, colnames(cts_coding)) #check if samples match
# # names(cts_coding) <- metadata_df[,"new_samples_name"][idx_samples]
# # 1. Check that every column in cts_coding has a matching "samples" entry:
# old_cols <- colnames(cts_coding)
# if (any(! old_cols %in% metadata_df$samples)) {
#   warning("Some columns in cts_coding have no entry in metadata_df$samples:",
#           paste(setdiff(old_cols, metadata_df$samples), collapse = ", "))
# }
# # 2. Create a vector of new names in the same order as cts_coding’s columns:
# new_cols <- metadata_df$new_samples_name[ match(old_cols, metadata_df$samples) ]
# # 3. Assign them back to cts_coding:
# colnames(cts_coding) <- new_cols
# cts_coding <- cts_coding[, !is.na(colnames(cts_coding))]
# metadata_rown_df <- metadata_df %>% column_to_rownames("new_samples_name")
metadata_rown_df <- metadata_df %>% column_to_rownames("new_samples_name")
} else{
#drop samples
# if(!is.null(drop_samples)){
# metadata_df <- metadata_df %>% dplyr::filter(!samples %in% drop_samples)
# }
metadata_rown_df <- metadata_df %>% column_to_rownames("samples")
}

message("\nmetadata_rown_df has ", nrow(metadata_rown_df), " rows and ", ncol(metadata_rown_df), " columns\n")
print(head(metadata_rown_df))

message("Checking if all samples in metadata_rown_df are present in cts_coding")
print(head(cts_coding))
##Just commented
#cts_coding <- cts_coding %>% dplyr::select(any_of(rownames(metadata_rown_df))) #%>% head() #drop samples from counts matrix
##################################################################
############ Plot QC
##################################################################
message("Plotting QC metrics")
df_qc <- merge(qc_metrics, metadata_rown_df, by.x = "sample", by.y = "samples")
pltQC <- ggplot(df_qc, aes(x=sample, y=input.reads, fill=condition)) + 
geom_bar(stat="identity") + theme_minimal() + 
theme(axis.text.x = element_text(angle = 45)) + 
labs(title = "Total reads per sample", y = "Total reads") + theme(legend.position = "top")
ggsave(pltQC, file = "figures/total_reads_per_sample.pdf", width = 12, height = 12)

qc_plots <- c("input.reads", "uniquely.mapped", "total.counts", "rRNA" ,"mito.coding" ,"coding", "norm.factor",  "eff.lib.size",   "jct_fwd", "jct_rev", "unassigned_rev","QC" ,"median.degeneracy", "coding.dedup")

pdf("figures/qc_plots.pdf")
for(i in qc_plots){
p <- ggplot(df_qc, aes(x=sample, y=.data[[i]], fill=condition)) + 
geom_bar(stat="identity") + theme_minimal() + 
theme(axis.text.x = element_text(angle = 45)) + 
labs(title = i, y = i) + theme(legend.position = "top")
print(p)
}
dev.off()
# qc_metrics %>% dplyr::left_join(metadata_rown_df, by = c("sample" , "samples"))

##exclude samples from dseq analysis
# cts_coding <- cts_coding %>% dplyr::select(!veh_2) #remove outlier samples `veh_2`
##########################
#function to make dseq object
#import `scripts/helper_functions.R` for make_dseq_obj
#create dseq object and extract contrasts
dseq_func_out <- make_dseq_obj(ref_col = ref_variable, add_extra_contrasts = TRUE)
dds <- dseq_func_out[["dSeqObj"]]
contrasts_ls <- dseq_func_out[["constrasts"]]

#from mike love, how to normalize
dds <- DESeq(dds) #this command estimates the size factors and dispersion estimates

message("writing normalized counts to file")
#write files for comparism with methylation
source(paste0(workflow_dir,"scripts/write_normalizedCounts.R"))
# /data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/scripts/write_normalizedCounts.R
################################################################################################
################################################################################################
##pre - filtering
# rowSums(counts(dds)) <= 10
keep <- rowSums(counts(dds) >= MIN_ReadsCounts) >= smallestGroupSize
# keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
dds <- dds[keep,] #filter dds
# head(counts(dds))
names(assays(dds))
names(mcols(dds))
# keep <- rowSums(counts(dds_2_filterCounts) >= 10) >= smallestGroupSize #keep only counts where rowSums is above 10

################################################################################################
############################### Step 2; Run deseq2 analysis
##############################################################################################
#run DESeq analysis - normalization and filtering
dds <- DESeq(dds)
# names(mcols(dds))
# res <- results(dds, alpha=0.05) #use alpha of 0.05
# summary(res)
# results_per_contrast <- results(dds, alpha=0.05, contrast="condition_PARPi_vs_CTRL") #use alpha of 0.05
# contrast_input="condition_PARPi_vs_CTRL"
#results_per_contrast <- results(dds, alpha=0.05, contrast=c("condition", "PARPi" ,"CTRL"))
# resultsNames(results_per_contrast)
# res_test <- lfcShrink(dds, 
#      coef = 2, #apeglm only works with `coef` argument`
#      type = "apeglm"
#      )	

# Transform counts for data visualization
rld <- rlog(dds, blind=blind_transform)
# rld_mat <- assay(rld)

# Plot PCA 
pltPCA_samples <- plotPCA(rld, intgroup=c("condition", "samples")) + ggtitle("all samples & conditions")
pltPCA_ConditionsOnly <- plotPCA(rld, intgroup=c("condition")) + ggtitle("all samples & conditions")
ggsave(pltPCA_samples, file = paste0("figures/PCA_by_condition_and_samples_",paste0(drop_samples, collapse="_"),"Removed.pdf"), width = 12, height = 12)
ggsave(pltPCA_ConditionsOnly, file = paste0("figures/PCA_by_conditionOnly_",paste0(drop_samples, collapse="_"),".pdf"), width = 12, height = 12)

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
rld_df <- rld_mat %>% as.data.frame() %>% rownames_to_column("gene.id")

#save a copy of the rlog matrix
fwrite(rld_df, file=paste0("data/DESeq2_rlog_Transform_Blind",blind_transform,".tsv"), sep="\t")


# Plot heatmaps
pltHeatmap_samples <- pheatmap(rld_cor, annotation = metadata_rown_df[,c("condition"), drop = FALSE], scale = "row", main = "Heatmap of rlog transformed counts")
ggsave(pltHeatmap_samples, file = paste0("figures/heatmap_conditions_samples_",paste0(drop_samples, collapse="_"),"_dropped.pdf"), width = 12, height = 12)


# Compute the distance matrix of samples
message("Computing distance matrix of samples")
sampleDists <- as.dist(1 - cor(rld_mat)) #what does, closet to 1 mean? similar samples
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, annotation = metadata_rown_df[,c("condition"), drop = FALSE], main = "sample distance of rlog transformed counts",filename = paste0("figures/sampleDistMatrix_samples_",paste0(drop_samples, collapse="_"),"_dropped.pdf"))



res_rld_mat_df <- rld_mat %>%
  as.data.frame() %>%
  rownames_to_column(var="gene.id") %>% 
  as_tibble()

res_rld_mat_df_pivtLonger <- res_rld_mat_df %>% pivot_longer(!c("gene.id"), names_to = "samples", values_to = "rlog_transformed_counts")
#head(res_rld_mat_df_pivtLonger)

#slotNames(dds)
# dds_copy <- dds
# # colnames(dds_copy)
# colDat <- colData(dds_copy) %>% as.data.frame()
# idx_samples <- colDat %>% dplyr::filter(condition %in% contrast_as_vect) %>% rownames_to_column("samples.id") %>% pull("samples.id")
# #idx_samples <- c(idx_samples, "Ctrl_4")
# dds_copy_subset <- dds_copy[,idx_samples]
# dds_copy_subset <- DESeq(dds_copy_subset)
# contrast_input <- "condition_AZCT_vs_DMSO" #uncomment to test run function

# metadata(dds_copy)
# elementMetadata(dds_copy)
# dseqObject <- dds_copy[,-idx]
################################################################################################################################################
######################################################################## LRT geneclustering ###################################################
#prepare data for clustering
dds2 <- DESeq(dds, test = "LRT", reduced = ~1)
names(mcols(dds2))
acrossGroups <- results(dds2)
acrossGroups <- acrossGroups[order(acrossGroups$pvalue), ]
acrossGroups[1:3, ]

rld_4dds2 <- rlog(dds2, blind=FALSE)
rld_4dds2_mat <- assay(rld_4dds2)

sigChanges <- rownames(acrossGroups)[acrossGroups$padj < 0.01 & !is.na(acrossGroups$padj)]
sigMat <- rld_4dds2_mat[rownames(rld_4dds2_mat) %in% sigChanges, ]
pheatmap(sigMat, scale = "row", show_rownames = FALSE, annotation = metadata_rown_df[,c("condition"), drop = FALSE],main = "results_LRT gene clustering",filename = "figures/geneClustering_LRT_full_vrs_reducedModel.pdf")
k <- pheatmap(sigMat, scale = "row", kmeans_k = 4, annotation = metadata_rown_df[,c("condition"), drop = FALSE],main = "results_LRT gene clustering; k=4",filename = "figures/geneClustering_LRT_full_vrs_reducedModel_kmeans4.pdf")

# names(k$kmeans)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
# clusterDF[1:10, , drop = FALSE]
# head(clusterDF)
OrderByCluster <- sigMat[order(clusterDF$Cluster), ]

pheatmap(OrderByCluster, scale = "row", annotation_row = clusterDF, show_rownames = FALSE,
    cluster_rows = FALSE,
    annotation = metadata_rown_df[,c("condition"), drop = FALSE], 
    # annotation_colors = ann_colors,
    filename = "figures/geneClustering_LRT_full_vrs_reducedModel_kmeans4_withClutersAnnot.pdf")


# find the optimal clusters
heatmap_matrixAdjusted <- cola::adjust_matrix(
    sigMat,
    sd_quantile = 0.50,
    max_na = 0.25,
    verbose = TRUE)
rowScaledMat_cola <- t(scale(t(heatmap_matrixAdjusted)))


rowScaledMat <- t(scale(t(sigMat)))
set.seed(123) #for reproducibility
clusterNum <- NbClust(rowScaledMat, distance = "euclidean", min.nc = 2, max.nc = 12,
    method = "kmeans", index = "silhouette")
clusterNum$Best.nc

library(factoextra)
set.seed(342)
pdf("figures/optimal_clusters_silhouette.pdf")
fviz_nbclust(
  rowScaledMat_cola,
  FUNcluster = kmeans,
  method     = "silhouette",
  k.max      = 12,
  nstart     = 25,
  iter.max   = 100
)
dev.off()

clusterNumCola <- NbClust(rowScaledMat_cola, distance = "euclidean", min.nc = 2, max.nc = 12,
    method = "kmeans", index = "silhouette")
clusterNumCola$Best.nc

# merge automatic clustering with manual clustering
clusterNumCola$Best.partition[1:10]
orderedCluster <- sort(clusterNumCola$Best.partition)
rowScaledMat_cola <- rowScaledMat_cola[match(names(orderedCluster), rownames(rowScaledMat_cola)), ]
# length(clusterNum$Best.partition)

# annotation_row = clusterDF
pheatmap(rowScaledMat_cola, scale = "row",
        cluster_rows = TRUE, show_rownames = FALSE,
    annotation = metadata_rown_df[,c("condition"), drop = FALSE], 
    filename = "figures/geneClustering_LRT_full_vrs_reducedModel_kmeans_withOptimumClustering.pdf")
################################################################################################################################################
# nrow(rld_4dds2_mat)
#defining the optimal clusters

### save the dds object for alll conditions
#saveRDS(dds, file = paste0("data/dds_coding_genes_",paste0(drop_samples, collapse="_"),"_dropped.rds"))




#####################################################################################################################
############################### Step 3; Run deseq2 analysis for each contrast
#####################################################################################################################
DEResults_ls <- list()
################################################################################################
############################### Step #; Run DE for each contrast
##############################################################################################
#run for all contrasts
#contrast_input="condition_AZA_vs_DMSO"
# DEResults_ls <- lapply(contrasts_ls[1:2], function(x) {run_multiple_contrasts(contrast_input=x, dseqObject = dds, shrink = FALSE, export_transformed_counts = TRUE)})

# names(DEResults_ls) <- contrasts_ls[1:2]

DEResults_ls <- lapply(contrasts_ls, function(x) {run_multiple_contrasts(contrast_input=x, dseqObject = dds, shrink = FALSE, export_transformed_counts = TRUE)})
names(DEResults_ls) <- contrasts_ls
# names(DEResults_ls[[1]])
#summary(DEResults_ls[[1]])
# data.frame(DEResults_ls[[1]]) %>% dplyr::filter(padj < 0.05) %>% dim()
# contrast_input=contrasts_ls[1]









################################################################################################
############################### Step #; GSEA
##############################################################################################

# #cnet plots
# ## convert gene ID to Symbol
# edox <- setReadable(msig_GSEA_obj, 'org.Hs.eg.db', 'ENTREZID')
# p1 <- cnetplot(edox, foldChange=foldchanges)
# ## categorySize can be scaled by 'pvalue' or 'geneNum'
# p2 <- cnetplot(edox, categorySize="pvalue", foldChange=foldchanges)
# # p3 <- cnetplot(edox, foldChange=foldchanges, circular = TRUE, colorEdge = TRUE) 
# #cnetPlot_cowp <- cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])
# ggsave(p1, file = paste0("figures/cnetplot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,".pdf"), width = ggwidth, height = ggheight)

# gseaObj_name <- paste0("msig_category_",category_tag) #what name should be given to the gsea object in list

# if(do_GSEA){
# out_gsea_fun <- list(rankedFC=foldchanges, gseaObj_name = msig_GSEA_obj)
# saveRDS(out_gsea_fun, file=paste0("data/GSEA_msigdbr_cat_",category_tag,"_",name_de_res,".rds"))
# }else if(do_gse_GO){
# out_gse <- list(rankedFC=foldchanges, gseGO=gse)
# saveRDS(out_gse, file=paste0("data/GSE_GO_",name_de_res,".rds"))
# }else(do_gse_GO & do_GSEA){
#     message("saving GSE & GSEA enrichment distribution plots")
#     #make list of results
# out_gsea_fun <- list(rankedFC=foldchanges, gseGO=gse, gseaObj_name = msig_GSEA_obj)
# saveRDS(out_gsea_fun, file=paste0("data/GSE_GO_plus_msigdbr_cat_",category_tag,"_",name_de_res,".rds"))
# return(out_gsea_fun)
# }

# return(out_gsea_fun)
#}

#paste0("figures/gsea_msigdbr_H_type_CKI_DMSO_EnrichDistrib.pdf")
#run gsea for all contrasts
# lapply(names_contrasts_ls, function(x) gsea_analysis(x))

source(paste0(workflow_dir,"scripts/run_GSEA_conditionSpecific.R"))

# dfres <- read_tsv("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/codingGenes_DE/data/condition_AZA_vs_DMSO/Dseq2Results_condition_AZA_vs_DMSO.tsv")
# h_gsea <- gsea_analysis("condition_AZA_vs_DMSO", results_df_with_annot=dfres, specie_type="Mus musculus",category_tag = "H", nCategory_2show = 50, ggwidth = 15, ggheight = 12)


# DEResults_ls[[x]]$results_with_Annot %>% filter(LFClabel %in% c("up", "down"))  %>% dplyr::filter(LFClabel %in% c("up", "down"))
h_gsea <- lapply(contrasts_ls, function(x) gsea_analysis(x, results_df_with_annot=DEResults_ls[[x]]$results_with_Annot, specie_type="Mus musculus", qval_filter = 0.05, category_tag = "MH",  subcategory_tag = NULL, nCategory_2show = 20, ggwidth = 15, ggheight = 12))
c2_gsea <- lapply(contrasts_ls, function(x) gsea_analysis(x, results_df_with_annot=DEResults_ls[[x]]$results_with_Annot, specie_type="Mus musculus",do_gse_GO=FALSE, nCategory_2show = 20, qval_filter = 0.05, category_tag = "C2",  subcategory_tag = NULL, ggwidth = 15, ggheight = 12))
names(h_gsea)

#save GSEA results
twosaveGsea <- list(h_gsea = h_gsea, c2_gsea = c2_gsea) 
save(twosaveGsea, file = paste0("data/gsea_results_coding_genes.RData"))

if(!file.exists("data/gsea_results_coding_genes.RData")){
    message("GSEA results not saved")
    message("Running GSEA analysis for select contrasts")

ls_results <- c("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/codingGenes_DE/data//condition_CKi_vs_DMSO/Dseq2Results_condition_CKi_vs_DMSO.tsv",                                 
"/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/codingGenes_DE/data//condition_QSTAT_vs_DMSO/Dseq2Results_condition_QSTAT_vs_DMSO.tsv",                             
"/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/codingGenes_DE/data//condition_QSTAT-CKi_vs_DMSO/Dseq2Results_condition_QSTAT-CKi_vs_DMSO.tsv")
lapply(ls_results, function(x) {
    conds <- basename(dirname(x))
  dfRes <- read_tsv(x)
  print(head(dfRes))
  c2_gsea <- gsea_analysis(results_df_with_annot=dfRes, name_de_res = conds, specie_type="Mus musculus", do_gse_GO=FALSE, qval_filter = 0.05, nCategory_2show = 20, category_tag = "C2", subcategory_tag = NULL,ggwidth = 15, ggheight = 12)  
  h1_gsea <- gsea_analysis(results_df_with_annot=dfRes, name_de_res = conds, specie_type="Mus musculus", do_gse_GO=FALSE, qval_filter = 0.05, nCategory_2show = 20, category_tag = "H", subcategory_tag = NULL,ggwidth = 15, ggheight = 12)
})
} else{
    message("GSEA results saved to data/gsea_results_coding_genes.RData")
}


# names(h_gsea[[1]][["rankedFC"]])

# c2_gsea[[1]]
# head(h_gsea[[1]]@result)
# slotnames(h_gsea)
# slotNames(h_gsea)



### Add volcano plots