#author; samuel ahuno
# combine multiple conditions deseq2 results
#multiple conditions
library("data.table")
library(tidyverse)
library(UpSetR)
library("ComplexHeatmap")
library(rlist)
library(dplyr)
library(NbClust)
library(GSVA)
library(optparse)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
require(DOSE)
library(msigdbr)
# library(devtools)
# install_github("ropensci/magick")

# # RUN THIS SCRIPT WITH THE FOLLOWING COMMAND
# Rscript /data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/scripts/combine_results_of_mutiple_conditons.R \
#   --specie_type "Mus musculus" \
#   --category_tag "C2" \
#   --fractEnriched 0.5 \
#   --do_gse_GO \
#   --nCategory_2show 10 \
#   --nShowBar 100 \
#   --ggwidth 13 \
#   --ggheight 11 \
#   --results_dir "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/codingGenesDExpr/data"
  
# Parse command line options
library(optparse)
option_list <- list(
  make_option(c("-s", "--specie_type"), type="character", default="Mus musculus",
              help="Species type, e.g., 'Homo sapiens' or 'Mus musculus'"),
  make_option(c("-c", "--category_tag"), type="character", default="C2",
              help="Category tag for enrichment analysis"),
  make_option(c("-f", "--fractEnriched"), type="double", default=0.5,
              help="Fraction enriched threshold"),
  make_option(c("-g", "--do_gse_GO"), action="store_true", default=TRUE,
              help="Whether to perform GSEA GO"),
  make_option(c("-n", "--nCategory_2show"), type="integer", default=10,
              help="Number of categories to show"),
  make_option(c("-b", "--nShowBar"), type="integer", default=100,
              help="Number of bars to show in plots"),
  make_option(c("--ggwidth"), type="double", default=13,
              help="Width of ggplot figures"),
  make_option(c("--ggheight"), type="double", default=11,
              help="Height of ggplot figures"),
  make_option(c("-r", "--results_dir"), type="character",
              default="data/",
              help="Results directory path")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Assign variables from parsed options
specie_type    <- opt$specie_type
category_tag   <- opt$category_tag
fractEnriched  <- opt$fractEnriched
do_gse_GO      <- opt$do_gse_GO
nCategory_2show<- opt$nCategory_2show
nShowBar       <- opt$nShowBar
ggwidth        <- opt$ggwidth
ggheight       <- opt$ggheight
results_dir    <- opt$results_dir

files_deseq_results1 <- list.files(results_dir, recursive = TRUE, full.names = TRUE, pattern = "Dseq2Results_") 
files_deseq_results <- files_deseq_results1[-grep("only_LTR_Test.tsv", files_deseq_results1)]

map_deseq_results <- map(files_deseq_results, ~fread(.x))
names_files <- sapply(files_deseq_results, function(x){gsub("Dseq2Results_condition_|.tsv", "", basename(x))})
upLFC_all <- lapply(map_deseq_results, function(x){x[!is.na(padj)][padj<0.05 & log2FoldChange > 1, ]})
# lapply(upLFC_all, dim)
upLFC_all_geneIDs_genes_only_dt <- lapply(upLFC_all, function(x){x[,"gene.id" ]})
upLFC_all_geneIDs_genes_only_vector <- lapply(upLFC_all, function(x){x[,"gene.id" ] %>% pull(gene.id)})

names(upLFC_all_geneIDs_genes_only_dt) <- names_files
names(upLFC_all_geneIDs_genes_only_vector) <- names_files
names(names_files) <- NULL #remove long file paths
dt_cbind_LFCup <- list.cbind(upLFC_all_geneIDs_genes_only_dt)

names_files[grep("vs_DMSO",names_files)]
names_files[grep("vs_DMSO",names_files)]
names_files[grep("CKi_vs",names_files)]

combinedTherapy <- c("SETDB1i-CKi_vs_SETDB1i" ,"SETDB1i-CKi_vs_CKi","QSTAT-CKi_vs_QSTAT" ,"QSTAT-CKi_vs_CKi","SETDB1i-CKi_vs_QSTAT-CKi")
monTherapy <- c("AZA_vs_DMSO", "CKi_vs_DMSO","QSTAT_vs_DMSO", "QSTAT-CKi_vs_DMSO", "SETDB1i_vs_DMSO", "SETDB1i-CKi_vs_DMSO")

list_drugs_vs_dmso <- names_files[grepl("vs_DMSO",names_files)]

intersections_all_conds <- fromList(upLFC_all_geneIDs_genes_only_vector)

upset_out <- upset(intersections_all_conds, order.by = "freq", nsets = 12)
pdf("upset_all_conditions_coding_genes_LFC_above1_padj_005.pdf",width = 10, height = 10)
print(upset_out)
dev.off()

#plot alll combo for supplementary
upset_monoCombo <- upset(intersections_all_conds, order.by = "freq", sets = c(monTherapy,combinedTherapy), set_size.show = TRUE)
pdf("upset_mono_combination_coding_genes_LFC_above1_padj_005_v2.pdf",width = 18, height = 10)
print(upset_monoCombo)
dev.off()

#easy on the eyes for presentation
upset_monoCombo_curated <- upset(intersections_all_conds, order.by = "freq", sets = list_drugs_vs_dmso, set_size.show = TRUE)
pdf("upset_mono_combination_coding_genes_LFC_above1_padj_005_curated.pdf",width = 18, height = 10)
print(upset_monoCombo_curated)
dev.off()


################################################################################################################
#clustering of gebes with significant LFC
################################################################################################################
sigLFC_all <- lapply(map_deseq_results, function(x){x[!is.na(padj)][padj<0.05 & abs(log2FoldChange) > 1, ]})

sigLFC_all_genes_only_dt <- lapply(sigLFC_all, function(x){x[,c("gene.id","log2FoldChange")]})
names(sigLFC_all_genes_only_dt) <- names_files
sigLFC_all_genes_only_wide_df <- rbindlist(sigLFC_all_genes_only_dt, fill=T, idcol = "condition") %>% 
  dcast(gene.id ~ condition, value.var = "log2FoldChange") %>% 
  as.data.frame() %>% 
  column_to_rownames("gene.id") 

#dmso is ref so remove cases where it was used as alt
sigLFC_all_genes_only_wide_df <- sigLFC_all_genes_only_wide_df %>% select(!starts_with("DMSO_vs")) #%>% head()


#remove NA to make it easy for plotting 
sigLFC_all_genes_only_wide_df_matrix <- data.matrix(sigLFC_all_genes_only_wide_df)
sigLFC_all_genes_only_wide_df_matrix[is.na(sigLFC_all_genes_only_wide_df_matrix)] = 0

#seperate this plot
matrixMonoCombo_DMSO_ref <- sigLFC_all_genes_only_wide_df_matrix[, list_drugs_vs_dmso]
matrixMonoCombo_DMSO_ref[is.na(matrixMonoCombo_DMSO_ref)] = 0
matrixMonoCombo_DMSO_ref_cleaned <- matrixMonoCombo_DMSO_ref[rowSums(matrixMonoCombo_DMSO_ref) > 0,]

constrasts_df <- data.frame(contrast = names(sigLFC_all_genes_only_wide_df))
DrugsAnnot <- constrasts_df %>% dplyr::filter(contrast %in% list_drugs_vs_dmso)%>% mutate(Combinations = case_when(str_detect(contrast, "-") ~ "combo",TRUE ~ "mono"), BroadTargets = case_when(str_detect(contrast, "QSTAT_vs|SETDB1i_vs") ~ "Chromatin", str_detect(contrast, "-CKi_vs_") ~ "Chromatin+MEK" ,str_detect(contrast, "AZA_") ~ "DNAme", str_detect(contrast, "CKi_vs_DMSO") ~ "MEK", TRUE ~ NA), Action = case_when(str_detect(contrast, "QSTAT_vs") ~ "HDACi", str_detect(contrast, "SETDB1i_vs") ~ "H3K9mei" ,str_detect(contrast, "QSTAT-CKi_vs_DMSO") ~ "HDACi+MEKi", str_detect(contrast, "SETDB1i-CKi_vs_DMSO") ~ "H3K9mei+MEKi",str_detect(contrast, "AZA_") ~ "DNMTi", str_detect(contrast, "CKi_vs_DMSO") ~ "MEKi", TRUE ~ NA))
DrugsAnnot <- DrugsAnnot %>% mutate(MEKi = case_when(str_detect(Action, "MEKi") ~ "Yes", TRUE ~ "No"))
# )
rownames(DrugsAnnot) <- DrugsAnnot$contrast
colnames(matrixMonoCombo_DMSO_ref_cleaned)
# sum(complete.cases(sigLFC_all_genes_only_wide_df_matrix))

pdf(file = paste0("clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_allContrasts.pdf"), width=12, height=9)
pheatmap(sigLFC_all_genes_only_wide_df_matrix, main = paste0("abs(LFC) > 1 & padj < 0.05; consensus coding genes[NA=0]"), scale = "row", annotation_names_row=FALSE, show_rownames = F)
dev.off()

# p(file = paste0("clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.pdf")) #, width=12, height=9
pheatmap::pheatmap(matrixMonoCombo_DMSO_ref_cleaned,annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")], main = paste0("abs(LFC) > 1 & padj < 0.05; consensus coding genes[NA=0]"), scale = "row", annotation_names_row=FALSE, show_rownames = F, colorRampPalette(c("#5D3A9B", "white", "#E66100"))(100), filename = paste0("clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.pdf"))
# dev.off()


# png(file = paste0("Kmeans_clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.png")) #, width=12, height=9
k <- pheatmap::pheatmap(matrixMonoCombo_DMSO_ref_cleaned, kmeans_k = 3, annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")], colorRampPalette(c("#5D3A9B", "white", "#E66100"))(100),main = paste0("abs(LFC) > 1 & padj < 0.05; consensus coding genes[NA=0]"), scale = "row", annotation_names_row=FALSE, show_rownames = F, filename = paste0("Kmeans_clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.png"))


# names(k$kmeans)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
OrderByCluster <- matrixMonoCombo_DMSO_ref_cleaned[order(clusterDF$Cluster), ]
pheatmap::pheatmap(OrderByCluster, annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")], scale = "row", annotation_row = clusterDF, show_rownames = FALSE,
colorRampPalette(c("#5D3A9B", "white", "#E66100"))(100),
    cluster_rows = FALSE, filename = paste0("orderedKmeans_clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.png"))


matrixMonoCombo_DMSO_ref_cleaned_scaled <- t(scale(t(matrixMonoCombo_DMSO_ref_cleaned)))
clusterNum <- NbClust(matrixMonoCombo_DMSO_ref_cleaned_scaled, distance = "euclidean", min.nc = 2, max.nc = 12, method = "kmeans", index = "silhouette")
message("optimal number of clusters: ", clusterNum$Best.nc[1], "; with value index: ", clusterNum$Best.nc[2])

clusterDF$`gene.id` <- rownames(clusterDF)
write.table(clusterDF, file = "orderedKmeans_clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.tsv", sep="\t", row.names = FALSE, col.names = TRUE)

#merge nbclusters with metadata
# clusterNum$Best.partition[1:10]
orderedCluster <- sort(clusterNum$Best.partition)
matrixMonoCombo_DMSO_ref_cleaned_scaled_mergedNBClust <- matrixMonoCombo_DMSO_ref_cleaned_scaled[match(names(orderedCluster), rownames(matrixMonoCombo_DMSO_ref_cleaned_scaled)), ]

pheatmap::pheatmap(matrixMonoCombo_DMSO_ref_cleaned_scaled_mergedNBClust, scale = "row", annotation_row = clusterDF, show_rownames = FALSE,
    cluster_rows = FALSE,
    annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")], 
    filename = paste0("kmeans_withOptimumClustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.pdf"))

# heatmapA <- DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")]
# column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))

# pdf(paste0("kmeans_withOptimumClustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref_heatmap_v2.pdf"))
# ComplexHeatmap::Heatmap(matrixMonoCombo_DMSO_ref_cleaned_scaled_mergedNBClust, name = "mat",
# # right_annotation = DrugsAnnot[,c("Combinations", "Action", "BroadTargets", "MEKi")] , 
# cluster_rows = FALSE, cluster_columns = TRUE, use_raster = TRUE, show_row_names = FALSE, 
# show_column_names = TRUE, col = colorRampPalette(c("#5D3A9B", "white", "#E66100"))(100))
# dev.off()

optim_orderedClusterDF <- data.frame(clusterid=factor(orderedCluster))

optim_orderedClusterDF$`gene.id` <- names(orderedCluster)
write.table(optim_orderedClusterDF, file = "Optim_orderedKmeans_clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.tsv", sep="\t", row.names = FALSE, col.names = TRUE)


# k <- pheatmap(sigMat, scale = "row", kmeans_k = 7, annotation = metadata_rown_df[,c("condition"), drop = FALSE] ,filename = paste0("figures/", contrast_input,"/geneClustering_LRT_full_vrs_reducedModel_kmeans7_",contrast_input,"_only.pdf"))
splitCluster <- split(optim_orderedClusterDF, optim_orderedClusterDF$clusterid)
# splitCluster[[1]]
# #subset on combination therapies
# pdf(file = paste0("clustering_sig_DE_NA_set_to_0_mono_combination_therapy.pdf"), width=12, height=9)
# pheatmap(sigLFC_all_genes_only_wide_df_matrix, main = paste0("abs(LFC) > & padj < 0.05 [coding genes]"), scale = "row", annotation_names_row=FALSE, show_rownames = F)
# dev.off()


################################################################################################################
###################################################################################################################
####### gene set enrichment analysis ###### 
################################################################################################################
#######################
# GSEA GO; to identify the clusters of genes that are enriched in the gene sets
# splitCluster[[1]]

#i need global geneSet
uniqueGlobalGeneSet <- unique(map_deseq_results[[1]]$entrez.gene.id)
uniqueGlobalGeneSet2 <- uniqueGlobalGeneSet[!is.na(uniqueGlobalGeneSet)]

## use clusters as Background gene set 
Background_clusters <- map_deseq_results[[1]][!is.na(match(map_deseq_results[[1]]$gene.id, optim_orderedClusterDF$`gene.id`)),]
Background_clusters_final <- Background_clusters[!is.na(Background_clusters$entrez.gene.id)]

AnnotClusters <- function(dataIn){
  # sigLFC_all[[1]]$entrez.gene.id
Keys_sample1 <- map_deseq_results[[1]][!is.na(match(map_deseq_results[[1]]$gene.id, dataIn$`gene.id`)),]
Keys_sample2 <- Keys_sample1[!is.na(Keys_sample1$entrez.gene.id)]

# optim_orderedClusterDF$`gene.id`
uniqueEntrezIDs <- unique(Keys_sample2$entrez.gene.id)
if(do_gse_GO){
if(specie_type == "Homo sapiens"){
    ego <- enrichGO(gene      = uniqueEntrezIDs,
                universe      = Background_clusters_final,
                ont           = "BP",
                pAdjustMethod = "BH",
                OrgDb = "org.Hs.eg.db", 
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,readable= TRUE)    

}else if(specie_type == "Mus musculus"){
  ego <- enrichGO(gene        = uniqueEntrezIDs,
                universe      = Background_clusters_final,
                ont           = "BP",
                pAdjustMethod = "BH",
                OrgDb = "org.Mm.eg.db", 
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,readable= TRUE)
}
return(ego)
}
}

#run annotation of clusters
# AnnotClusters(dataIn=splitCluster[[1]])
AnnotattedClusters <- lapply(splitCluster, AnnotClusters)
AnnotattedClusters_DT_list <- lapply(AnnotattedClusters,function(x){as.data.table(x)} )
# head(ego, 20) %>% arrange(qvalue) %>% 
cbindlist(AnnotattedClusters_DT_list)

AnnotattedClusters_DT <- rbindlist(AnnotattedClusters_DT_list,fill=FALSE, idcol="clusterID")
AnnotattedClusters_DT[,ID_GO:=paste(ID, Description)]
AnnotattedClusters_DT[,ID_GO:=str_wrap(ID_GO, width = 10)]

AnnotattedClusters_DT_2plot <- AnnotattedClusters_DT %>% group_by(clusterID) %>% slice_head(n = 10)

GO_bP_plot <-  ggplot(data=AnnotattedClusters_DT_2plot, aes(x = -log10(qvalue), y=ID_GO)) + geom_col() + facet_wrap(~clusterID, scales = "free") + theme_minimal() + theme(axis.text.x=element_text(angle=90)) + labs(title = "GO BP enrichment")
# %>% write.table(file = "GSEA_GO_enrichedSets.tsv", sep="\t", row.names = FALSE, col.names = TRUE)
ggsave(GO_bP_plot, filename = "GO_OverRepresentation_BP_enrichedSets.pdf", width = 18, height = 10, units = "in", dpi = 300)

