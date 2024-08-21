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

# library(devtools)
# install_github("ropensci/magick")

results_dir <- "/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/sandbox/data/"
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
# duplicated(dt_cbind_LFCup)

# sum(duplicated(dt_cbind))
# table(dt_cbind_LFCup$`SETDB1i_vs_QSTAT-CKi.gene.id`)
# dt_cbind_LFCup %>% group_by(`SETDB1i_vs_QSTAT-CKi.gene.id`) %>% summarize(n=n())

# dt_cbind
# duplicated_rows_dplyr <- dt_cbind_LFCup %>%
#   group_by_all() %>%
#   dplyr::filter(n() > 1) %>%
#   ungroup()
# Reduce(function(d1, d2) cbind(d1, d2), upLFC_all_geneIDs_only)

# list_to_matrix(upLFC_all_geneIDs_only)


# map_deseq_results[[1]][,"gene.id" ]
# map_deseq_results[[1]] %>% dplyr::filter(padj < 0.05 & log2FoldChange > 1)
# lapply(upLFC_all_geneIDs_genes_only_dt, dim)

names_files[grep("vs_DMSO",names_files)]
names_files[grep("vs_DMSO",names_files)]
names_files[grep("CKi_vs",names_files)]

combinedTherapy <- c("SETDB1i-CKi_vs_SETDB1i" ,"SETDB1i-CKi_vs_CKi","QSTAT-CKi_vs_QSTAT" ,"QSTAT-CKi_vs_CKi","SETDB1i-CKi_vs_QSTAT-CKi")
monTherapy <- c("AZA_vs_DMSO", "CKi_vs_DMSO","QSTAT_vs_DMSO", "QSTAT-CKi_vs_DMSO", "SETDB1i_vs_DMSO", "SETDB1i-CKi_vs_DMSO")

list_drugs_vs_dmso <- names_files[grepl("vs_DMSO",names_files)]
# list_drugs_vs_dmso <- names(sigLFC_all_genes_only_wide_df)[grepl("vs_DMSO",names(sigLFC_all_genes_only_wide_df))]



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
# ggsave(upset_out, filename = "upset_all_conditions_genes_LFC_above1_padj_0.05.png",width = 10, height = 10)

#easy on the eyes for presentation
upset_monoCombo_curated <- upset(intersections_all_conds, order.by = "freq", sets = list_drugs_vs_dmso, set_size.show = TRUE)
pdf("upset_mono_combination_coding_genes_LFC_above1_padj_005_curated.pdf",width = 18, height = 10)
print(upset_monoCombo_curated)
dev.off()


################################################################################################################
#clustering of gebes with significant LFC
################################################################################################################


sigLFC_all <- lapply(map_deseq_results, function(x){x[!is.na(padj)][padj<0.05 & abs(log2FoldChange) > 1, ]})
sigLFC_all[[1]]
sigLFC_all_genes_only_dt <- lapply(sigLFC_all, function(x){x[,c("gene.id","log2FoldChange")]})
names(sigLFC_all_genes_only_dt) <- names_files
sigLFC_all_genes_only_wide_df <- rbindlist(sigLFC_all_genes_only_dt, fill=T, idcol = "condition") %>% 
  dcast(gene.id ~ condition, value.var = "log2FoldChange") %>% 
  as.data.frame() %>% 
  column_to_rownames("gene.id") 

#dmso is ref so remove cases when used as alt
sigLFC_all_genes_only_wide_df <- sigLFC_all_genes_only_wide_df %>% select(!starts_with("DMSO_vs")) #%>% head()



sigLFC_all_genes_only_wide_df_matrix <- data.matrix(sigLFC_all_genes_only_wide_df)
sigLFC_all_genes_only_wide_df_matrix[is.na(sigLFC_all_genes_only_wide_df_matrix)] = 0

#seperate this plot
matrixMonoCombo_DMSO_ref <- sigLFC_all_genes_only_wide_df_matrix[, list_drugs_vs_dmso]
matrixMonoCombo_DMSO_ref[is.na(matrixMonoCombo_DMSO_ref)] = 0
matrixMonoCombo_DMSO_ref_cleaned <- matrixMonoCombo_DMSO_ref[rowSums(matrixMonoCombo_DMSO_ref) > 0,]

constrasts_df <- data.frame(contrast = names(sigLFC_all_genes_only_wide_df))
DrugsAnnot <- constrasts_df %>% dplyr::filter(contrast %in% list_drugs_vs_dmso)%>% mutate(Combinations = case_when(str_detect(contrast, "-") ~ "combo",TRUE ~ "mono"), BroadTargets = case_when(str_detect(contrast, "QSTAT_vs|SETDB1i_vs") ~ "Chromatin", str_detect(contrast, "-CKi_vs_") ~ "Chromatin+MEK" ,str_detect(contrast, "AZA_") ~ "DNAme", str_detect(contrast, "CKi_vs_DMSO") ~ "MEK", TRUE ~ NA), Action = case_when(str_detect(contrast, "QSTAT_vs") ~ "HDACi", str_detect(contrast, "SETDB1i_vs") ~ "H3K9mei" ,str_detect(contrast, "QSTAT-CKi_vs_DMSO") ~ "HDACi+MEKi", str_detect(contrast, "SETDB1i-CKi_vs_DMSO") ~ "H3K9mei+MEKi",str_detect(contrast, "AZA_") ~ "DNMTi", str_detect(contrast, "CKi_vs_DMSO") ~ "MEKi", TRUE ~ NA))
# )
rownames(DrugsAnnot) <- DrugsAnnot$contrast
colnames(matrixMonoCombo_DMSO_ref_cleaned)
# sum(complete.cases(sigLFC_all_genes_only_wide_df_matrix))

pdf(file = paste0("clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_allContrasts.pdf"), width=12, height=9)
pheatmap(sigLFC_all_genes_only_wide_df_matrix, main = paste0("abs(LFC) > 1 & padj < 0.05; consensus coding genes[NA=0]"), scale = "row", annotation_names_row=FALSE, show_rownames = F)
dev.off()

# p(file = paste0("clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.pdf")) #, width=12, height=9
pheatmap::pheatmap(matrixMonoCombo_DMSO_ref_cleaned,annotation = DrugsAnnot[,c("Combinations", "Action")], main = paste0("abs(LFC) > 1 & padj < 0.05; consensus coding genes[NA=0]"), scale = "row", annotation_names_row=FALSE, show_rownames = F, colorRampPalette(c("#5D3A9B", "white", "#E66100"))(100), filename = paste0("clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.pdf"))
# dev.off()


# png(file = paste0("Kmeans_clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.png")) #, width=12, height=9
k <- pheatmap::pheatmap(matrixMonoCombo_DMSO_ref_cleaned, kmeans_k = 3, annotation = DrugsAnnot[,c("Combinations", "Action")], colorRampPalette(c("#5D3A9B", "white", "#E66100"))(100),main = paste0("abs(LFC) > 1 & padj < 0.05; consensus coding genes[NA=0]"), scale = "row", annotation_names_row=FALSE, show_rownames = F, filename = paste0("Kmeans_clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.png"))


# names(k$kmeans)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
OrderByCluster <- matrixMonoCombo_DMSO_ref_cleaned[order(clusterDF$Cluster), ]
pheatmap::pheatmap(OrderByCluster, annotation = DrugsAnnot[,c("Combinations", "Action")], scale = "row", annotation_row = clusterDF, show_rownames = FALSE,
colorRampPalette(c("#5D3A9B", "white", "#E66100"))(100),
    cluster_rows = FALSE, filename = paste0("orderedKmeans_clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.png"))


matrixMonoCombo_DMSO_ref_cleaned_scaled <- t(scale(t(matrixMonoCombo_DMSO_ref_cleaned)))
clusterNum <- NbClust(matrixMonoCombo_DMSO_ref_cleaned_scaled, distance = "euclidean", min.nc = 2, max.nc = 12,method = "kmeans", index = "silhouette")
message("optimal number of clusters: ", clusterNum$Best.nc[1], "; with value index: ", clusterNum$Best.nc[2])

clusterDF$`gene.id`` <- rownames(clusterDF)
write.table(clusterDF, file = "orderedKmeans_clustering_sig_coding_DE_NA_set_to_0_mono_combination_therapy_DMSO_ref.tsv", sep="\t", row.names = FALSE, col.names = TRUE)



# k <- pheatmap(sigMat, scale = "row", kmeans_k = 7, annotation = metadata_rown_df[,c("condition"), drop = FALSE] ,filename = paste0("figures/", contrast_input,"/geneClustering_LRT_full_vrs_reducedModel_kmeans7_",contrast_input,"_only.pdf"))

# #subset on combination therapies
# pdf(file = paste0("clustering_sig_DE_NA_set_to_0_mono_combination_therapy.pdf"), width=12, height=9)
# pheatmap(sigLFC_all_genes_only_wide_df_matrix, main = paste0("abs(LFC) > & padj < 0.05 [coding genes]"), scale = "row", annotation_names_row=FALSE, show_rownames = F)
# dev.off()