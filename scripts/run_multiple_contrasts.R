
addLFC_Groups <- function(dataIn){
d <- dataIn %>% 
  dplyr::mutate(groups = case_when(
    padj < padj_val & log2FoldChange < -1 ~ "sigLFC_Dwn",
    padj < padj_val & log2FoldChange > 1 ~ "sigLFC_Up",
    padj > padj_val & log2FoldChange < -1 ~ "nonSigLFC_Dwn",
    padj > padj_val & log2FoldChange > 1 ~ "nonSigLFC_Up",
    padj < padj_val & abs(log2FoldChange) < 1 ~ "sigNoLFC",
    padj > padj_val & abs(log2FoldChange) < 1 ~ "nonSigNoLFC",
    TRUE ~ NA))
    return(d)
}



run_multiple_contrasts <- function(contrast_input, dseqObject = dds, shrink = FALSE, padj_val = 0.05, export_transformed_counts = FALSE, ntopVarGenes=20, foldChangeCutoff=1){
  message("analyzing ", contrast_input)
dir.create(paste0("data/", contrast_input, "/"), recursive = TRUE, showWarnings = TRUE)
dir.create(paste0("figures/", contrast_input, "/"), recursive = TRUE, showWarnings = TRUE)

#contrast_label <- paste(contrast_input, collapse = "_") #used for plot labels
str_vec <- gsub("_vs_", "_", contrast_input); 
contrast_as_vect <- unlist(str_split(str_vec, pattern="_"))

#print(contrast_as_vect)
##################################################################
#export Transform counts for data visualization
message("using constrast")
if(export_transformed_counts){
dseqObject_copy <- dseqObject
colDat <- colData(dseqObject_copy) %>% as.data.frame()
idx_samples <- colDat %>% dplyr::filter(condition %in% contrast_as_vect) %>% rownames_to_column("samples.id") %>% pull("samples.id")
dds_copy_subset <- dseqObject_copy[,idx_samples]
dds_copy_subset$condition <- droplevels(dds_copy_subset$condition)
# colnames(dds_copy_subset)
# names(mcols(dds_copy_subset))
#drop, genes with no padj
# dds_copy_subset <- dds_copy_subset[!is.na(dds_copy_subset$padj), ]

#relevel with last element in contrast vector
dds_copy_subset$condition <- relevel(dds_copy_subset$condition, ref = contrast_as_vect[3]) 
dds_copy_subset <- DESeq(dds_copy_subset) #call `DESeq` to reset the analysis with only the needed conditions and sample

# Transform counts for data visualization
rld <- rlog(dds_copy_subset, blind=TRUE)
rld_mat <- assay(rld)
rld_df <- rld_mat %>% as.data.frame() %>% rownames_to_column("gene.id")
#save a copy of the rlog matrix
write_tsv(rld_df, file=paste0("data/",contrast_input,"/DESeq2_rlog_Transform_BlindTRUE_",contrast_input,".tsv"))

#pairwise clustering
#prepare data for clustering

dds2 <- DESeq(dds_copy_subset, test = "LRT", reduced = ~1)
acrossGroups <- results(dds2)
acrossGroups <- acrossGroups[order(acrossGroups$pvalue), ]
rld_subsetSample <- rlog(dds2, blind=FALSE)
rld_subsetSample_mat <- assay(rld_subsetSample)

sigChanges <- rownames(acrossGroups)[acrossGroups$padj < 0.05 & !is.na(acrossGroups$padj)]
sigMat <- rld_subsetSample_mat[rownames(rld_subsetSample_mat) %in% sigChanges, ]
pheatmap(sigMat, scale = "row", show_rownames = FALSE, annotation = metadata_rown_df[,c("condition"), drop = FALSE] ,filename = paste0("figures/", contrast_input,"/geneClustering_LRT_full_vrs_reducedModel_",contrast_input,"_only.pdf"))
k <- pheatmap(sigMat, scale = "row", kmeans_k = 7, annotation = metadata_rown_df[,c("condition"), drop = FALSE] ,filename = paste0("figures/", contrast_input,"/geneClustering_LRT_full_vrs_reducedModel_kmeans7_",contrast_input,"_only.pdf"))

# names(k$kmeans)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
# clusterDF[1:10, , drop = FALSE]
# head(clusterDF)

OrderByCluster <- sigMat[order(clusterDF$Cluster), ]

# Specify colors
# ann_colors = list(Cluster = brewer.pal(n = length(levels(clusterDF$Cluster)), name = 'Dark2'))


pheatmap(OrderByCluster, scale = "row", annotation_row = clusterDF, show_rownames = FALSE,
    cluster_rows = FALSE,
    annotation = metadata_rown_df[,c("condition"), drop = FALSE], 
    # annotation_colors = ann_colors,
    filename = paste0("figures/", contrast_input,"/geneClustering_LRT_full_vrs_reducedModel_kmeans7_withClutersAnnot_",contrast_input,"_only.pdf"))

# find the optimal clusters
rowScaledMat <- t(scale(t(sigMat)))
clusterNum <- NbClust(rowScaledMat, distance = "euclidean", min.nc = 2, max.nc = 12,
    method = "kmeans", index = "silhouette")
message("optimal number of clusters: ", clusterNum$Best.nc[1], "; with value index: ", clusterNum$Best.nc[2])

#save results as tsv
DeSeq2Results_df_condOnly <- as.data.frame(acrossGroups) %>% rownames_to_column("gene.id")
#merge gene annotations and deseq results
DeSeq2Results_df_condOnly_annot <- left_join(DeSeq2Results_df_condOnly, annot %>% dplyr::select(c(gene.id, gene.symbol, description,  entrez.gene.id)), by="gene.id") %>% dplyr::filter(!is.na(padj)) 
#save deseq results with annotation as table
write.table(DeSeq2Results_df_condOnly_annot, file = paste0("data/",contrast_input,"/Dseq2Results_",contrast_input,"_only_LTR_Test.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)
}

# "figures/geneClustering_LRT_full_vrs_reducedModel_kmeans7.pdf"

##################################################################
#heatmaps of gene clustering
rld <- rlog(dseqObject, blind=FALSE) #use, dseq2 object with all samples/conditions/contrasts, doeen't matter since it's experiment-wide investigation
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), ntopVarGenes) #top variable genes across all samples
rld_top  <- assay(rld)[topVarGenes, ]
rld_top  <- rld_top - rowMeans(rld_top)
anno <- as.data.frame(colData(rld)[, c("samples","condition")])
row_annotations_df <- annot %>% dplyr::filter(gene.id %in% rownames(rld_top)) %>% dplyr::select(c(gene.id, gene.symbol, chr)) #%>% dplyr::mutate(gene.id = factor(gene.id, levels=rownames(rld_top)), rname = rownames(rld_top)) #order(.$gene.id)

#how to actually order the rows
row_annotations_df$gene.id <- factor(row_annotations_df$gene.id, levels=rownames(rld_top))
# %>% rownames_to_column("gene.id") %>% as.data.frame()
row_annotations_df <- row_annotations_df[order(row_annotations_df$gene.id)]
rownames(row_annotations_df) <- row_annotations_df$gene.id

#careful genes with multiple transcript might caause it to fail
if(sum(duplicated(row_annotations_df$gene.symbol)) < 1){
  message("no duplicated gene symbols, assign gene symbols to rownames")
  rownames(rld_top) <- row_annotations_df$gene.symbol
  rownames(row_annotations_df) <- row_annotations_df$gene.symbol

# heatmapGeneClustering <- pheatmap(rld_top, annotation_col = anno, annotation_row = row_annotations_df, annotation_names_row = TRUE)
heatmapGeneClustering <- pheatmap(rld_top, annotation_col = anno)
# ggsave(heatmapGeneClustering, file = paste0("figures/heatmap_genes_and_samples.pdf"), width = 12, height = 12)
ggsave(heatmapGeneClustering, file = paste0("figures/",contrast_input,"/heatmap_samples_and_top",ntopVarGenes,"genes.pdf"), width = 12 , height = 9)

}else{
heatmapGeneClustering <- pheatmap(rld_top, annotation_col = anno, annotation_row = row_annotations_df, annotation_names_row = TRUE)
# ggsave(heatmapGeneClustering, file = paste0("figures/heatmap_genes_and_samples.pdf"), width = 12, height = 12)
ggsave(heatmapGeneClustering, file = paste0("figures/",contrast_input,"/heatmap_samples_and_top",ntopVarGenes,"genes.pdf"), width = 12 , height = 9)
}


##################################################################

#run DESeq analysis - normalization and filtering
#use alpha of 0.05
results_per_contrast <- results(dseqObject, alpha=padj_val, contrast=contrast_as_vect)
#messsage(resultsNames(results_per_contrast))
#Shrink the log2 fold changes to be more accurate
if(shrink == TRUE){
  res <- lfcShrink(dds,
     coef = contrast_as_vect,
     type = "apeglm")
}
#contrast=contrast_as_vect
#save copy of each contrast results to disk
saveRDS(results_per_contrast, file=paste0("data/",contrast_input,"/Dseq2ResultsObject_",contrast_input,"_padjust.rds"))

DeSeq2Results_df <- as.data.frame(results_per_contrast) %>% rownames_to_column("gene.id")
#merge gene annotations and deseq results
DeSeq2Results_df_annot <- left_join(DeSeq2Results_df, annot %>% dplyr::select(c(gene.id, gene.symbol, description,  entrez.gene.id)), by="gene.id") %>% dplyr::filter(!is.na(padj)) 
#save deseq results with annotation as table
write.table(DeSeq2Results_df_annot, file = paste0("data/",contrast_input,"/Dseq2Results_",contrast_input,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)


#################################
###### manual volcano plots
# DeSeq2Results_df_annot_DEgroups %>% dplyr::group_by(groups) %>% dplyr::summarize(n=n())
DeSeq2Results_df_annot_DEgroups <- DeSeq2Results_df_annot %>% 
  dplyr::mutate(groups = case_when(
    padj < padj_val & log2FoldChange < -1 ~ "sigLFC_Dwn",
    padj < padj_val & log2FoldChange > 1 ~ "sigLFC_Up",
    padj > padj_val & log2FoldChange < -1 ~ "nonSigLFC_Dwn",
    padj > padj_val & log2FoldChange > 1 ~ "nonSigLFC_Up",
    padj < padj_val & abs(log2FoldChange) < 1 ~ "sigNoLFC",
    padj > padj_val & abs(log2FoldChange) < 1 ~ "nonSigNoLFC",
    TRUE ~ NA))
    # TRUE ~ "Unclassified"))
    # groups = paste0(groups, " (", n, ")"),
    # Counts = paste0(groups, "_", gene.symbol)),
    # Counts = factor(groups, levels = c("sigDown", "sigUp", "nonSigDown", "nonSigUp", "NoFC", "Unclassified"))
    # )

#grouped like enhanced volcano
dt_withGroupsEnhanced <- DeSeq2Results_df_annot %>% 
  dplyr::mutate(groups = case_when(
    padj < padj_val & abs(log2FoldChange) > 1 ~ "p-value and log2 FC",
    padj < padj_val & abs(log2FoldChange) < 1 ~ "p-value",
    padj > padj_val & abs(log2FoldChange) < 1 ~ "Log2 FC",
    TRUE ~ "NS"),
    Annot = case_when(
    padj < padj_val & abs(log2FoldChange) > 1 ~ "p-value and log2 FC",
    TRUE ~ "NS"),
    )

dt_withGroupsEnhanced <- dt_withGroupsEnhanced %>% add_count(Annot) %>% mutate(counts = paste0(Annot, " (", n, ")")) #%>% head()
dt_withGroupsEnhanced <- dt_withGroupsEnhanced %>% add_count(groups) %>% mutate(groups_counts = paste0(groups, " (", n, ")")) #%>% head()


# NS Log2 FC p−value p − value and log2 FC

DeSeq2Results_pvals <- DeSeq2Results_df_annot %>% 
  dplyr::mutate(groups = case_when(
    pvalue < padj_val & log2FoldChange < -1 ~ "sigLFC_Dwn",
    pvalue < padj_val & log2FoldChange > 1 ~ "sigLFC_Up",
    pvalue > padj_val & log2FoldChange < -1 ~ "nonSigLFC_Dwn",
    pvalue > padj_val & log2FoldChange > 1 ~ "nonSigLFC_Up",
    pvalue < padj_val & abs(log2FoldChange) < 1 ~ "sigNoLFC",
    pvalue > padj_val & abs(log2FoldChange) < 1 ~ "nonSigNoLFC",
    TRUE ~ NA))

DeSeq2Results_df_annot_DEgroups <- DeSeq2Results_df_annot_DEgroups %>% add_count(groups) %>% mutate(counts = paste0(groups, " (", n, ")")) #%>% head()
DeSeq2Results_pvals <- DeSeq2Results_pvals %>% add_count(groups) %>% mutate(counts = paste0(groups, " (", n, ")")) #%>% head()

nGroupsAssigned <- length(unique(DeSeq2Results_df_annot_DEgroups$groups))
nGroupsAssigned_pvals <- length(unique(DeSeq2Results_pvals$groups))

# nGroupsAssigned_pval <- length(unique(DeSeq2Results_df_annot_DEgroups$groups))


stats_df_padj <- DeSeq2Results_df_annot_DEgroups %>% dplyr::group_by(groups) %>% dplyr::summarize(n=n()) %>% ungroup() %>% mutate(groups = paste0(groups, " (", n, ")"), x = 5, y = seq.int(from=65,to=80,length.out=nGroupsAssigned))
stats_df_pval <- DeSeq2Results_pvals %>% dplyr::group_by(groups) %>% dplyr::summarize(n=n()) %>% ungroup() %>% mutate(groups = paste0(groups, " (", n, ")"), x = 5, y = seq.int(from=65,to=80,length.out=nGroupsAssigned_pvals))

# labels_padj <- paste(stats_df$groups, collapse = "\n")
labels_padj <- paste(stats_df_padj$groups, collapse = "; ")
labels_pval <- paste(stats_df_pval$groups, collapse = "; ")
# labels_pval <- paste(stats_df_pval$groups, collapse = "\n")
message("\nprinting labels_pval\n")
message(stats_df_pval)
message("\n ......... \n")
message("\n", names(dt_withGroupsEnhanced), "\n")


# table(DeSeq2Results_df_annot_DEgroups$groups)

plot_volcano <- ggplot(data=DeSeq2Results_df_annot_DEgroups) + 
geom_point(aes(x=log2FoldChange,y=-log10(padj), color=counts), alpha = 0.3) + 
geom_vline(xintercept = c(-1, 1), linetype="dotted") + 
geom_hline(yintercept = -log10(padj_val), linetype="dotted") +
theme_minimal() + 
labs(title = paste0("Volcano plot of ", contrast_input)) + 
theme(axis.text.x = element_text(angle = 45))
ggsave(plot_volcano, file = paste0("figures/",contrast_input,"/volcano_",contrast_input,"_LFC_padj_manual.pdf"), width = 9 , height = 7,dpi = 300)


plot_volcano_manual2 <- ggplot(data=dt_withGroupsEnhanced) + 
geom_point(aes(x=log2FoldChange,y=-log10(padj), color=groups_counts), alpha = 0.3) + 
geom_vline(xintercept = c(-1, 1), linetype="dotted") + 
geom_hline(yintercept = -log10(padj_val), linetype="dotted") +
theme_minimal() + 
labs(title = paste0("Volcano plot of ", contrast_input)) + 
theme(axis.text.x = element_text(angle = 45))
ggsave(plot_volcano_manual2, file = paste0("figures/",contrast_input,"/volcano_",contrast_input,"_LFC_padj_minimal_Annot_manual.pdf"), width = 9 , height = 7,dpi = 300)


significant_dt_withGroupsEnhanced <- dt_withGroupsEnhanced %>% dplyr::filter(padj < padj_val & abs(log2FoldChange) > foldChangeCutoff)
#signifcant and not significant genes annotations
plot_volcano_manual3 <- ggplot(data=dt_withGroupsEnhanced) + 
geom_point(aes(x=log2FoldChange,y=-log10(padj), color=counts), alpha = 0.3) + 
geom_vline(xintercept = c(-1, 1), linetype="dotted") + 
geom_hline(yintercept = -log10(padj_val), linetype="dotted") +
geom_label_repel(data=significant_dt_withGroupsEnhanced, aes(x=log2FoldChange, y=-log10(padj)), label=significant_dt_withGroupsEnhanced$`gene.symbol`, col='black', size=4) +
theme_minimal() + 
labs(title = paste0(contrast_input)) + 
scale_color_manual("Groups",values=c("grey","blue")) +
theme(axis.text.x = element_text(angle = 45)) + theme(legend.position = "top", legend.title = element_blank()) #
ggsave(plot_volcano_manual3, file = paste0("figures/",contrast_input,"/volcano_",contrast_input,"_LFC_padj_Sign_n_NonSig_manual.pdf"), width = 9 , height = 7,dpi = 300)



# DeSeq2Results_df_annot_DEgroups %>% dplyr::group_by(groups) %>% dplyr::summarize(n=n())
####do volcano plots
plt_volc_pval <- EnhancedVolcano(DeSeq2Results_df_annot,
    lab = DeSeq2Results_df_annot$gene.symbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    subtitleLabSize = 9,
    subtitle = labels_pval,
    pCutoff = padj_val,
    FCcutoff = foldChangeCutoff,
    title = paste0(contrast_input," (p < ", padj_val, ")"),
    pointSize = 1, labSize = 5,
    colAlpha = 0.2) +
    theme(plot.subtitle = element_blank())

# plt_volc_pval <- plt_volc_pval + labs(subtitle = labels_pval)
#annotate("text", x = 5, y = 60, label = labels_pval)
plt_volc_pval <- plt_volc_pval + geom_text(data=stats_df_pval, aes(x=x, y=y, label=groups), size=3)
ggsave(plt_volc_pval, file = paste0("figures/",contrast_input,"/volcano_",contrast_input,"_detailed_Annotation_pvals.pdf"))

#with adjusted p - values
#recreate data frame for ploting
DeSeq2Results_df_annot_padjustOnly <- DeSeq2Results_df_annot %>% dplyr::select(!pvalue) %>% mutate(pvalue = padj) 
plt_volc_padjust <- EnhancedVolcano(DeSeq2Results_df_annot_padjustOnly,
    lab = DeSeq2Results_df_annot_padjustOnly$gene.symbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    subtitle = labels_padj,
    subtitleLabSize = 9,
    pCutoff = padj_val,
    FCcutoff = foldChangeCutoff,
    title = paste0(contrast_input," (BH <", padj_val, ")"),
    widthConnectors = 0.75,
    labSize = 4.0,
    drawConnectors = TRUE,
    ylab = bquote('-' ~Log[10]~ 'p adjusted'),
    pointSize = 1, 
    colAlpha = 0.2) 
    # theme(plot.subtitle = labels_padj)
    # theme(plot.subtitle = element_blank())

message("\nplotting volcano plot with padjust\n")
message(labels_padj)
#add labels to the plot
# plt_volc_padjust <- plt_volc_padjust + labs(subtitle = labels_padj)

# plt_volc_padjust <- plt_volc_padjust + geom_text(data=stats_df_padj, aes(x=x, y=y, label=groups), size=3)

# labs(subtitle = labels_padj)
#annotate("text", x = 5, y = 60, label = labels_padj) #+ theme(plot.title = element_text(color="red", size=14, face="bold.italic"))
ggsave(plt_volc_padjust, file = paste0("figures/",contrast_input,"/volcano_",contrast_input,"_detailed_Annotation_padjust.pdf"))

pvals_df <- pivot_longer(DeSeq2Results_df_annot, cols = c(pvalue, padj)) #%>% head()
plt_pval  <- ggplot(data = pvals_df %>% filter(baseMean > 1)) + 
geom_histogram(aes(value)) + facet_wrap(~name) + 
theme_minimal() + 
labs(title = paste0("p & padj values histgram ", contrast_input)) + 
theme(axis.text.x = element_text(angle = 45))
ggsave(plt_pval, file = paste0("figures/",contrast_input,"/pNpadjustedValueHistogram_",contrast_input,".pdf"), width = 9 , height = 7,dpi = 300)

####Identifying significant genes
#Subset to return genes with padj < 0.05
sigLRT_genes <- DeSeq2Results_df_annot %>% 
  dplyr::filter(padj < padj_val)

sig_genes_merged_with_rlog <- sigLRT_genes %>% dplyr::select(c(gene.symbol, gene.id)) %>% left_join(res_rld_mat_df_pivtLonger, by="gene.id")


sig_genes_ids_merged_with_rlog <- sigLRT_genes %>% dplyr::select(c(gene.symbol, gene.id)) %>% left_join(res_rld_mat_df, by="gene.id")
rlog_4Clustering <- sig_genes_ids_merged_with_rlog %>% dplyr::select(!c(gene.id))
rownames(rlog_4Clustering) <- make.names(rlog_4Clustering[,1], unique = TRUE)
rlog_4Clustering <- rlog_4Clustering %>% dplyr::select(!"gene.symbol")
rlog_4Clustering_matrix <- data.matrix(rlog_4Clustering)

annot_coldf <- as.data.frame(colData(dseqObject)[,"condition", drop=FALSE])

pdf(file = paste0("figures/",contrast_input,"/Clustering_sig_DE",contrast_input,"_padjust",padj_val,".pdf"), width=12, height=9)
pheatmap(rlog_4Clustering_matrix, main = paste0("Clustering of Sig DGE ", contrast_input), scale = "row",annotation_names_row=FALSE, annotation_col=annot_coldf)
dev.off()
outresults = list(DeSeq2Results_df_condOnly_annot=results_per_contrast, results_per_contrast=results_per_contrast)
return(outresults)
}