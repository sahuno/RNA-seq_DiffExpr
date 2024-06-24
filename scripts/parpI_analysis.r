#name: samuel ahuno
#load packages
library(DESeq2) #main desdeq analysis
library(tidyverse)
library(EnhancedVolcano) # used for making volcano plots
library(viridis)
library(magrittr)
library("pheatmap")
library(data.table)
library("ggrepel")
library(ggfortify)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
library(cowplot)
library(ggfortify)
library(clusterProfiler)

########################################################################################
# set prroject  parameters
########################################################################################
#set wehere count data are;
source_dir = "/work/greenbaum/projects/RNA_seq_pipeline/projects/Project_13727_B/" 
#proj_name = "2023_BRCA_PARP_DESeq2" #project name

message("workig dir is - ", getwd())
#create folders in not aleady there
dir.create(paste0("figures/"), recursive = TRUE, showWarnings = TRUE)
dir.create(paste0("data/"), recursive = TRUE, showWarnings = TRUE)


blind_transform <- FALSE #should the rlog transformation be blind
drop_samples <- c("veh_2") #samples to drop from analysis
########################################################################################
## #read in data, extract coding genes counts 
########################################################################################
#load datasets from rna-seq results folder - sasha
annot <- fread(paste0(source_dir,'annot.tsv')) 
counts_annot <- fread(paste0(source_dir,'CT/counts_annot.tsv'))
cts <- fread(paste0(source_dir,'CT/counts.tsv'))

#get protein coding genes only
gs.coding_ids <- annot[gene.type=="protein_coding",.(gene.id)]
cts_coding_df <- cts[gs.coding_ids, on = "gene.id"] %>% as.data.frame()
cts_coding <- cts_coding_df %>% column_to_rownames(var = "gene.id") # #convert colnames to rowids

#### repeat landscape
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
metadata_df <- read.delim('/work/greenbaum/users/ahunos/apps/RNA-seq_DiffExpr/data/sample_metadata.tsv', header = TRUE)
#metadata_df; needs 3 columns; samples-sample names matching rna-seq sample list, condition - treatment/ctrl, new_samples_name - new sample names
#rename  samples = 

#Assign new sample names
cond_assign_new_sample_names = ncol(metadata_df) > 2 & any(str_detect(names(metadata_df), "new_samples_name"))

if(cond_assign_new_sample_names){

  #drop samples
if(!is.null(drop_samples)){
metadata_df <- metadata_df %>% dplyr::filter(!samples %in% drop_samples)
cts_coding <- cts_coding %>% dplyr::select(!all_of(drop_samples))
}

#assign new sample names
idx_samples <- match(metadata_df$samples, colnames(cts_coding)) #check if samples match
names(cts_coding) <- metadata_df[,"new_samples_name"][idx_samples]

metadata_rown_df <- metadata_df %>% column_to_rownames("new_samples_name")

}else{
#drop samples
# if(!is.null(drop_samples)){
# metadata_df <- metadata_df %>% dplyr::filter(!samples %in% drop_samples)
# }

metadata_rown_df <- metadata_df %>% column_to_rownames("samples")
}



##exclude samples from dseq analysis
# cts_coding <- cts_coding %>% dplyr::select(!veh_2) #remove outlier samples `veh_2`
##########################
#function to make dseq object
make_dseq_obj <- function(ref_col = NULL){
if(any(str_detect(names(metadata_rown_df), "condition"))){
    message("using condition column to assign samples to groups")

metadata_rown_df$condition <- factor(metadata_rown_df$condition)

dds <- DESeqDataSetFromMatrix(countData = cts_coding,
                              colData = metadata_rown_df,
                              design = ~ condition)
                              
#reset levels
message("setting CTRL as reference group")
dds$condition <- relevel(dds$condition, ref = ref_col)

#create contrast and run for all for all contrasts
#ref_col = "CTRL"
cond_levels <- levels(dds$condition) #get condition levels
#cond_levels <- c(cond_levels, "treat_2")
cond_case_level <- cond_levels[!grepl(ref_col, cond_levels)] #take non-base levels in conditions
contrasts_ls <- paste0("condition_", cond_case_level, "_vs_", ref_col)


out_dseq <- list(constrasts = contrasts_ls, dSeqObj = dds)
return(out_dseq)
}

}

#create dseq object and extract contrasts
dseq_func_out <- make_dseq_obj(ref_col = "CTRL")
dds <- dseq_func_out[["dSeqObj"]]
contrasts_ls <- dseq_func_out[["constrasts"]]

################################################################################################
##pre - filtering
keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
dds <- dds[keep,] #filter dds
# head(counts(dds))



################################################################################################
############################### Step 2; Run deseq2 analysis
##############################################################################################
#run DESeq analysis - normalization and filtering
dds <- DESeq(dds)
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

#save normalizeds counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)


# Transform counts for data visualization
rld <- rlog(dds, 
	    blind=blind_transform)

# Plot PCA 
pltPCA_samples <- plotPCA(rld, 
	intgroup=c("condition", "samples"))
ggsave(pltPCA_samples, file = paste0("figures/PCA_by_condition_and_samples.pdf"), width = 12, height = 12)


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

#save a copy of the rlog matrix
fwrite(rld_mat %>% data.frame(), file=paste0("data/DESeq2_rlog_Transform_Blind",blind_transform,".tsv"), sep="\t")



# Plot heatmap
pltHeatmap_samples <- pheatmap(rld_cor, 
	 annotation = metadata_rown_df[,c("condition"), drop = FALSE])
ggsave(pltHeatmap_samples, file = paste0("figures/heatmap_conditions.pdf"), width = 12, height = 12)

res_rld_mat_df <- rld_mat %>%
  data.frame() %>%
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
# contrast_input <- "PARPi_vs_CTRL" #uncomment to test run function

# metadata(dds_copy)
# elementMetadata(dds_copy)
# dseqObject <- dds_copy[,-idx]

#####################################################################################################################
############################### Step 3; Run deseq2 analysis for each contrast
#####################################################################################################################
DEResults_ls <- list()
run_multiple_contrasts <- function(contrast_input, dseqObject = dds, shrink = FALSE, padj_val = 0.05, export_transformed_counts = FALSE){

dir.create(paste0("data/", contrast_input, "/"), recursive = TRUE, showWarnings = TRUE)
dir.create(paste0("figures/", contrast_input, "/"), recursive = TRUE, showWarnings = TRUE)

#contrast_label <- paste(contrast_input, collapse = "_") #used for plot labels
str_vec <- gsub("_vs_", "_", contrast_input); 
contrast_as_vect <- unlist(str_split(str_vec, pattern="_"))

# dseqObject <- dseqObject[,-idx]


# export Transform counts for data visualization
if(export_transformed_counts){
dseqObject_copy <- dseqObject
colDat <- colData(dseqObject_copy) %>% as.data.frame()
idx_samples <- colDat %>% dplyr::filter(condition %in% contrast_as_vect) %>% rownames_to_column("samples.id") %>% pull("samples.id")
dds_copy_subset <- dseqObject_copy[,idx_samples]
dds_copy_subset <- DESeq(dds_copy_subset)


rld <- rlog(dds_copy_subset, #export data
	    blind=TRUE)

rld_mat <- assay(rld)
rld_df <- rld_mat %>% as.data.frame() %>% rownames_to_column("gene.id")
# head(rld_df)
#save a copy of the rlog matrix
write_tsv(rld_df, file=paste0("data/",contrast_input,"/DESeq2_rlog_Transform_BlindTRUE_",contrast_input,".tsv"))
}

# dt_t <- fread("/juno/work/greenbaum/users/ahunos/apps/RNA-seq_DiffExpr/data/condition_PARPi_vs_CTRL/DESeq2_rlog_Transform_BlindTRUE.tsv")
# head(dt_t)
#run DESeq analysis - normalization and filtering

results_per_contrast <- results(dseqObject, alpha=padj_val, contrast=contrast_as_vect) #use alpha of 0.05
#messsage(resultsNames(results_per_contrast))
# Shrink the log2 fold changes to be more accurate
if(shrink == TRUE){
res <- lfcShrink(dds, 
    #  contrast=contrast_as_vect, 
     coef = contrast_as_vect, #apeglm only works with `coef` argument`
     type = "apeglm"
     )	 
     # The coef will be dependent on what your contras was. and should be identical to what
}



# pdf("MA_plot_lfcShrink_dmso_vrs_all_treated.pdf")
# plotMA(res, ylim=c(-2,2))
# dev.off()




saveRDS(results_per_contrast, file=paste0("data/",contrast_input,"/Dseq2ResultsObject_",contrast_input,"_padjust.rds"))

DeSeq2Results_df <- as.data.frame(results_per_contrast) %>% rownames_to_column("gene.id")
DeSeq2Results_df_annot <- left_join(DeSeq2Results_df, annot %>% dplyr::select(c(gene.id, gene.symbol, description,  entrez.gene.id)), by="gene.id") #merge gene annotations and deseq results
write.table(DeSeq2Results_df_annot, file = paste0("data/",contrast_input,"/Dseq2Results_",contrast_input,".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE) #save deseq results with annotation as table

#### do volcano plots
plt_volc_pval <- EnhancedVolcano(DeSeq2Results_df_annot,
    lab = DeSeq2Results_df_annot$gene.symbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoff = padj_val,
    title = paste0(contrast_input," (p < ", padj_val, ")"),
    pointSize = 1, labSize = 5,
    colAlpha = 0.2) +
    theme(plot.subtitle = element_blank())
ggsave(plt_volc_pval, file = paste0("figures/",contrast_input,"/volcano_",contrast_input,"_pvals.pdf"))

#with adjusted p - values
DeSeq2Results_df_annot_padjustOnly <- DeSeq2Results_df_annot %>% dplyr::select(!pvalue) %>% mutate(pvalue = padj) #recreate data frame for ploting
plt_volc_padjust <- EnhancedVolcano(DeSeq2Results_df_annot_padjustOnly,
    lab = DeSeq2Results_df_annot_padjustOnly$gene.symbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoff = padj_val,
    title = paste0(contrast_input," (BH <", padj_val, ")"),
    widthConnectors = 0.75,
    labSize = 4.0,
    drawConnectors = TRUE,
    ylab = bquote('-' ~Log[10]~ 'p adjusted'),
    pointSize = 1, 
    #labSize = 5,
    colAlpha = 0.2) + theme(plot.subtitle = element_blank())
    
ggsave(plt_volc_padjust, file = paste0("figures/",contrast_input,"/volcano_",contrast_input,"_padjust.pdf"))

#### Identifying significant genes
# Subset to return genes with padj < 0.05
sigLRT_genes <- DeSeq2Results_df_annot %>% 
  dplyr::filter(padj < padj_val)

# Get number of significant genes
# nrow(sigLRT_genes)
# print(head(sigLRT_genes))
# sigLRT_genes$gene.symbol

sig_genes_merged_with_rlog <- sigLRT_genes %>% dplyr::select(c(gene.symbol, gene.id)) %>% left_join(res_rld_mat_df_pivtLonger, by="gene.id")


# # Make heatmap of significant genes
heatmap_DE_sig_plot <- ggplot(sig_genes_merged_with_rlog, 
                          aes(x=samples, y=gene.symbol, fill=rlog_transformed_counts)) + 
                            geom_raster() + scale_fill_viridis() + 
                            labs(title = paste0("Heatmap of Differentially Gene Expression ", contrast_input)) + 
                              theme(axis.text.x=element_text(angle=65, hjust=1), 
                                    legend.position = "top")

ggsave(heatmap_DE_sig_plot, file = paste0("figures/",contrast_input,"/Heatmap_sig_DE",contrast_input,"_padjust",padj_val,".pdf"))


sig_genes_ids_merged_with_rlog <- sigLRT_genes %>% dplyr::select(c(gene.symbol, gene.id)) %>% 
                                            left_join(res_rld_mat_df, by="gene.id") #%>% 
                                                                        #dplyr::select(gene.id)
#head(sig_genes_ids_merged_with_rlog)
rlog_4Clustering <- sig_genes_ids_merged_with_rlog %>% dplyr::select(!c(gene.id))
rownames(rlog_4Clustering) <- make.names(rlog_4Clustering[,1], unique = TRUE)
rlog_4Clustering <- rlog_4Clustering %>% dplyr::select(!"gene.symbol")
rlog_4Clustering_matrix <- data.matrix(rlog_4Clustering)

annot_coldf <- as.data.frame(colData(dds)[,"condition", drop=FALSE])


pdf(file = paste0("figures/",contrast_input,"/Clustering_sig_DE",contrast_input,"_padjust",padj_val,".pdf"), width=12, height=9)
pheatmap(rlog_4Clustering_matrix, main = paste0("Clustering of Sig DGE ", contrast_input), annotation_col=annot_coldf)
dev.off()


return(results_per_contrast) #return results for specific contrasts
}


################################################################################################
############################### Step #; Run DE for each contrast
##############################################################################################
#run for all contrasts
DEResults_ls <- lapply(contrasts_ls, function(x) run_multiple_contrasts(contrast_input=x, dseqObject = dds, shrink = FALSE, export_transformed_counts = TRUE))
names(DEResults_ls) <- contrasts_ls
#summary(DEResults_ls[[1]])
# data.frame(DEResults_ls[[1]]) %>% dplyr::filter(padj < 0.05) %>% dim()










################################################################################################
############################### Step #; GSEA
##############################################################################################
gsea_analysis <- function(name_de_res, specie_type = "Homo sapiens", category_tag = "H", nCategory_2show = 10, nShowBar=100,
                          fractEnriched=0.5, do_gse_GO = TRUE, do_GSEA=TRUE, ggwidth = 9, ggheight = 7){
  require(DOSE)
  library(msigdbr)

# name_de_res="condition_PARP_CTRL"
dir.create(paste0("figures/", name_de_res,"/enrichment_Plots/GO/"), recursive = TRUE, showWarnings = TRUE)



de_res_df <- DEResults_ls[[name_de_res]] %>% as.data.frame() %>% rownames_to_column(var = "gene.id")

  DEResults_type_trt_ctrl_annt  <- left_join(de_res_df, annot %>% dplyr::select(c(ensembl.gene.id, gene.id, gene.symbol, description, entrez.gene.id)), by="gene.id") #merge gene annotations and deseq results
#gene set enrichment analysis
## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(DEResults_type_trt_ctrl_annt, entrez.gene.id != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez.gene.id) == F), ]

## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrez.gene.id

## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

if(do_gse_GO){
    #do gene set enrichment using GO
gse <- gseGO(geneList = foldchanges, 
             ont ="ALL", 
             keyType = "ENTREZID", 
              eps = 0,
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "none")
dotplt_gse_GO_out <- dotplot(gse, showCategory=nCategory_2show, split=".sign") + facet_grid(~factor(.sign, levels=c('suppressed', 'activated'))) +
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 28, face="bold"))
    #theme(plot.title = element_text(size=22)) 
ggsave(dotplt_gse_GO_out, file = paste0("figures/",name_de_res,"/dotplot_gse_GO","_",name_de_res,"_top_",nCategory_2show,".pdf"), width = ggwidth, height = ggheight)


## Output images for all significant KEGG pathways
func_gse_GO_plot2  <- function(x) {
gse_GO_plot2 <-  enrichplot::gseaplot2(gse, geneSetID=gse$ID[x], title = gse$ID[x])
ggsave(gse_GO_plot2, file = paste0("figures/",name_de_res,"/enrichment_Plots/GO/gseaPlot_gse_GO","_",gse$ID[x],"_" ,name_de_res,".pdf"))

}
purrr::map(1:length(gse$ID), func_gse_GO_plot2)


message("extracting GSE results table from gseGO object")
gse_results <- gse@result

write.table(gse_results, 
file=paste0("data/",name_de_res,"/GSE_GO_enrichedSets.tsv"),
sep="\t", row.names = FALSE, col.names = TRUE)

#rank gene sets by normalized  enrichment scores and q-values
gse_results_sig_qval <- gse_results %>% dplyr::filter(qvalues < 0.05) %>% arrange(desc(abs(NES))) #%>% head() #sigGSEA

##do lollipop plots of enrichment
lollipop_gse <- gse_results_sig_qval %>% arrange(desc(abs(NES))) %>% 
                top_frac(fractEnriched) %>% 
                ggplot(aes(x=ID, y = NES)) +  
geom_segment(aes(x=ID, xend=ID, y=0, yend=NES)) + geom_point()  + 
coord_flip() +
  labs(title = paste0("normalized enrichment scores (NES) for sig GO ", " gene sets (q < 0.05)"),
  x= "ID", y="normalized enrichment scores (NES)") +
   geom_vline(xintercept = 0, color = "blue") +
    theme(plot.title = element_text(hjust = -1))

ggsave(lollipop_gse, file = paste0("figures/",name_de_res,"/lollipop_GSE_Go_top_",fractEnriched,".pdf"), width = ggwidth, height = ggheight)

}


if(do_GSEA){
    message("performing GSEA using MSigDb ", category_tag," Specie - ", specie_type)
    dir.create(paste0("figures/", name_de_res,"/enrichment_Plots/GSEA_MSigDb_",category_tag,"/"), recursive = TRUE, showWarnings = TRUE)

#get data base to use
m_t2g <- msigdbr(species = specie_type, category = category_tag) %>% 
  dplyr::select(gs_name, entrez_gene) #hall mark genes

# Run GSEA analysis
msig_GSEA_obj <- GSEA(foldchanges, TERM2GENE = m_t2g, verbose = TRUE, eps = 0)

#, category = "C2"
#make dotplot
dotplt_gsea_out <- dotplot(msig_GSEA_obj, showCategory=nCategory_2show, split=".sign") + facet_grid(~factor(.sign, levels=c('suppressed', 'activated'))) +
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 28, face="bold"))
ggsave(dotplt_gsea_out, file = paste0("figures/",name_de_res,"/dotplot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_",nCategory_2show,".pdf"), width = ggwidth, height = ggheight)


message("GSEA enrichment plots")
## Output images for all significant KEGG pathways
get_gsea_plot2  <- function(x) {
#    pathview(gene.data = foldchanges, pathway.id = gseaKEGG_results$ID[x], species = "hsa", 
#        limit = list(gene = 2, cpd = 1))
gsea_plot2 <-  enrichplot::gseaplot2(msig_GSEA_obj, geneSetID=msig_GSEA_obj$ID[x], title = msig_GSEA_obj$ID[x])
ggsave(gsea_plot2, file = paste0("figures/",name_de_res,"/enrichment_Plots/GSEA_MSigDb_",category_tag,"/","gseaPlot_gsea_msigdbr_cat_",category_tag,"_",msig_GSEA_obj$ID[x],"_" ,name_de_res,".pdf"))

}
purrr::map(1:length(msig_GSEA_obj$ID), get_gsea_plot2)

message("GSEA ridge plots")
#ridge plot; enrichment distribution
enrich_dist_plot_out <- ridgeplot(msig_GSEA_obj) + labs(x = "enrichment distribution")
ggsave(enrich_dist_plot_out, file = paste0("figures/",name_de_res,"/ridgePlot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_",".pdf"), width=18, height=27)

message("extracting GSEA results table from GSEA object")
msig_GSEA_results <- msig_GSEA_obj@result
write.table(msig_GSEA_results, 
file=paste0("data/",name_de_res,"/GSEA_msigdbr_cat_",category_tag,"_enrichedSets.tsv"),
sep="\t", row.names = FALSE, col.names = TRUE)

#rank gene sets by normalized  enrichment scores and q-values
msig_GSEA_results_sig_qval <- msig_GSEA_results %>% dplyr::filter(qvalues < 0.05) %>% arrange(desc(abs(NES))) #%>% head() #sigGSEA
print(head(msig_GSEA_results_sig_qval)) ##test results


message("performing lollipop plots of sig enriched gene sets")
#do lollipop plots of GSEA enrichment
lollipop_gsea <- msig_GSEA_results_sig_qval %>% arrange(desc(abs(NES))) %>% top_frac(fractEnriched) %>% ggplot(aes(x=ID, y = NES)) +  
geom_segment(aes(x=ID, xend=ID, y=0, yend=NES)) + geom_point()  + 
coord_flip() +
  labs(title = paste0("normalized enrichment scores (NES) for sig ",category_tag, " gene sets (q < 0.05)"),
  x= "ID", y="normalized enrichment scores (NES)") +
   geom_vline(xintercept = 0, color = "blue") +
    theme(plot.title = element_text(hjust = -1))

ggsave(lollipop_gsea, file = paste0("figures/",name_de_res,"/lollipop_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_frac_",nCategory_2show,".pdf"), width = ggwidth, height = ggheight)

}
}
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
h_gsea <- lapply(contrasts_ls, function(x) gsea_analysis(x, category_tag = "H", nCategory_2show = 100, ggwidth = 15, ggheight = 12))
c2_gsea <- lapply(contrasts_ls, function(x) gsea_analysis(x, do_gse_GO=FALSE,category_tag = "C2", ggwidth = 15, ggheight = 12))
# names(h_gsea)
# names(h_gsea[[1]][["rankedFC"]])

# c2_gsea[[1]]
# head(h_gsea[[1]]@result)
# slotnames(h_gsea)
# slotNames(h_gsea)
