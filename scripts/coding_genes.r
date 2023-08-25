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
metadata_df <- read.delim('data/sample_metadata.tsv', header = TRUE)
#metadata_df; needs 3 columns; samples-sample names matching rna-seq sample list, condition - treatment/ctrl, new_samples_name - new sample names
#rename  samples = 

#Assign new sample names
cond_assign_new_sample_names = ncol(metadata_df) > 2 & any(str_detect(names(metadata_df), "new_samples_name"))
if(cond_assign_new_sample_names){
idx_samples <- match(metadata_df$samples, colnames(cts_coding)) #check if samples match
names(cts_coding) <- metadata_df[,"new_samples_name"][idx_samples]

metadata_rown_df <- metadata_df %>% column_to_rownames("new_samples_name")

}

metadata_rown_df <- metadata_df %>% column_to_rownames("samples")


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
res <- results(dds, alpha=0.05) #use alpha of 0.05
summary(res)
results_per_contrast <- results(dds, alpha=0.05, contrast="condition_PARPi_vs_CTRL") #use alpha of 0.05
contrast_input="condition_PARPi_vs_CTRL"

#####################################################################################################################
############################### Step 3; Run deseq2 analysis for each contrast
#####################################################################################################################
DEResults_ls <- list()
run_multiple_contrasts <- function(contrast_input, dseqObject = dds, shrink = FALSE, padj = 0.05){

#contrast_label <- paste(contrast_input, collapse = "_") #used for plot labels
str_vec <- gsub("_vs_", "_", contrast_input); 
contrast_as_vect <- unlist(str_split(str_vec, pattern="_"))

#run DESeq analysis - normalization and filtering

results_per_contrast <- results(dseqObject, alpha=padj, contrast=contrast_as_vect) #use alpha of 0.05
#messsage(resultsNames(results_per_contrast))
# Shrink the log2 fold changes to be more accurate
if(shrink == TRUE){
results_per_contrast <- lfcShrink(dds, 
     contrast=contrast_as_vect, 
     type = "apeglm")	 
     # The coef will be dependent on what your contras was. and should be identical to what
}

dir.create(paste0("data/", contrast_input, "/"), recursive = TRUE, showWarnings = TRUE)
dir.create(paste0("figures/", contrast_input, "/"), recursive = TRUE, showWarnings = TRUE)

saveRDS(results_per_contrast, file=paste0("data/",contrast_input,"/Dseq2ResultsObject_",contrast_input,"_padjust.rds"))

DeSeq2Results_df <- as.data.frame(results_per_contrast) %>% rownames_to_column("gene.id")
# head(DeSeq2Results_df)
#DeSeq2Results_df <- DeSeq2Results_df %>% rownames_to_column(var = "ensembl.gene.id")
DeSeq2Results_df_annot <- left_join(DeSeq2Results_df, annot %>% dplyr::select(c(gene.id, gene.symbol, description,  entrez.gene.id)), by="gene.id") #merge gene annotations and deseq results
# head(DeSeq2Results_df_annot)
# #sum(is.na(DeSeq2Results_df_annot$gene.symbol))
# head(DeSeq2Results_df_annot)

#### do volcano plots
plt_volc_pval <- EnhancedVolcano(DeSeq2Results_df_annot,
    lab = DeSeq2Results_df_annot$gene.symbol,
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoff = padj,
    title = paste0(contrast_input," (p < ", padj, ")"),
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
    pCutoff = padj,
    title = paste0(contrast_input," (BH <", padj, ")"),
    widthConnectors = 0.75,
    labSize = 4.0,
    drawConnectors = TRUE,
    ylab = bquote(~Log[10]~ 'p adjusted'),
    pointSize = 1, 
    #labSize = 5,
    colAlpha = 0.2) + theme(plot.subtitle = element_blank())
    
ggsave(plt_volc_padjust, file = paste0("figures/",contrast_input,"/volcano_",contrast_input,"_padjust.pdf"))



#### Identifying significant genes
# Subset to return genes with padj < 0.05
sigLRT_genes <- DeSeq2Results_df_annot %>% 
  dplyr::filter(padj < padj)

# Get number of significant genes
nrow(sigLRT_genes)

# Compare to numbers we had from Wald test
# nrow(sigOE)
# nrow(sigKD)
return(results_per_contrast) #return results for specific contrasts
}


#run for all contrasts
DEResults_ls <- lapply(contrasts_ls[1], function(x) run_multiple_contrasts(contrast_input=x, dseqObject = dds))
names_contrasts_ls <- unlist(lapply(contrasts_ls, function(x) paste0(x, collapse = "_")))
names(DEResults_ls) <- names_contrasts_ls












################################################################################################
############################### Step #; GSEA
##############################################################################################
gsea_analysis <- function(name_de_res, specie_type = "Homo sapiens", category_tag = "H", nCategory_2show = 10, nShowBar=100,fractEnriched=0.5,
                            ggwidth = 9, ggheight = 7){
  require(DOSE)
  library(msigdbr)

# name_de_res="condition_PARP_CTRL"
dir.create(paste0("figures/", name_de_res,"/enrichment_Plots/"), recursive = TRUE, showWarnings = TRUE)



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

#do gene set enrichment using GO
gse <- gseGO(geneList=foldchanges, 
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
ggsave(gse_GO_plot2, file = paste0("figures/",name_de_res,"/enrichment_Plots/gseaPlot_gse_GO","_",gse$ID[x],"_" ,name_de_res,".pdf"))

}
purrr::map(1:length(gse$ID), func_gse_GO_plot2)


message("extracting GSE results table from gseGO object")
gse_results <- gse@result
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

ggsave(lollipop_gse, file = paste0("figures/lollipop_gse","_",name_de_res,"_top_frac_",nCategory_2show,".pdf"), width = ggwidth, height = ggheight)



message("performing GSEA using MSigDb ", category_tag," Specie - ", specie_type)
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
ggsave(gsea_plot2, file = paste0("figures/",name_de_res,"/enrichment_Plots/gseaPlot_gsea_msigdbr_cat_",category_tag,"_",msig_GSEA_obj$ID[x],"_" ,name_de_res,".pdf"))

}
purrr::map(1:length(msig_GSEA_obj$ID), get_gsea_plot2)

message("GSEA ridge plots")
#ridge plot; enrichment distribution
enrich_dist_plot_out <- ridgeplot(msig_GSEA_obj) + labs(x = "enrichment distribution")
ggsave(enrich_dist_plot_out, file = paste0("figures/",name_de_res,"/ridgePlot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_",".pdf"), width=18, height=27)

message("extracting GSEA results table from GSEA object")
msig_GSEA_results <- msig_GSEA_obj@result
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


# #cnet plots
# ## convert gene ID to Symbol
# edox <- setReadable(msig_GSEA_obj, 'org.Hs.eg.db', 'ENTREZID')
# p1 <- cnetplot(edox, foldChange=foldchanges)
# ## categorySize can be scaled by 'pvalue' or 'geneNum'
# p2 <- cnetplot(edox, categorySize="pvalue", foldChange=foldchanges)
# # p3 <- cnetplot(edox, foldChange=foldchanges, circular = TRUE, colorEdge = TRUE) 
# #cnetPlot_cowp <- cowplot::plot_grid(p1, p2, ncol=2, labels=LETTERS[1:2])
# ggsave(p1, file = paste0("figures/cnetplot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,".pdf"), width = ggwidth, height = ggheight)


#make list of results
gseaObj_name <- paste0("msig_category_",category_tag) #what name should be given to the gsea object in list
out_gsea_fun <- list(rankedFC=foldchanges,gseGO=gse, gseaObj_name = msig_GSEA_obj)
saveRDS(out_gsea_fun, file=paste0("data/GSE_GO_plus_msigdbr_cat_",category_tag,"_",name_de_res,".rds"))

# return(out_gsea_fun)
return(out_gsea_fun)
}

#paste0("figures/gsea_msigdbr_H_type_CKI_DMSO_EnrichDistrib.pdf")
#run gsea for all contrasts
# lapply(names_contrasts_ls, function(x) gsea_analysis(x))
h_gsea <- lapply(names_contrasts_ls, function(x) gsea_analysis(x, category_tag = "H", nCategory_2show = 100))
c2_gsea <- lapply(names_contrasts_ls, function(x) gsea_analysis(x, category_tag = "C2", ggwidth = 15, ggheight = 12))
names(h_gsea)
names(h_gsea[[1]][["rankedFC"]])

# c2_gsea[[1]]
# head(h_gsea[[1]]@result)
# slotnames(h_gsea)
# slotNames(h_gsea)

edox <- setReadable(c2_gsea[[1]], 'org.Hs.eg.db', 'ENTREZID')

tst <- h_gsea[[1]][["rankedFC"]]
#do gene set enrichment using GO
gse <- gseGO(geneList=tst, 
             ont ="ALL", 
             keyType = "ENTREZID", 
              eps = 0,
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "none")
             
            #  org.Hs.eg.db