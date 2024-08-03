################################################################################################
############################### Step #; GSEA
##############################################################################################
# name_de_res = "condition_AZCT_vs_DMSO"
gsea_analysis <- function(name_de_res, results_df_with_annot=DeSeq2Results_df_annot, specie_type = "Mus musculus", category_tag = "C2", nCategory_2show = 10, nShowBar=100,
                          fractEnriched=0.5, do_gse_GO = TRUE, do_GSEA=TRUE, ggwidth = 9, ggheight = 7){
  require(DOSE)
  library(msigdbr)

# name_de_res="condition_PARP_CTRL"
dir.create(paste0("figures/", name_de_res,"/enrichment_Plots/GO/"), recursive = TRUE, showWarnings = TRUE)



# de_res_df <- DEResults_ls[[name_de_res]] %>% as.data.frame() %>% rownames_to_column(var = "gene.id")
# DEResults_type_trt_ctrl_annt  <- left_join(de_res_df, annot %>% dplyr::select(c(ensembl.gene.id, gene.id, gene.symbol, description, entrez.gene.id)), by="gene.id") #merge gene annotations and deseq results

DEResults_type_trt_ctrl_annt <- results_df_with_annot
DEResults_type_trt_ctrl_annt %>% dplyr::filter(is.na(padj))
DEResults_type_trt_ctrl_annt %>% dplyr::filter(is.na(pvalue))

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
if(do_gse_GO){
if(specie_type == "Homo sapiens"){
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
}else if(specie_type == "Mus musculus"){
  gse <- gseGO(geneList = foldchanges, 
             ont ="ALL", 
             keyType = "ENTREZID", 
              eps = 0,
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Mm.eg.db", 
             pAdjustMethod = "none")
}
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
gse_results_sig_qval <- gse_results %>% dplyr::filter(qvalue < 0.05) %>% arrange(desc(abs(NES))) #%>% head() #sigGSEA

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
msig_GSEA_results_sig_qval <- msig_GSEA_results %>% dplyr::filter(qvalue < 0.05) %>% arrange(desc(abs(NES))) #%>% head() #sigGSEA
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