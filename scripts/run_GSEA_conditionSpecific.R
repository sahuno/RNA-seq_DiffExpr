################################################################################################
############################### Step #; GSEA
##############################################################################################
# name_de_res = "condition_AZCT_vs_DMSO"
gsea_analysis <- function(name_de_res, results_df_with_annot=DeSeq2Results_df_annot, 
                          specie_type = "Mus musculus", 
                          category_tag = "M2", subcategory_tag = "GCP", nCategory_2show = 10, nShowBar=100,
                          qval_filter = 0.05,
                          fractEnriched=0.5, do_gse_GO = TRUE, do_GSEA=TRUE, ggwidth = 9, ggheight = 7, ggwidthC2 = 7, ggheightC2 = 18){
  require(DOSE)
  library(msigdbr)

# Input validation
  # Check for valid species
  valid_species <- msigdbr_species()$species_name
  if(!specie_type %in% valid_species) {
    stop(paste0("Invalid species '", specie_type, "'. Valid options are: ", 
                paste(valid_species, collapse = ", ")))
  }
  
  # Check for valid collections
  if(do_GSEA) {
    available_collections <- msigdbr_collections()
    valid_collections <- unique(available_collections$gs_collection)
    
    if(!category_tag %in% valid_collections) {
      stop(paste0("Invalid collection '", category_tag, "'. Valid options are: ", 
                  paste(valid_collections, collapse = ", ")))
    }
    
    # Check for valid subcategory if provided
    if(!is.null(subcategory_tag)) {
      valid_subcategories <- available_collections %>% 
        filter(gs_collection == category_tag) %>% 
        pull(gs_subcollection) %>% 
        unique() %>% 
        na.omit()
      
      if(!subcategory_tag %in% valid_subcategories) {
        stop(paste0("Invalid subcategory '", subcategory_tag, 
                    "' for collection '", category_tag, 
                    "'. Valid options are: ", 
                    paste(valid_subcategories, collapse = ", ")))
      }
    }
    
    # Show what gene sets are being used
    collection_info <- available_collections %>% 
      filter(gs_collection == category_tag)
    
    if(!is.null(subcategory_tag)) {
      collection_info <- collection_info %>% 
        filter(gs_subcollection == subcategory_tag)
    }
    
    message(paste0("Using collection: ", category_tag, 
                   ifelse(!is.null(subcategory_tag), 
                          paste0(" (subcategory: ", subcategory_tag, ")"), ""),
                   " for species: ", specie_type))
    message(paste0("Number of gene sets available: ", 
                   sum(collection_info$num_genesets)))
  }
  
  # Create output directories
  # name_de_res="condition_PARP_CTRL"
  dir.create(paste0("figures/", name_de_res,"/enrichment_Plots/GO/"), 
             recursive = TRUE, showWarnings = FALSE)




# de_res_df <- DEResults_ls[[name_de_res]] %>% as.data.frame() %>% rownames_to_column(var = "gene.id")
# DEResults_type_trt_ctrl_annt  <- left_join(de_res_df, annot %>% dplyr::select(c(ensembl.gene.id, gene.id, gene.symbol, description, entrez.gene.id)), by="gene.id") #merge gene annotations and deseq results

DEResults_type_trt_ctrl_annt <- results_df_with_annot
DEResults_type_trt_ctrl_annt %>% dplyr::filter(is.na(padj))
DEResults_type_trt_ctrl_annt %>% dplyr::filter(is.na(pvalue))

#gene set enrichment analysis
## Remove any NA values (reduces the data by quite a bit)
res_entrez <- dplyr::filter(DEResults_type_trt_ctrl_annt, !is.na(entrez.gene.id) & entrez.gene.id != "NA")
message(paste0("Retained ", nrow(res_entrez), " genes with valid Entrez IDs out of ", 
               nrow(DEResults_type_trt_ctrl_annt), " total genes"))

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
             pvalueCutoff = qval_filter, 
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
             pvalueCutoff = qval_filter, 
             verbose = TRUE, 
             OrgDb = "org.Mm.eg.db", 
             pAdjustMethod = "none")
}
dotplt_gse_GO_out <- dotplot(gse, showCategory=nCategory_2show, split=".sign") + 
                      facet_grid(~factor(.sign, levels=c('suppressed', 'activated'))) +
                      labs(title = paste0(name_de_res, " GSE GO gene sets (q < 0.05)")) +
                      theme(strip.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 28, face="bold"))
    #theme(plot.title = element_text(size=22)) 
ggsave(dotplt_gse_GO_out, file = paste0("figures/",name_de_res,"/dotplot_gse_GO","_",name_de_res,"_top_",nCategory_2show,".pdf"), width = ggwidth, height = ggheight)

## Output images for all significant KEGG pathways
func_gse_GO_plot2  <- function(x) {
gse_GO_plot2 <-  enrichplot::gseaplot2(gse, geneSetID=gse$ID[x], title = gse$ID[x])
ggsave(gse_GO_plot2, file = paste0("figures/",name_de_res,"/enrichment_Plots/GO/gseaPlot_gse_GO","_",gse$ID[x],"_" ,name_de_res,".pdf"))
}
purrr::map(1:length(gse$ID), possibly(func_gse_GO_plot2, otherwise = NULL))


message("extracting GSE results table from gseGO object")
gse_results <- gse@result
write.table(gse_results, 
file=paste0("data/",name_de_res,"/GSE_GO_enrichedSets.tsv"),
sep="\t", row.names = FALSE, col.names = TRUE)

#rank gene sets by normalized  enrichment scores and q-values
gse_results_sig_qval <- gse_results %>% dplyr::filter(qvalue < qval_filter) %>% arrange(desc(abs(NES))) #%>% head() #sigGSEA

##do lollipop plots of enrichment
if(category_tag == "C2"){
  fractEnriched <- 0.01
  message("Using fractEnriched = 0.01 for C2 category")
  lollipop_gse <- gse_results_sig_qval %>% arrange(desc(abs(NES))) %>% 
                top_frac(fractEnriched) %>% 
                ggplot(aes(x=ID, y = NES)) +  
geom_segment(aes(x=ID, xend=ID, y=0, yend=NES)) + geom_point()  + 
coord_flip() +
  labs(title = paste0("normalized enrichment scores (NES) for sig GO ", " gene sets (q < ", qval_filter, ") ",name_de_res),
  x= "ID", y="normalized enrichment scores (NES)") +
   geom_vline(xintercept = 0, color = "blue") +
    theme(plot.title = element_text(hjust = -1))

ggsave(lollipop_gse, file = paste0("figures/",name_de_res,"/lollipop_GSE_Go_top_",fractEnriched,".pdf"), width = ggwidth, height = ggheight)
} else {
lollipop_gse <- gse_results_sig_qval %>% arrange(desc(abs(NES))) %>% 
                top_frac(fractEnriched) %>% 
                ggplot(aes(x=ID, y = NES)) +  
geom_segment(aes(x=ID, xend=ID, y=0, yend=NES)) + geom_point()  + 
coord_flip() +
  labs(title = paste0("normalized enrichment scores (NES) for sig GO ", " gene sets (q < ", qval_filter, ") ",name_de_res),
  x= "ID", y="normalized enrichment scores (NES)") +
   geom_vline(xintercept = 0, color = "blue") +
    theme(plot.title = element_text(hjust = -1))

ggsave(lollipop_gse, file = paste0("figures/",name_de_res,"/lollipop_GSE_Go_top_",fractEnriched,".pdf"), width = ggwidth, height = ggheight)
}
}


if(do_GSEA){

  message("performing GSEA using MSigDb collection=", category_tag, 
            ifelse(!is.null(subcategory_tag), 
                   paste0(", subcategory=", subcategory_tag), ""),
            ", Species: ", specie_type)
    
    dir.create(paste0("figures/", name_de_res,"/enrichment_Plots/GSEA_MSigDb_",
                      category_tag, 
                      ifelse(!is.null(subcategory_tag), 
                             paste0("_", subcategory_tag), ""),"/"), 
               recursive = TRUE, showWarnings = FALSE)

#get data base to use; subcategory
message("Retrieving MSigDb gene sets for category: ", category_tag, " and species: ", specie_type)
# m_t2g <- msigdbr::msigdbr(species = "Mus musculus", collection = "M2"); msigdbr::msigdbr(species = "mouse", collection = "M2", subcollection = "CGP")
# msigdbr::msigdbr(species = "Mus musculus", collection = "M2")
# msigdbr::msigdbr(species = "Mus musculus", collection = "M2", subcollection = "CGP")

# msigdbr::msigdbr(species = "Mus musculus", collection = "C2")


# m_t2g <- msigdbr::msigdbr(species = "Mus musculus", collection = "C2")
# m_t2g <- msigdbr(species = "Mus musculus", collection = "C2") %%
# msigdbr_collections(db_species = "Mm")

# Get database to use - UPDATED SYNTAX
    if(!is.null(subcategory_tag)) {
      m_t2g <- msigdbr(species = specie_type, 
                       collection = category_tag, 
                       subcollection = subcategory_tag) %>% 
        dplyr::select(gs_name, ncbi_gene)
    } else {
      m_t2g <- msigdbr(species = specie_type, 
                       collection = category_tag) %>% 
        dplyr::select(gs_name, ncbi_gene)
    }
    
    # Check if we got any gene sets
    n_genesets <- length(unique(m_t2g$gs_name))
    message(paste0("Retrieved ", n_genesets, " gene sets"))
    
    if(n_genesets == 0) {
      warning(paste0("No gene sets found for collection ", category_tag, 
                     " in species ", specie_type))
      return(NULL)
    }
# Run GSEA analysis
msig_GSEA_obj <- GSEA(foldchanges, TERM2GENE = m_t2g, verbose = TRUE, eps = 0)

if(length(msig_GSEA_obj$ID) == 0) {
    warning("No enriched gene sets found for ", name_de_res)
    return(NULL)
  }
#, category = "C2"
#make dotplot
message("creating GSEA dotplot")
dotplt_gsea_out <- dotplot(msig_GSEA_obj, showCategory=nCategory_2show, split=".sign") + 
                    facet_grid(~factor(.sign, levels=c('suppressed', 'activated'))) +
                    labs(title = paste0("GSEA MSigDb ", category_tag, " gene sets (q < ", qval_filter, ") for ", name_de_res) ) +
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
purrr::map(1:length(msig_GSEA_obj$ID), possibly(get_gsea_plot2, otherwise = NULL))

message("GSEA ridge plots")
#ridge plot; enrichment distribution
enrich_dist_plot_out <- ridgeplot(msig_GSEA_obj) + labs(x = "enrichment distribution")
enrich_dist_plot_out <- enrich_dist_plot_out + 
                          labs(title = paste0("GSEA MSigDb ", category_tag, " gene sets (q < ", qval_filter, ") for ", name_de_res)) +
  theme(plot.title = element_text(hjust = 0.5, size = 28, face="bold"))

ggsave(enrich_dist_plot_out, file = paste0("figures/",name_de_res,"/ridgePlot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_",".pdf"), width=18, height=27)

message("extracting GSEA results table from GSEA object")
msig_GSEA_results <- msig_GSEA_obj@result
write.table(msig_GSEA_results, 
file=paste0("data/",name_de_res,"/GSEA_msigdbr_cat_",category_tag,"_enrichedSets.tsv"),
sep="\t", row.names = FALSE, col.names = TRUE)

#rank gene sets by normalized  enrichment scores and q-values
msig_GSEA_results_sig_qval <- msig_GSEA_results %>% dplyr::filter(qvalue < qval_filter) %>% arrange(desc(abs(NES))) #%>% head() #sigGSEA
print(head(msig_GSEA_results_sig_qval)) ##test results


message("performing lollipop plots of sig enriched gene sets")
#do lollipop plots of GSEA enrichment
if(category_tag == "C2"){
  fractEnriched <- 0.04
  message(paste0("Using fractEnriched = ",fractEnriched," for C2 category"))

  topFrac_DF <- msig_GSEA_results_sig_qval %>% arrange(desc(abs(NES))) %>% top_frac(fractEnriched)
  lollipop_gsea <- topFrac_DF %>% ggplot(aes(x=ID, y = NES)) +  
geom_segment(aes(x=ID, xend=ID, y=0, yend=NES)) + geom_point()  + 
coord_flip() +
  labs(title = paste0("normalized enrichment scores (NES) for sig ",category_tag, " gene sets (q < ", qval_filter, ") ", name_de_res),
  x= "ID", y="normalized enrichment scores (NES)") +
   geom_vline(xintercept = 0, color = "blue") +
    theme(plot.title = element_text(hjust = -1))

ggsave(lollipop_gsea, file = paste0("figures/",name_de_res,"/lollipop_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_frac_",nCategory_2show,".pdf"), width = ggwidthC2, height = ggheightC2)
ggsave(lollipop_gsea, file = paste0("figures/",name_de_res,"/lollipop_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_frac_",nCategory_2show,".png"), width = ggwidthC2, height = ggheightC2)

} else {
  topFrac_DF <-  msig_GSEA_results_sig_qval %>% arrange(desc(abs(NES))) %>% top_frac(fractEnriched)  
lollipop_gsea <- topFrac_DF %>% ggplot(aes(x=ID, y = NES)) +  
geom_segment(aes(x=ID, xend=ID, y=0, yend=NES)) + geom_point()  + 
coord_flip() +
  labs(title = paste0("normalized enrichment scores (NES) for sig ",category_tag, " gene sets (q < ", qval_filter, ") ", name_de_res),
  x= "ID", y="normalized enrichment scores (NES)") +
   geom_vline(xintercept = 0, color = "blue") +
    theme(plot.title = element_text(hjust = -1))

ggsave(lollipop_gsea, file = paste0("figures/",name_de_res,"/lollipop_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_frac_",nCategory_2show,".pdf"), width = ggwidth, height = ggheight)
ggsave(lollipop_gsea, file = paste0("figures/",name_de_res,"/lollipop_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_frac_",nCategory_2show,".png"), width = ggwidth, height = ggheight)

}
}
}