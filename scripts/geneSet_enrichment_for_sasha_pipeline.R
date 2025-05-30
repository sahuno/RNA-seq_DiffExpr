library(tidyverse)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
require(DOSE)
library(msigdbr)

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
  make_option(c("-d", "--de_file"), type="character",
              default="/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/DE/DE.xlsx",
              help="Path to DE Excel file"),
  make_option(c("-N", "--name_de_res"), type="character", default="triplicates_epigenetics_diyva",
              help="Name identifier for DE results")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Assign variables from parsed options
specie_type     <- opt$specie_type
category_tag    <- opt$category_tag
fractEnriched   <- opt$fractEnriched
do_gse_GO       <- opt$do_gse_GO
nCategory_2show <- opt$nCategory_2show
nShowBar        <- opt$nShowBar
ggwidth         <- opt$ggwidth
ggheight        <- opt$ggheight
de_AzaMouse     <- opt$de_file
name_de_res     <- opt$name_de_res

# Bold console messages
if (!requireNamespace("crayon", quietly = TRUE)) {
  install.packages("crayon")
}
library(crayon)

de <- read.xlsx(de_AzaMouse, sheet=2)
de_sigChanges <- rownames(de)[de$padj < 0.05 & !is.na(de$padj)]
res_entrez <- dplyr::filter(de, entrez.gene.id != "NA")
## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez.gene.id) == F), ]

## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrez.gene.id

## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

message(crayon::bold(paste0("performing GSEA using MSigDb ", category_tag, " Specie - ", specie_type)))
dir.create(paste0("figures/", name_de_res,"/enrichment_Plots/GSEA_MSigDb_",category_tag,"/"), recursive = TRUE, showWarnings = TRUE)
dir.create(paste0("data/", name_de_res,"/enrichment_Plots/GSEA_MSigDb_",category_tag,"/"), recursive = TRUE, showWarnings = TRUE)


#get data base to use
m_t2g <- msigdbr(species = specie_type, category = category_tag) %>% 
  dplyr::select(gs_name, entrez_gene) #hall mark genes

# Run GSEA analysis
msig_GSEA_obj <- GSEA(foldchanges, TERM2GENE = m_t2g, verbose = TRUE, eps = 0)

#make dotplot
dotplt_gsea_out <- dotplot(msig_GSEA_obj, showCategory=nCategory_2show, split=".sign") + facet_grid(~factor(.sign, levels=c('suppressed', 'activated'))) +
  theme(strip.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 28, face="bold"))

ggsave(dotplt_gsea_out, file = paste0("figures/",name_de_res,"/dotplot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_top_",nCategory_2show,".pdf"), width = ggwidth, height = ggheight)

message("GSEA enrichment plots")
## Output images for all significant KEGG pathways
get_gsea_plot2  <- function(x) {
gsea_plot2 <-  enrichplot::gseaplot2(msig_GSEA_obj, geneSetID=msig_GSEA_obj$ID[x], title = msig_GSEA_obj$ID[x])
ggsave(gsea_plot2, file = paste0("figures/",name_de_res,"/enrichment_Plots/GSEA_MSigDb_",category_tag,"/","gseaPlot_gsea_msigdbr_cat_",category_tag,"_",msig_GSEA_obj$ID[x],"_" ,name_de_res,".pdf"))
}
purrr::map(1:length(msig_GSEA_obj$ID), get_gsea_plot2)

message("GSEA ridge plots")
enrich_dist_plot_out <- ridgeplot(msig_GSEA_obj) + labs(x = "enrichment distribution")
ggsave(enrich_dist_plot_out, file = paste0("figures/",name_de_res,"/ridgePlot_gsea_msigdbr_cat_",category_tag,"_",name_de_res,"_",".pdf"), width=18, height=27)

message("extracting GSEA results table from GSEA object")
msig_GSEA_results <- msig_GSEA_obj@result
write.table(msig_GSEA_results, 
file=paste0("data/",name_de_res,"/GSEA_msigdbr_cat_",category_tag,"_enrichedSets.tsv"),
sep="\t", row.names = FALSE, col.names = TRUE)

#rank gene sets by normalized  enrichment scores and q-values
msig_GSEA_results_sig_qval <- msig_GSEA_results %>% dplyr::filter(qvalue < 0.05) %>% arrange(desc(abs(NES))) #%>% head() #sigGSEA
# print(head(msig_GSEA_results_sig_qval)) ##test results

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


#######################
# GSEA GO

# if(do_gse_GO){
# if(specie_type == "Homo sapiens"){
#   gse <- gseGO(geneList = foldchanges, 
#              ont ="ALL", 
#              keyType = "ENTREZID", 
#               eps = 0,
#              minGSSize = 3, 
#              maxGSSize = 800, 
#              pvalueCutoff = 0.05, 
#              verbose = TRUE, 
#              OrgDb = "org.Hs.eg.db", 
#              pAdjustMethod = "none")
# }else if(specie_type == "Mus musculus"){
#   gse <- gseGO(geneList = foldchanges, 
#              ont ="ALL", 
#              keyType = "ENTREZID", 
#               eps = 0,
#              minGSSize = 3, 
#              maxGSSize = 800, 
#              pvalueCutoff = 0.05, 
#              verbose = TRUE, 
#              OrgDb = "org.Mm.eg.db", 
#              pAdjustMethod = "none")
# }
# dotplt_gse_GO_out <- dotplot(gse, showCategory=nCategory_2show, split=".sign") + facet_grid(~factor(.sign, levels=c('suppressed', 'activated'))) +
#   theme(strip.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 28, face="bold"))
#     #theme(plot.title = element_text(size=22)) 
# ggsave(dotplt_gse_GO_out, file = paste0("figures/",name_de_res,"/dotplot_gse_GO","_",name_de_res,"_top_",nCategory_2show,".pdf"), width = ggwidth, height = ggheight)
# }

# ## Output images for all significant KEGG pathways
# func_gse_GO_plot2  <- function(x) {
# gse_GO_plot2 <-  enrichplot::gseaplot2(gse, geneSetID=gse$ID[x], title = gse$ID[x])
# ggsave(gse_GO_plot2, file = paste0("figures/",name_de_res,"/enrichment_Plots/GO/gseaPlot_gse_GO","_",gse$ID[x],"_" ,name_de_res,".pdf"))
# }
# purrr::map(1:length(gse$ID), func_gse_GO_plot2)
