#cgas_sting pathway enrichment analyis
library(GSVA)
library(data.table)

### Get c gas genes
C_GAS_gene_list_KEGG<-fread("/juno/work/greenbaum/users/ahunos/apps/RNA-seq_DiffExpr/_GAS_gene_list_KEGG")
C_GAS_gene_list_KEGG<-C_GAS_gene_list_KEGG[-c(1),]
C_GAS_gene_list_KEGG<-data.frame(C_GAS_gene_list_KEGG)
colnames(C_GAS_gene_list_KEGG)<-"genes"

# find /path/to/search -user username
 
dseqObj <- readRDS("/work/greenbaum/users/ahunos/sandbox/diff_expr_misc/data/condition_AZCT_vs_DMSO/Dseq2ResultsObject_condition_AZCT_vs_DMSO_padjust.rds")
##Calculate ssGSES for each sample using normalized counts 
gsva_result_C_gas<-gsva(as.matrix(normalized_counts),list(C_GAS_gene_list_KEGG$genes),method="gsva",kcdf="Gaussian",verbose=FALSE )

gsva_result_C_gas_mat<-data.frame(t(gsva_result_C_gas))
gsva_result_C_gas_mat_Z_sore_transformed<-t(gsva_result_C_gas_mat %>% mutate_at(colnames(gsva_result_C_gas_mat), scale))
colnames(gsva_result_C_gas_mat_Z_sore_transformed)<-colnames(gsva_result_C_gas)
