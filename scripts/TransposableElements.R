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


#/juno/work/greenbaum/projects/TRI_EPIGENETIC/RNASeq_DE_TriEpi
########################################################################################
# set prroject  parameters
########################################################################################
#set wehere count data are;
source_dir = "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/" 
#proj_name = "2023_BRCA_PARP_DESeq2" #project name

setwd("/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/sandbox/")
message("workig dir is - ", getwd())
#create folders in not aleady there
dir.create(paste0("figures/"), recursive = TRUE, showWarnings = TRUE)
dir.create(paste0("data/"), recursive = TRUE, showWarnings = TRUE)
dir.create(paste0("figures/geneWiseNormalizedCounts"), recursive = TRUE, showWarnings = TRUE)


blind_transform <- TRUE #should the rlog transformation be blind
drop_samples <- NULL #samples to drop from analysis
ref_variable = "DMSO"
smallestGroupSize <- 3
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
metadata_df <- read.csv(file = paste0(source_dir,'metadata_triplicates.csv'),sep="," ,header = TRUE)
#metadata_df; needs 3 columns; samples-sample names matching rna-seq sample list, condition - treatment/ctrl, new_samples_name - new sample names
unique(metadata_df$condition_long)
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

#add gene annotations to summarize gene expression
message("adding gene annotations to dseq object")
gene_ids_from_rownames <- data.frame(gene.id = rownames(cts_coding))
merge_annot <- merge(gene_ids_from_rownames, annot, by="gene.id", all.x = TRUE)
merge_annot$basepairs <- merge_annot$length #add gene length for fpkm normalization

mcols(dds) <- DataFrame(mcols(dds), merge_annot)

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
dseq_func_out <- make_dseq_obj(ref_col = ref_variable)
dds <- dseq_func_out[["dSeqObj"]]
contrasts_ls <- dseq_func_out[["constrasts"]]