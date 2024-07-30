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
library(RColorBrewer)
library(NbClust)

options("width"=200)

#/juno/work/greenbaum/projects/TRI_EPIGENETIC/RNASeq_DE_TriEpi
########################################################################################
# set prroject  parameters
########################################################################################
#set wehere count data are;
source_dir = "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/RNA_seq/" 
#proj_name = "2023_BRCA_PARP_DESeq2" #project name

#BLOG: fullname path to RNA-seq data
#/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/triplicates_mouse


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
#doing this cos the file is not in the `source_dir` directory
qc_metrics <- fread("/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/triplicates_mouse/qc.tsv")

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



######
df_qc <- merge(qc_metrics, metadata_rown_df, by.x = "sample", by.y = "samples")
pltQC <- ggplot(df_qc, aes(x=sample, y=input.reads, fill=condition)) + 
geom_bar(stat="identity") + theme_minimal() + 
theme(axis.text.x = element_text(angle = 45)) + 
labs(title = "Total reads per sample", y = "Total reads") + theme(legend.position = "top")
ggsave(pltQC, file = "figures/total_reads_per_sample.pdf", width = 12, height = 12)

qc_plots <- c("input.reads", "uniquely.mapped", "total.counts", "rRNA" ,"mito.coding" ,"coding", "norm.factor",  "eff.lib.size",   "jct_fwd", "jct_rev", "unassigned_rev","QC" ,"median.degeneracy", "coding.dedup")

pdf("figures/qc_plots.pdf")
for(i in qc_plots){
  # print(i)
p <- ggplot(df_qc, aes(x=sample, y=.data[[i]], fill=condition)) + 
geom_bar(stat="identity") + theme_minimal() + 
theme(axis.text.x = element_text(angle = 45)) + 
labs(title = i, y = i) + theme(legend.position = "top")
print(p)
#ggsave(pltQC, file = "figures/total_reads_per_sample.pdf", width = 12, height = 12)
}
dev.off()
# qc_metrics %>% dplyr::left_join(metadata_rown_df, by = c("sample" , "samples"))

##exclude samples from dseq analysis
# cts_coding <- cts_coding %>% dplyr::select(!veh_2) #remove outlier samples `veh_2`
##########################
#function to make dseq object
make_dseq_obj <- function(ref_col = NULL, add_extra_contrasts = FALSE){
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
# expanded_grid[!(expanded_grid$Condition1 == expanded_grid$Condition2),]
# expanded_grid$sameCond <- expanded_grid$Condition1 == expanded_grid$Condition2


#cond_levels <- c(cond_levels, "treat_2")
cond_case_level <- cond_levels[!grepl(ref_col, cond_levels)] #take non-base levels in conditions
contrasts_ls <- paste0("condition_", cond_case_level, "_vs_", ref_col)

#add additional contrasts using pairwise comparisons
if(add_extra_contrasts){
expanded_grid <- expand.grid(Condition1 = cond_levels, Condition2 = cond_levels)
sameCond <- (expanded_grid$Condition1 == expanded_grid$Condition2)
expanded_grid <- expanded_grid[!sameCond,]
expanded_grid$contrast <- paste0("condition_", expanded_grid$Condition1, "_vs_", expanded_grid$Condition2)
contrasts_ls <- c(contrasts_ls, expanded_grid$contrast)
}

out_dseq <- list(constrasts = contrasts_ls, dSeqObj = dds)
return(out_dseq)
}
}

#create dseq object and extract contrasts
dseq_func_out <- make_dseq_obj(ref_col = ref_variable, add_extra_contrasts = TRUE)
dds <- dseq_func_out[["dSeqObj"]]
contrasts_ls <- dseq_func_out[["constrasts"]]

#from mike love, how to normalize
dds <- DESeq(dds) #this command estimates the size factors and dispersion estimates
# dds <- estimateSizeFactors(dds) #not necessary anymore 
dds_2_filterCounts <- dds
keep <- rowSums(counts(dds_2_filterCounts) >= 10) >= smallestGroupSize #keep only counts where rowSums is above 10
dds_prefiltered <- dds_2_filterCounts[keep,]
dds_prefiltered <- DESeq(dds_prefiltered) #necessary to re-run DESeq after filtering

dds_normalizedCounts <- counts(dds, normalized=TRUE) #non filtered counts
dds_normalizedCounts_dt <- as.data.table(dds_normalizedCounts) #save a copy to add gene annotations
dt_normalizedCounts_pre_filteredCounts <- counts(dds_prefiltered, normalized=TRUE)

#filter genes with at least 10 normalized counts
dds_prefiltered[rowSums(counts(dds_prefiltered)) < 1,] #check the filtered dataset
dds_2_filterCounts[rowSums(counts(dds_2_filterCounts)) < 1,]#yes, there were genes in original matrix with 0 counts

zeroCounts_FilteredCountMatirx <- apply(counts(dds_prefiltered), 2, function(x) sum(x == 0)) #check for zero values in columns
zeroCounts_orginalCountMatirx <- apply(counts(dds_2_filterCounts), 2, function(x) sum(x == 0)) #check for zero values in columns

## save zero counts data and plot
GenesWithZeroCounts_data = data.frame(samples = rep(names(zeroCounts_FilteredCountMatirx),2), 
  counts = c(zeroCounts_FilteredCountMatirx, zeroCounts_orginalCountMatirx),
  type = c(rep("filtered", 21), rep("original", 21)))

plot_zeroCounts <- ggplot() + geom_bar(data = GenesWithZeroCounts_data, aes(x = samples, y = counts, fill = type), stat = "identity", position=position_dodge()) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45)) + 
  labs(title = "Number of genes with zero counts per sample") + 
  theme(legend.position = "top")

ggsave(filename = paste0("nGenesWithZeroCounts_perSample.png"), plot_zeroCounts, width = 7, height = 10, units = "in", dpi = 300)

# GenesWithZeroCounts_data
# dt_normalizedCounts_pre_filteredCounts
# dt_normalizedCounts_pre_filteredCounts[,min(NormalizedCounts)]
# rowSums(dt_normalizedCounts_pre_filteredCounts >= 10)
# dt_normalizedCounts_pre_filteredCounts[rowSums(dt_normalizedCounts_pre_filteredCounts) >= 10]

# fpkm(dds_prefiltered)
# Transform counts for data visualization, add to tables to export
rld_filtered <- rlog(dds_prefiltered, blind=TRUE)
vst_filtered <- vst(dds_prefiltered, blind=TRUE)
rld_filtered_mat <- assay(rld_filtered)
vst_filtered_mat <- assay(vst_filtered)
fpkm_filtered <- fpkm(dds_prefiltered)


vst_filteredCounts_long <- vst_filtered_mat %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
  pivot_longer(!geneID, names_to = "samples", values_to = "VST")

rlog_filteredCounts_long <- rld_filtered_mat %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
  pivot_longer(!geneID, names_to = "samples", values_to = "Rlog")

fpkm_filteredCounts_long <- fpkm_filtered %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
  pivot_longer(!geneID, names_to = "samples", values_to = "FPKM")



# fwrite(as.data.table(dds_normalizedCounts), "data/dds_normalizedCounts_all_conditions.tsv")

# head(assays(dds)[["mu"]]) # get the mean of the normalized counts
# assays(dds)
# other normalizations
# cts_coding
# mcols(dds)
# names(mcols(dds))
cpm_out <- edgeR::cpm(dds, log = FALSE, prior.count = 2) #prior.count = 2 is added to avoid log(0) error
rpkm_out <- edgeR::rpkm(dds, gene.length = mcols(dds)$gene.length,normalized.lib.sizes = TRUE, log = FALSE, prior.count = 2) #gene.length = NULL, 

cpm_filteredCounts <- edgeR::cpm(dds_prefiltered, log = FALSE, prior.count = 2) #prior.count = 2 is added to avoid log(0) error
rpkm_filteredCounts <- edgeR::rpkm(dds_prefiltered, gene.length = mcols(dds_prefiltered)$gene.length, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 2) #gene.length = NULL, 

#convert to long format for easy plotting
cpm_filteredCounts_long <- cpm_filteredCounts %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
  pivot_longer(!geneID, names_to = "samples", values_to = "CPM")

rpkm_filteredCounts_long <- rpkm_filteredCounts %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
  pivot_longer(!geneID, names_to = "samples", values_to = "RPKM")


# joined_RPKM_CPM_df <- inner_join(cpm_filteredCounts_long, rpkm_filteredCounts_long, by = c("geneID" , "samples")) #%>% left_join(dds_normalizedCounts_long_2SAVE, by = c("geneID" , "samples"))

#merge all normalized counts and transformations for visulazation
rnaSeqMatrices_list <- list(rlog = rlog_filteredCounts_long, vst = vst_filteredCounts_long,
cpm = cpm_filteredCounts_long, rpkm=rpkm_filteredCounts_long, fpkm = fpkm_filteredCounts_long)
# rpkm_out2 <- edgeR::rpkm(dds, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 2) #gene.length = NULL, 

rnaSeqMatrices_long_df <- reduce(rnaSeqMatrices_list, left_join, by = c("geneID" , "samples"))

# map(rnaSeqMatrices, function(x){
#   reduce(left_join, by = c("geneID" , "samples"))
# })
# head(rpkm_out); head(rpkm_out2)
# ENSMUSG00000000001.4  56.242237 59.322645 55.262150 50.71634920 51.304924
# ENSMUSG00000000001.4  47.91211228 49.993275 50.27821 43.006207 53.719418
# ENSMUSG00000000001.4  55.13339537 55.12542128 53.322225 54.92939288 54.562924
#cpm_out2 <- edgeR::cpm(counts(dds), log = FALSE, prior.count = 2)
# rowRanges(dds)

dds_normalizedCounts_long_2SAVE <- dt_normalizedCounts_pre_filteredCounts %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
  pivot_longer(!geneID, names_to = "samples", values_to = "NormalizedCounts") %>%  
  left_join(as.data.frame(metadata_rown_df) %>% mutate(samples = str_replace_all(samples, "\\.", "-"))) %>% 
  left_join(annot %>% dplyr::select(c(`gene.id`,`gene.symbol`, length)), by = c("geneID" = "gene.id")) #%>%
#mutate(geneID = paste0(`gene.symbol`,"_",geneID))

dds_normalizedCounts_long_2SAVE %>% filter(NormalizedCounts == 0) #%>% nrow() #no zero values
dds_normalizedCounts_long_2SAVE <- left_join(dds_normalizedCounts_long_2SAVE, rnaSeqMatrices_long_df, by = c("geneID" , "samples")) %>% mutate(geneID = paste0(`gene.symbol`,"_",geneID)) 
# dds_normalizedCounts_long_2SAVE <- left_join(dds_normalizedCounts_long_2SAVE, joined_RPKM_CPM_df, by = c("geneID" , "samples")) %>% mutate(geneID = paste0(`gene.symbol`,"_",geneID)) 

#double check, no zero values
dds_normalizedCounts_long_2SAVE %>% summarize(min(NormalizedCounts), max(NormalizedCounts), mean(NormalizedCounts), sd(NormalizedCounts))
# joined_RPKM_CPM_df
fwrite(as.data.table(dds_normalizedCounts_long_2SAVE), "data/dds_normalizedCounts_with_annot_all_samples.tsv")
# fwrite(as.data.table(dds_normalizedCounts_long_2SAVE), "data/dds_normalizedCounts_with_annot_all_samples.tsv")



################################################################################################
#plot gene expressions
#1. all coding genes per sample; sample page
#first make  a long format of the data for easy ggplots; this is log format
dds_normalizedCounts_long <- dds_normalizedCounts %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
  pivot_longer(!geneID, names_to = "samples", values_to = "NormalizedCounts") %>% 
      mutate(NormalizedCounts = log(0.1 + NormalizedCounts))

#add sample metadata 
dds_normalizedCounts_long_metadata <- dds_normalizedCounts_long %>% 
  left_join(as.data.frame(metadata_rown_df) %>% 
    mutate(samples = str_replace_all(samples, "\\.", "-"))) 
#save a copy of the normalized counts with metadata
data2Save <- dds_normalizedCounts_long_metadata %>% left_join(annot %>% dplyr::select(c(`gene.id`,`gene.symbol`)), by = c("geneID" = "gene.id")) %>%
mutate(geneID = paste0(`gene.symbol`,"_",geneID))
# fwrite(as.data.table(data2Save), "data/dds_normalizedCounts_with_annot_all_samples.tsv")


#### 


#head(dds_normalizedCounts) 
#filter genes with at least 10 normalized counts
dds_normalizedCounts_long_min10 <- dds_normalizedCounts %>% as.data.frame() %>% 
  rownames_to_column(var = "geneID") %>% 
    pivot_longer(!geneID, names_to = "samples", values_to = "NormalizedCounts") %>% 
      dplyr::filter(NormalizedCounts >= 10) %>% 
      mutate(NormalizedCounts = log(NormalizedCounts))
      
stats_df <-  dds_normalizedCounts_long_min10 %>% 
      group_by(samples) |>
      summarize(n=n(), med= round(median(NormalizedCounts), 1), mu= round(mean(NormalizedCounts), 1), 
        sd=round(sd(NormalizedCounts), 1)) |> mutate(stats = paste0("med-", med," mu-" ,mu, " sd-", sd)) |> dplyr::select(!c(n,med,  mu, sd))
stats_df$stats <- str_wrap(stats_df$stats, width = 10)  # Adjust 'width' as needed




dds_normalizedCounts_long |> group_by(samples) |>
      summarize(n=n(), min(NormalizedCounts), med= round(median(NormalizedCounts), 1), mu= round(mean(NormalizedCounts), 1), 
        sd=round(sd(NormalizedCounts), 1))

plot_NormalizedCounts <- ggplot(dds_normalizedCounts_long, aes(NormalizedCounts)) + geom_freqpoly() + 
  facet_wrap(~samples) + labs(title = "log(0.1+Normalized read counts) per sample including 0 Expr") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45))
ggsave(plot_NormalizedCounts, file = "figures/freqpoly_NormalizedCounts.pdf") 

#plot samples gene expression with min10 genes expressed
plot_NormalizedCounts_min10 <- ggplot(dds_normalizedCounts_long_min10, aes(NormalizedCounts)) + 
  geom_freqpoly() + facet_wrap(~samples) + 
  geom_text(data=stats_df, mapping = aes(x = -Inf, y = -Inf, label = stats),
  hjust   = -0.1, vjust   = -1, size=3 ) + 
  labs(title = "log(>=10 normalized read counts) per sample", y = "counts") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45))
ggsave(plot_NormalizedCounts_min10, file = "figures/freqpoly_NormalizedCounts_min10.pdf") 



#dds_normalizedCounts_long_metadata %>% dplyr::group_by(condition, samples) %>% summarize(n())

plot_NormalizedCounts_groups <- ggplot(dds_normalizedCounts_long_metadata, 
  aes(NormalizedCounts, color = samples)) + 
    geom_freqpoly() +
    facet_wrap(~condition) +
  scale_fill_brewer(palette = "Dark2", name = "condition") + 
  labs(title = "log(0.1 + normalized read counts) per sample", y = "counts") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) 
ggsave(plot_NormalizedCounts_groups, file = "figures/freqpoly_NormalizedCounts_conditions.pdf") 
#aes(alpha = samples)

plot_NormalizedCounts_groupsall <- ggplot(dds_normalizedCounts_long_metadata, 
  aes(NormalizedCounts, color = condition)) + 
    geom_freqpoly(aes(alpha = samples)) + 
  scale_fill_brewer(palette = "Dark2", name = "condition") + 
  labs(title = "log(0.1+Normalized read counts)") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))
ggsave(plot_NormalizedCounts_groupsall, file = "figures/freqpoly_NormalizedCounts_all_conditions.pdf") 


##plot each genes expression 
data_filtered <- dds_normalizedCounts_long_metadata %>% dplyr::filter(NormalizedCounts >= log(10)) 
data_filtered <- data_filtered %>% left_join(annot %>% dplyr::select(c(`gene.id`,`gene.symbol`)), by = c("geneID" = "gene.id")) %>%
mutate(geneID = paste0(`gene.symbol`,"_",geneID))





plot_genes <- function(geneID){
  data_in <- data_filtered %>% dplyr::filter(geneID == {{geneID}})
print(geneID)
plt_point <- ggplot(data = data_in) + 
geom_point(aes(x=samples, y=NormalizedCounts, color=condition)) + labs(title = geneID)+
scale_colour_brewer(palette="Dark2")+
labs(title = paste0("",geneID),
y = "log(0.1 + >10 Normalized Counts)")+
theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

plt_hist <- ggplot(data = data_in %>% 
dplyr::filter(geneID == geneID)) + 
geom_histogram(aes(NormalizedCounts)) + 
labs(title = paste0("Per Gene dist across samples"), x = "log(Normalized Counts)") +
#labs(title = "ENSMUSG00000000001.4; all samples") +
theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

#merge plots for visualization
plt_point_hist <- plot_grid(plt_point, plt_hist)
ggsave(plt_point_hist, file = paste0("figures/geneWiseNormalizedCounts/",geneID,"_point_hist.pdf"), width = 9 , height = 5) 
}

#plot for all genes
#geneList <- head(unique(data_filtered$geneID))
# geneList <- unique(data_filtered$geneID)
# walk(geneList, function(x){
# plot_genes(geneID=x)
# })

# data_filtered %>% dplyr::filter(geneID == "ENSMUSG00000000049.11")

#plot means of each gene
perGene_Stats <- data_filtered %>% group_by(geneID) %>% 
        summarize(n=n(), med= round(median(NormalizedCounts), 1), 
                  mu= round(mean(NormalizedCounts), 1), 
                  sd=round(sd(NormalizedCounts), 1)) |> 
        mutate(stats = paste0("med-", med," mu-" ,mu, " sd-", sd)) 

#convert to longer for scater plots        
perGene_Stats_longer <- perGene_Stats  |> 
        pivot_longer(!c(geneID,stats))
        #dplyr::select(!c(n, med, mu, sd))
perGene_Stats_longer$stats <- str_wrap(perGene_Stats_longer$stats, width = 10)  # Adjust 'width' as needed

#perGene_Stats_longer  |> filter(name == "n") %>% summarize(min(value))


plot_aggregate_conditions_normalizedCounts <- ggplot(perGene_Stats_longer, aes(value)) + 
geom_histogram() + 
labs(title  = "per gene statistics across samples and conditions", subtitle = "log(panels) !n") + 
facet_wrap(~name, scales = "free") + 
theme_minimal()
ggsave(plot_aggregate_conditions_normalizedCounts, file = paste0("figures/","histogram_aggregate_conditions_normalizedCounts_point.pdf"), width = 9 , height = 5) 


plot_sd_mean_per_gene_across_samples <- ggplot(perGene_Stats, aes(x=mu, y = sd)) + 
geom_point() + labs(title  = "mu vrs sd per gene statistics across samples and conditions",
x = "log(Gene mean across samples)",
y  = "log(Gene sd across samples)") + 
geom_smooth() +
theme_minimal()
ggsave(plot_sd_mean_per_gene_across_samples, file = paste0("figures/","mean_vrs_sd_perGene_across_samples_conditions.pdf"), width = 9 , height = 5) 


plot_sd_median_per_gene_across_samples <- ggplot(perGene_Stats, aes(x=med, y = sd)) + 
geom_point() + labs(title  = "median vrs sd per gene statistics across samples and conditions",
x = "log(Gene median across samples)",
y  = "log(Gene sd across samples)") + 
geom_smooth() +
theme_minimal()
ggsave(plot_sd_median_per_gene_across_samples, file = paste0("figures/","median_vrs_sd_perGene_across_samples_conditions.pdf"), width = 9 , height = 5) 




################################################################################################
################################################################################################
##pre - filtering
# rowSums(counts(dds)) <= 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
# keep <- rowSums(counts(dds)) >= 10 #keep only counts where rowSums is above 10
dds <- dds[keep,] #filter dds
# head(counts(dds))
names(assays(dds))

# keep <- rowSums(counts(dds_2_filterCounts) >= 10) >= smallestGroupSize #keep only counts where rowSums is above 10

################################################################################################
############################### Step 2; Run deseq2 analysis
##############################################################################################
#run DESeq analysis - normalization and filtering
dds <- DESeq(dds)
names(mcols(dds))
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

# Transform counts for data visualization
rld <- rlog(dds, blind=blind_transform)
# rld_mat <- assay(rld)

# Plot PCA 
pltPCA_samples <- plotPCA(rld, intgroup=c("condition", "samples"))
pltPCA_ConditionsOnly <- plotPCA(rld, intgroup=c("condition"))
ggsave(pltPCA_samples, file = paste0("figures/PCA_by_condition_and_samples.pdf"), width = 12, height = 12)
ggsave(pltPCA_ConditionsOnly, file = paste0("figures/PCA_by_conditionOnly.pdf"), width = 12, height = 12)

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
rld_df <- rld_mat %>% as.data.frame() %>% rownames_to_column("gene.id")

#save a copy of the rlog matrix
fwrite(rld_df, file=paste0("data/DESeq2_rlog_Transform_Blind",blind_transform,".tsv"), sep="\t")


# Plot heatmaps
pltHeatmap_samples <- pheatmap(rld_cor, annotation = metadata_rown_df[,c("condition"), drop = FALSE], scale = "row")
ggsave(pltHeatmap_samples, file = paste0("figures/heatmap_conditions.pdf"), width = 12, height = 12)


# Compute the distance matrix of samples
sampleDists <- as.dist(1 - cor(rld_mat)) #what does, closet to 1 mean? similar samples
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, annotation = metadata_rown_df[,c("condition"), drop = FALSE] ,filename = "figures/sampleDistMatrix.pdf")



res_rld_mat_df <- rld_mat %>%
  as.data.frame() %>%
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
# contrast_input <- "condition_AZCT_vs_DMSO" #uncomment to test run function

# metadata(dds_copy)
# elementMetadata(dds_copy)
# dseqObject <- dds_copy[,-idx]
################################################################################################################################################
######################################################################## LRT geneclustering ###################################################
#prepare data for clustering
dds2 <- DESeq(dds, test = "LRT", reduced = ~1)
acrossGroups <- results(dds2)
acrossGroups <- acrossGroups[order(acrossGroups$pvalue), ]
acrossGroups[1:3, ]

rld_4dds2 <- rlog(dds2, blind=FALSE)
rld_4dds2_mat <- assay(rld_4dds2)

sigChanges <- rownames(acrossGroups)[acrossGroups$padj < 0.01 & !is.na(acrossGroups$padj)]
sigMat <- rld_4dds2_mat[rownames(rld_4dds2_mat) %in% sigChanges, ]
pheatmap(sigMat, scale = "row", show_rownames = FALSE, annotation = metadata_rown_df[,c("condition"), drop = FALSE] ,filename = "figures/geneClustering_LRT_full_vrs_reducedModel.pdf")
k <- pheatmap(sigMat, scale = "row", kmeans_k = 7, annotation = metadata_rown_df[,c("condition"), drop = FALSE] ,filename = "figures/geneClustering_LRT_full_vrs_reducedModel_kmeans7.pdf")

# names(k$kmeans)
clusterDF <- as.data.frame(factor(k$kmeans$cluster))
colnames(clusterDF) <- "Cluster"
clusterDF[1:10, , drop = FALSE]
head(clusterDF)

OrderByCluster <- sigMat[order(clusterDF$Cluster), ]

# Specify colors
# ann_colors = list(Cluster = brewer.pal(n = length(levels(clusterDF$Cluster)), name = 'Dark2'))


pheatmap(OrderByCluster, scale = "row", annotation_row = clusterDF, show_rownames = FALSE,
    cluster_rows = FALSE,
    annotation = metadata_rown_df[,c("condition"), drop = FALSE], 
    # annotation_colors = ann_colors,
    filename = "figures/geneClustering_LRT_full_vrs_reducedModel_kmeans7_withClutersAnnot.pdf")


# find the optimal clusters
rowScaledMat <- t(scale(t(sigMat)))
clusterNum <- NbClust(rowScaledMat, distance = "euclidean", min.nc = 2, max.nc = 12,
    method = "kmeans", index = "silhouette")
clusterNum$Best.nc


# merge automatic clustering with manual clustering
clusterNum$Best.partition[1:10]
orderedCluster <- sort(clusterNum$Best.partition)
sigMat <- sigMat[match(names(orderedCluster), rownames(sigMat)), ]
# length(clusterNum$Best.partition)

pheatmap(sigMat, scale = "row", annotation_row = clusterDF, show_rownames = FALSE,
    cluster_rows = FALSE,
    annotation = metadata_rown_df[,c("condition"), drop = FALSE], 
    filename = "figures/geneClustering_LRT_full_vrs_reducedModel_kmeans_withOptimumClustering.pdf")


################################################################################################################################################
# nrow(rld_4dds2_mat)
#defining the optimal clusters





#####################################################################################################################
############################### Step 3; Run deseq2 analysis for each contrast
#####################################################################################################################



DEResults_ls <- list()
source("../scripts/run_multiple_contrasts.R") #split this script into a function
################################################################################################
############################### Step #; Run DE for each contrast
##############################################################################################
#run for all contrasts
#contrast_input="condition_AZCT_vs_DMSO"

DEResults_ls <- lapply(contrasts_ls, function(x) {run_multiple_contrasts(contrast_input=x, dseqObject = dds, shrink = FALSE, export_transformed_counts = TRUE)})
names(DEResults_ls) <- contrasts_ls
#summary(DEResults_ls[[1]])
# data.frame(DEResults_ls[[1]]) %>% dplyr::filter(padj < 0.05) %>% dim()
# contrast_input=contrasts_ls[1]









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
# h_gsea <- lapply(contrasts_ls, function(x) gsea_analysis(x, specie_type="Mus musculus",category_tag = "H", nCategory_2show = 100, ggwidth = 15, ggheight = 12))
# c2_gsea <- lapply(contrasts_ls, function(x) gsea_analysis(x, specie_type="Mus musculus",do_gse_GO=FALSE, category_tag = "C2", ggwidth = 15, ggheight = 12))
# names(h_gsea)
# names(h_gsea[[1]][["rankedFC"]])

# c2_gsea[[1]]
# head(h_gsea[[1]]@result)
# slotnames(h_gsea)
# slotNames(h_gsea)

