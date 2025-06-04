###

# dds <- estimateSizeFactors(dds) #not necessary anymore 
dds_2_filterCounts <- dds
keep <- rowSums(counts(dds_2_filterCounts) >= MIN_ReadsCounts) >= smallestGroupSize #keep only counts where rowSums is above 10
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
  type = c(rep("filtered", length(colnames(dds_prefiltered))), rep("original", length(colnames(dds_prefiltered)))))

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

compute_tpm <- function(counts, gene_lengths) {
  # Convert raw counts to RPK (Reads Per Kilobase)
  rpk <- counts / gene_lengths
  # Compute scaling factor (per sample)
  scaling_factors <- colSums(rpk, na.rm = TRUE)
  # Calculate TPM
  tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6
  return(tpm)
}

tpm_matrix_filtered <- compute_tpm(counts(dds_prefiltered), gene_lengths=mcols(dds_prefiltered)$length)
tpm_filteredCounts_long <- tpm_matrix_filtered %>% 
  as.data.frame() %>% rownames_to_column(var = "geneID") %>% 
  pivot_longer(!geneID, names_to = "samples", values_to = "TPM")

#tpmTest <- counts(dds)/mcols(dds)$length
#scaling_factors2 <- colSums(tpmTest, na.rm = TRUE)

#merge all normalized counts and transformations for visulazation
rnaSeqMatrices_list <- list(rlog = rlog_filteredCounts_long, vst = vst_filteredCounts_long,
cpm = cpm_filteredCounts_long, rpkm=rpkm_filteredCounts_long, fpkm = fpkm_filteredCounts_long, tpm = tpm_filteredCounts_long)
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

