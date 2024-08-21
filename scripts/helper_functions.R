#utils dunctions
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



## plot individul gene expression 
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