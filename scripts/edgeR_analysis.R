library(edgeR)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

rawdata <- read.delim("/Users/ahunos/myWork/research/Greenbaum/data/preprocessed/counts_annot.tsv", check.names=FALSE, stringsAsFactors=FALSE)
y <- DGEList(counts=rawdata[,11:ncol(rawdata)], genes=rawdata[,1:11])
head(rawdata)
names(rawdata)

rawdata %>% filter(gene.type %in% "protein_coding") %>% filter(str_detect(description, "X-linked") )
rawdata %>% filter(gene.type %in% "protein_coding") %>% filter(str_detect(description, "x") )

#add metadata
metadata <- read.delim("/Users/ahunos/myWork/research/Greenbaum/data/preprocessed/metadata_triplicates.csv", sep = ",",check.names=FALSE, stringsAsFactors=FALSE)
y$samples$condition <- metadata$condition
unique(metadata$condition_long)
#idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Mm.egREFSEQ)
#y <- y[idfound,]
dim(y)


##filtering and normalization
names(y$genes)[4] <-  "Symbol" #rename "gene.symbol"

o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$entrez.gene.id)
#d <- duplicated(y$genes$Symbol)
sum(duplicated(y$genes$entrez.gene.id))
y <- y[!d,]
nrow(y)

#add lib size
y$samples$lib.size <- colSums(y$counts)
#y$genes$entrez.gene.id

#remove all genes without ids
y <- y[!is.na(y$genes$entrez.gene.id),]
rownames(y$counts) <- rownames(y$genes) <- y$genes$gene.id
#y$genes$entrez.gene.id <- NULL
y <- normLibSizes(y)

d0 <- calcNormFactors(y)
dim(d0)

##filter cpm below 1.
#rows should be greater than 1
cutoff <- 1
drop_idx <- which(apply(cpm(d0), MARGIN=1, max) < cutoff)
cpm_filtered_genes <- d0[-drop_idx,] 
dim(cpm_filtered_genes) # number of genes left
plotMDS(cpm_filtered_genes)

###export count per millions for use in other downstream analysis
cpms_sample <- cpm(cpm_filtered_genes, log=FALSE)
rpkm_sample <-  rpkm(cpm_filtered_genes, gene.length = cpm_filtered_genes$genes$length, normalized.lib.sizes = TRUE, log = TRUE)

rpkm_sample_ls <- lapply(seq_len(ncol(rpkm_sample)), function(i) {
  df <- data.frame(rpkm = rpkm_sample[,i])
  df <-  df %>% mutate(tiles = ntile(., 10))})
names(rpkm_sample_ls) <- colnames(rpkm_sample)
#head(rpkm_sample_ls[[1]])

#rpkm_sample_ls[["R.0.1"]]
##imap(names(rpkm_sample_ls), ~paste0(.x, .y))
iwalk(names(rpkm_sample_ls), ~write_csv(rpkm_sample_ls[[.x]], paste0(.x, '.csv')))


library(data.table)
rpkm_sample_df <- as.data.frame(rpkm_sample)
rpkm_sample_df <- rpkm_sample_df %>% rownames_to_column(var = "geneID")
#rpkm_sample_df$geneID <- rownames(rpkm_sample_df)
fwrite(cpms_sample, file="CPM_tri_epigenetic.csv")
fwrite(rpkm_sample_df, file="rpkm_sample_tri_epigenetic.csv")
#fwrite(rpkm_sample, file="rpkm_sample_tri_epigenetic.csv")

head(rpkm_sample_df)

rpkm_sample_long_df <- rpkm_sample_df %>% pivot_longer(!geneID)
rpkm_sample_long_df_stats <- rpkm_sample_long_df %>% group_by(geneID) %>% summarise(n = dplyr::n(), meanRpKM = mean(value, na.rm=TRUE),medianRpKM = median(value, na.rm=TRUE) ,sdRpKM = sd(value, na.rm=TRUE)) %>% ungroup()
rpkm_sample_df_stats_long <- rpkm_sample_long_df_stats  %>% pivot_longer(!geneID)

rpkm_sample_long_df_stats  %>% arrange(sdRpKM, meanRpKM)
rpkm_sample_df_stats_long %>% dplyr::filter(name == "sdRpKM") %>% arrange(desc(value))
rpkm_sample_df_stats_long %>% dplyr::filter(!name %in% c("n")) %>% arrange(value)
rpkm_sample_df_stats_long %>% dplyr::filter(name == "sdRpKM") %>% summarise(min = min(value, na.rm=TRUE)) 
rpkm_sample_df_stats_long %>% dplyr::filter(name == "sdRpKM") %>% summarise(max = max(value, na.rm=TRUE)) 

 # arrange(desc(value))

pp <- ggplot(rpkm_sample_df_stats_long, aes(value)) + 
  geom_histogram() + labs(title = "rpkm distribution across samples") +
  facet_wrap(~name, scales = "free") + theme_minimal()
ggsave(pp, filename = "/Users/ahunos/myWork/research/Greenbaum/OneOnOne_04242024/distribution_rpkm.pdf")


ppp_point <- ggplot(rpkm_sample_long_df_stats, aes(meanRpKM, sdRpKM)) + 
  geom_point() + labs(title = "mean vrs sd across samples") +
  theme_minimal()
ggsave(ppp_point, filename = "/Users/ahunos/myWork/research/Greenbaum/OneOnOne_04242024/mean_sd_rpkm.png", bg = "white")


####################################
##do rna-seq with EdgeR
####################################
design <- model.matrix(~0 + metadata$condition)
designEdgeR <- model.matrix(~metadata$condition)
rownames(designEdgeR) <- colnames(y)

rownames(design) <- colnames(y)
plotMDS(d, col = as.numeric(group))


### edgeR DE
d <- calcNormFactors(d)
dge.estDisp.out <- estimateDisp(d, design = designEdgeR)
fit <- glmFit(dge.estDisp.out, design = designEdgeR)
lrt <- glmLRT(fit, 3)
topTags(lrt)



########################################################################
########################################################################
########################################################################
## 3. Voom transformation and calculation of variance weights
y <- voom(d, design=design, plot = T)
fit <- lmFit(y, design=design)
head(coef(fit))
cpm(fit)
head(cpm(fit))[1:6, 1:6]
view(head(cpm(fit))[1:6, 1:6])
#compare to previous
tmp <- voom(d0, design=design, plot = T)


fit_eBayes <- eBayes(fit)
plotSA(fit_eBayes, main="Final model: Mean-variance trend")
head(y$E)

