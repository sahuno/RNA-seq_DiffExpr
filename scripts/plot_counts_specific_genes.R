library(data.table)
library(tidyverse)

cntsDT <- fread("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/codingGenes_DE/data/dds_normalizedCounts_with_annot_all_samples.tsv")

cntsDT
