name: samuel ahuno
date: august 25 2023

# differential expression analysis from counts and annotation table

How to run?
Run DESEQ2 for coding genes
Rscript scripts/run_multiple_contrasts

```
Rscript RNA-seq_DiffExpr/scripts/coding_genes.r \
  --source_dir "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/" \
  --workflow_dir "/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/" \
  --blind_transform TRUE \
  --drop_samples "R.S.2,R.C.3" \
  --ref_variable "DMSO" \
  --min_read_counts 50 \
  --smallest_group_size 3 \
  --qc_metrics "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/triplicates_mouse/qc.tsv"


```
