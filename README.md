name: samuel ahuno
date: august 25 2023

# differential expression analysis from counts and annotation table

How to run?

1. coding_genes.r (Main DE analysis)
   ↓ (produces individual contrast results)
2. run_GSEA_conditionSpecific.R → Individual condition GSEA (this script!)
3. combine_results_of_mutiple_conditons.R (This script - Meta-analysis)
   ↓ (combines and compares all contrasts)
4. Final integrated results


Run DESEQ2 for coding genes
writes results into current directory

mkdir -p codingGenes_DE
```
Rscript RNA-seq_DiffExpr/scripts/coding_genes.r \
  --source_dir "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025" \
  --workflow_dir "/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/" \
  --blind_transform TRUE \
  --drop_samples "R.S.2,R.C.3" \
  --ref_variable "DMSO" \
  --min_read_counts 50 \
  --smallest_group_size 3 \
  --qc_metrics "/data1/greenbab/projects/methylRNA/Methyl2Expression/data/preprocessed/triplicates_mouse/qc.tsv"
```


# meta-analysis of all DE DESEQ2 results
```
Rscript /data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/scripts/combine_results_of_mutiple_conditons.R \
  --specie_type "Mus musculus" \
  --category_tag "C2" \
  --fractEnriched 0.5 \
  --do_gse_GO \
  --nCategory_2show 10 \
  --nShowBar 100 \
  --ggwidth 13 \
  --ggheight 11 \
  --results_dir "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/codingGenesDExpr/data"
```  


gene set enrichment analysis with DE data directly from sahas pipeline
```
/data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/scripts/geneSet_enrichment_for_sasha_pipeline.R
```