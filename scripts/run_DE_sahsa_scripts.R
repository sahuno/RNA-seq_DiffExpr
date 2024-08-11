#how to run shasha de code
# dir of files; /data1/greenbab/users/ahunos/apps/rnaseq
run_de_analysis.R
--fc-threshold 0.05
--fdr-cutoff 0.05
--num-genes-go 50
--out-dir 
--min-row-sum
--targets

############################################################
############################################################
#shasha pipelines
module load rna-seq
module load singularity/3.7.1

rnsq-run-container.sh

#init project
rnsq-init.sh
#since i don't have raw files in same dir as the project, i will create a symbolic link to the raw files
ln -s /data1/greenbab/projects/public_Epigenetic_therapy/HCT116_AZA_RNA/fastq /data1/greenbab/projects/public_Epigenetic_therapy/HCT116_AZA_RNA/raw

#run the pipeline
rnsq-mkrawsamplelist.sh
rnsq-mksamplelist.sh
rnsq-run-mapping.sh --time=36:00:00

### downstream analysis
rnsq-join_counts.sh
run_heatmap.R #heatmap
rnsq-write_expr.R
rnsq-de.sh --targets sample_list_with_phenotype.tsv
############################################################
############################################################
############################################################
############################################################

phen <- "/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/ref_alt_table.tsv"
targets <- read.delim("sample_list_with_phenotype.tsv")
targets <- read.delim(phen)

ct.all <- ct[,targets$sample]
read.del



#########
#for human samples
rnsq-de.sh --targets sample_list_with_REFnALTphenotype.tsv
