# Index BAM files
# Usage:
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/RNA-seq_DiffExpr/scripts/indexBams.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_BulkONT_Cancer_wf/config/cluster_profiles/slurm --jobs 10 --cores all --keep-going --forceall -np

SAMPLES = [
    "R-SC-1_IGO_14608_B_19",
    "R-S-3_IGO_14608_B_12",
    "R-S-2_IGO_14608_B_11",
    "R-S-1_IGO_14608_B_10",
    "R-QC-3_IGO_14608_B_18",
    "R-QC-2_IGO_14608_B_17",
    "R-QC-1_IGO_14608_B_16",
    "R-Q-3_IGO_14608_B_9",
    "R-Q-2_IGO_14608_B_8",
    "R-Q-1_IGO_14608_B_7",
    "R-C-3_IGO_14608_B_15",
    "R-C-2_IGO_14608_B_14",
    "R-C-1_IGO_14608_B_13"
]

# SAMPLES = [
#     "R-A-2_IGO_14608_B_5",
#     "R-A-1_IGO_14608_B_4",
#     "R-0-3_IGO_14608_B_3",
#     "R-0-2_IGO_14608_B_2",
#     "R-0-1_IGO_14608_B_1"
# ]

rule all:
    input:
        expand("/data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/pipeline/{sample}.bam.bai", sample=SAMPLES)

rule index_bam:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    threads: 12
    shell:
        "samtools index -@ {threads} {input}"

