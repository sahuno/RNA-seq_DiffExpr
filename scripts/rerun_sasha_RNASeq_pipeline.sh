############################################################
#author: samuel ahuno
#date: 2025-03-12
#purpose: notes on how to rerun the sasha RNASeq pipeline
#project: triplicates_epigenetics_diyva
############################################################
#shasha pipelines
module load rna-seq
module load singularity/3.7.1

#working direcyory
cd /data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025

#link raw reads (fastq files) to the project
ln -s /data1/greenbab/junobackup/projects/divya_epigenetics/mouse_triplicates/rna/Project_14608_B/raw/ /data1/greenbab/projects/triplicates_epigenetics_diyva/RNA/rerun_RNASeq_11032025/raw

##start the container
rnsq-run-container.sh

#init project, sets the root directory to the currrent directory
rnsq-init.sh

#run the pipeline
rnsq-mkrawsamplelist.sh
rnsq-mksamplelist.sh

###edited config.sh

## now run mapping, this should be run outside the container 
rnsq-run-mapping.sh --time=48:00:00
