#make sample metadata 
#note; we assume you already have sample metadata file saved some where
#read coldata 
samples_df <- read.delim(paste0(source_dir,'sample_list'), header = FALSE)
names(samples_df) <- "samples"

#create samepl;e names the match those in counts data sets
trt <- paste0(rep("treat_", 3), c(1,2,3))
ctrls <- paste0(rep("veh_", 3), c(1,2,3))

#make condition column
cond_trt <- paste0(rep("PARPi", 3))
cond_ctrls <- paste0(rep("CTRL", 3))

#new sample names
new_cases <- paste0(rep("Parp_", 3), c(1,2,3))
new_ctrls <- paste0(rep("Ctrl_", 3), c(4,5, 6))

metadata_df <- data.frame(samples = c(trt, ctrls), condition= c(cond_trt, cond_ctrls), new_samples_name= c(new_cases, new_ctrls) ,stringsAsFactors = FALSE)
metadata_df

write.table(metadata_df, file = paste0( "data/sample_metadata.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE)

