library(data.table)
library(bigreadr)

ukb_data_location <- "/data2/morgante_lab/data/ukbiobank/genotypes/imputed/"
ids <- bigreadr::fread2(paste0(ukb_data_location, "ukb22828_c1_b0_v3_s487271.sample"))
ids <- ids[-1, ]

ids1 <- data.table::fread("../data/misc/ukb_cleaned_bc_ind_ids.txt", data.table=FALSE, 
                          colClasses = "character", header=FALSE)

for(fold in 1:5){
  samples <- paste0("../data/misc/ukb_cleaned_bc_ind_ids_plink_", fold, ".txt")
  samples_to_remove <- data.table::fread(samples, data.table=FALSE, colClasses = "character", 
                                         header=FALSE)
  
  if(ncol(samples_to_remove) > 1){
    samples_to_remove <- samples_to_remove[, 2]
  }
  
  #Keep only bc selected individuals
  ids_sel <- as.character(ids[as.character(ids$ID_2) %in% ids1[,2], 2])
  
  #Keep only training individuals  
  ids_sel_train <- ids_sel[!(ids_sel %in% samples_to_remove)]
  
  #Write out file
  write.table(ids_sel_train, paste0("../data/LDstore_input_files/ukb_bc_training_samples_LDstore_", fold, ".incl"),
              col.names = FALSE, row.names = FALSE, quote=FALSE, sep=" ")
}

