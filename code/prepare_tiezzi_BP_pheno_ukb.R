###Load libraries
library(data.table)

###Set seed
set.seed(1)

###Load data cleaned up by Francesco for the GxE project and the ids of the imputed genotype data
load("/data2/morgante_lab/ukbiobank_projects/GxE/datasets/data2_20230413.RData")
bt <- readRDS("/data2/morgante_lab/ukbiobank_projects/mr_mash_rss/data/raw/ukb_blood_traits.rds")
load("/data2/morgante_lab/ukbiobank_projects/GxE/G/eigenG_20230413.RData")
ids <- fread("/data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c1_b0_v3_s487271.sample",
             showProgress = FALSE, data.table = FALSE)
ids <- ids[-1, ]

###Select needed variables
eigenvec <- eigenG[[2]][1:21] ##first 20 PCs
rm(eigenG)

dttt <- merge(dtt, bt, by="ID", all.x=FALSE, all.y=FALSE, sort=FALSE)

dttt <- dttt[, c("ID", "Basophill_perc", "Eosinophill_perc", "Haemoglobin_conc", "HLR_perc", 
               "Lymphocyte_perc", "MCV", "Monocyte_perc", "MSCV", "Neutrophill_perc", "Platelet_dw", 
               "Platelet_count", "Platelet_crit", "RBC_count", "RBC_dw", "Reticulocyte_perc", 
               "WBC_count", "AOP", "Sex_SI")]
dttt <- dttt[complete.cases(dttt), ]

dat <- merge(dttt, eigenvec, by="ID", all.x=FALSE, all.y=FALSE, sort=FALSE)

###Subset individuals with imputed genotype data
dat <- dat[which(as.character(dat$ID) %in% as.character(ids[,2])), ]

###Assign fold
dat$fold <- sample(rep_len(c(1:5), nrow(dat)))

###Prepare output
indiv_list <- cbind(dat$ID, dat$ID)

###Write out output
saveRDS(dat, "../data/phenotypes/ukb_tiezzi_cleaned_BP_covar_pheno.rds")

fwrite(x=indiv_list, file="../data/misc/ukb_tiezzi_cleaned_BP_ind_ids.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)
                                 
