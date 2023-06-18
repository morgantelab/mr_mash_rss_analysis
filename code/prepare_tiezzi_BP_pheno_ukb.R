###Load libraries
library(data.table)

###Set seed
set.seed(1)

###Load data cleaned up by Francesco for the GxE project
load("/data2/morgante_lab/ukbiobank_projects/GxE/datasets/data3_20230511.RData")
load("/data2/morgante_lab/ukbiobank_projects/GxE/G/eigenG_20230413.RData")

###Select needed variables
eigenvec <- eigenG[[2]][1:21] ##first 20 PCs
rm(eigenG)
dttt <- dttt[, c("ID", "PP0a", "DP0a", "SP0a", "AOP", "Sex_SI")]
dat <- merge(dttt, eigenvec, by="ID", all.x=FALSE, all.y=FALSE, sort=FALSE)

###Assign fold
dat$fold <- sample(rep_len(c(1:5), nrow(dat)))

###Prepare output
indiv_list <- cbind(dat$ID, dat$ID)

###Write out output
saveRDS(dat, "../data/phenotypes/ukb_tiezzi_cleaned_BP_covar_pheno.rds")

fwrite(x=indiv_list, file="../data/misc/ukb_tiezzi_cleaned_BP_ind_ids.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)
                                 
