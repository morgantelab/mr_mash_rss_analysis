library(data.table)
library(dplyr)
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--n_fold"), type="integer")
outparse <- parse_args(parser)

n_fold <- outparse$n_fold

set.seed(1)

###Read in BC data
dat <- readRDS("../data/phenotypes/ukb_cleaned_bc_covar_pheno.rds")

dat_sampled <- dat %>% group_by(fold) %>% slice_sample(n=n_fold) %>% ungroup() %>% 
  arrange(as.numeric(id)) %>% as.data.frame()

###Prepare output
indiv_list <- cbind(dat_sampled$id, dat_sampled$id)

###Write out output
saveRDS(dat_sampled, "../data/phenotypes/ukb_cleaned_bc_sampled_covar_pheno.rds")

fwrite(x=indiv_list, file="../data/misc/ukb_cleaned_bc_sampled_ind_ids.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)
