###Load libraries needed
library(bigsnpr)
library(data.table)
library(RhpcBLASctl)

###Set MKL threads
RhpcBLASctl::blas_set_num_threads(1)

###Set seed
set.seed(1)

###Load data
geno_imp <- snp_attach("/scratch1/fabiom/ukb_bc_geno_imp_HM3.rds")
pheno_dat <- readRDS("../data/phenotypes/ukb_cleaned_bc_covar_pheno.rds")

###Extract variant names
variant_names <- data.frame(geno_imp$map$marker.ID)
variant_names_rsid <- data.frame(geno_imp$map$rsid)

###Extract individuals name by fold
for(i in 1:5){
  ids <- pheno_dat[which(pheno_dat$fold == i), "id"]
  dat <- data.frame(ids, ids)
  
  fwrite(x=dat, file=paste0("../data/misc/ukb_cleaned_bc_ind_ids_plink_", i, ".txt"), 
         sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)
}

###Write out output
fwrite(x=variant_names, file="../data/genotypes/ukb_hm3_variant_ids.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)

fwrite(x=variant_names_rsid, file="../data/genotypes/ukb_hm3_variant_rsids.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)
