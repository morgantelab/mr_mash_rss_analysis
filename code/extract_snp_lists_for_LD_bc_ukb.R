###Load libraries needed
library(bigsnpr)
library(RhpcBLASctl)

###Set MKL threads
RhpcBLASctl::blas_set_num_threads(1)

###Set seed
set.seed(1)

###Load genotype data
geno_imp <- snp_attach("/scratch1/fabiom/ukb_bc_geno_imp_HM3.rds")

###Extract variant names
variant_names <- data.frame(geno_imp$map$marker.ID)
variant_names_rsid <- data.frame(geno_imp$map$rsid)

###Write out output
fwrite(x=variant_names, file="../data/genotypes/ukb_hm3_variant_ids.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)

fwrite(x=variant_names_rsid, file="../data/genotypes/ukb_hm3_variant_rsids.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)
