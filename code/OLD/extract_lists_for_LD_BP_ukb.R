###Load libraries needed
library(bigsnpr)
library(data.table)
library(RhpcBLASctl)

###Set MKL threads
RhpcBLASctl::blas_set_num_threads(1)

###Set seed
set.seed(1)

###Load genotype data
geno_imp <- snp_attach("/scratch1/fabiom/ukb_geno_imp_HM3_tiezzi.rds")
geno_array_fam <- fread("/data2/morgante_lab/data/ukbiobank/genotypes/array/ukb22418_all_auto_b0_v2_s488243_caucasian_white_british_unrel.fam", 
                        showProgress=FALSE, header=FALSE, data.table = FALSE)[, 1:2]
pheno_dat <- readRDS("../data/phenotypes/ukb_tiezzi_cleaned_BP_covar_pheno.rds")

###Extract variant names
variant_names <- data.frame(geno_imp$map$marker.ID)
variant_names_rsid <- data.frame(geno_imp$map$rsid)

###Get IDs of individuals not used for summary stats computation
geno_array_fam_filt <- geno_array_fam[-which(as.character(geno_array_fam[, 2]) %in% as.character(pheno_dat$ID)), ]

###Write out output
fwrite(x=variant_names, file="../data/genotypes/ukb_hm3_variant_ids.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)

fwrite(x=variant_names_rsid, file="../data/genotypes/ukb_hm3_variant_rsids.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)

fwrite(x=geno_array_fam_filt, file="../data/LD_matrices/ukb_tiezzi_BP_inds_for_LD.txt", 
       sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)
