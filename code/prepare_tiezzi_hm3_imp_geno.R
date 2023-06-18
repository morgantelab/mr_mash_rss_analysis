library(bigsnpr)
library(bigreadr)
library(dplyr)
library(vctrs)
library(glue)
library(RhpcBLASctl)
library(parallel)

RhpcBLASctl::blas_set_num_threads(1)

###Load data
map_hapmap3 <- as.data.frame(readRDS("/data2/morgante_lab/ukbiobank_projects/mr_mash_rss/data/misc/map_hm3.rds"))
ids <- bigreadr::fread2("/data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c1_b0_v3_s487271.sample")
sel_ids <- bigreadr::fread2("/data2/morgante_lab/ukbiobank_projects/mr_mash_rss/data/misc/ukb_tiezzi_cleaned_BP_ind_ids.txt")

###Get list of HM3 SNPs to use
cl <- makeCluster(22)
clusterExport(cl, "map_hapmap3")
list_snp_id <- parLapply(1:22, function(chr) {
  mfi <- paste0("/data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb_mfi_chr", chr, "_v3.txt")
  infos_chr <- fread2(mfi, showProgress = FALSE)
  joined <- dplyr::inner_join(cbind(chr = chr, infos_chr), map_hapmap3[, c("chr", "rsid")],
                              by = c("chr" = "chr", "V2" = "rsid"))
  with(joined[!vctrs::vec_duplicate_detect(joined$V2), ],
       paste(chr, V3, V4, V5, sep = "_"))
})
stopCluster(cl)

###Keep only individuals retained by Francesco
ids <- ids[-1, ]
inds_to_keep <- which(ids[,2] %in% sel_ids[,2])

###Read in the genotype data and save the backing file
rds <- bigsnpr::snp_readBGEN(
  bgenfiles   = glue::glue("/data2/morgante_lab/data/ukbiobank/genotypes/imputed/ukb22828_c{chr}_b0_v3.bgen", chr = 1:22), 
  backingfile = "/data2/morgante_lab/ukbiobank_projects/mr_mash_rss/tmp/ukb_geno_imp_HM3_tiezzi",
  list_snp_id = list_snp_id,
  ind_row     = inds_to_keep,
  ncores      = 22
)  
  
rds$save()
