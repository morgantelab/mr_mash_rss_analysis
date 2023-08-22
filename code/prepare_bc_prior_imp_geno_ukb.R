library(data.table)
library(bigsnpr)
library(bigreadr)
library(dplyr)
library(vctrs)
library(glue)
library(RhpcBLASctl)
library(parallel)

RhpcBLASctl::blas_set_num_threads(1)

set.seed(1)

ukb_data_location <- "/data2/morgante_lab/data/ukbiobank/genotypes/imputed/"

###Load data
filenames <- list.files("../data/blood_cell_finemapping_regions", pattern="*.txt")

it <- 0
for(nam in filenames){
  it <- it+1
  
  dat <- fread(paste0("../data/blood_cell_finemapping_regions/", nam), showProgress=FALSE, 
               colClasses = "character", header=TRUE, data.table = FALSE)
  dat$chr <- readr::parse_number(unlist(strsplit(unlist(strsplit(nam, "_"))[2], ".", fixed=TRUE))[1])
  colnames(dat)[1] <- "rsid"
  if(it == 1){
    dat_full <- dat
  } else {
    dat_full <- bind_rows(dat_full, dat)
  }
}

ids <- bigreadr::fread2(paste0(ukb_data_location, "ukb22828_c1_b0_v3_s487271.sample"))
sel_ids <- bigreadr::fread2("../data/misc/ukb_cleaned_bc_ind_ids.txt")

###Get list of HM3 SNPs to use
cl <- makeCluster(22)
clusterExport(cl, c("dat_full", "ukb_data_location"))
list_snp_id <- parLapply(cl, 1:22, function(chr) {
  mfi <- paste0(ukb_data_location, "ukb_mfi_chr", chr, "_v3.txt")
  infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
  joined <- dplyr::inner_join(cbind(chr = chr, infos_chr), dat_full[, c("chr", "rsid")],
                              by = c("chr" = "chr", "V2" = "rsid"))
  with(joined[!vctrs::vec_duplicate_detect(joined$V2), ],
       paste(chr, V3, V4, V5, sep = "_"))
})
stopCluster(cl)

###Keep only selected individuals
ids <- ids[-1, ]
inds_to_keep <- which(as.character(ids[,2]) %in% as.character(sel_ids[,2]))

###Read in the genotype data and save the backing file
rds <- bigsnpr::snp_readBGEN(
  bgenfiles   = glue::glue(paste0(ukb_data_location, "ukb22828_c{chr}_b0_v3.bgen"), chr = 1:22), 
  backingfile = "/scratch1/fabiom/ukb_bc_geno_imp_prior",
  list_snp_id = list_snp_id,
  ind_row     = inds_to_keep,
  ncores      = 22
)  
