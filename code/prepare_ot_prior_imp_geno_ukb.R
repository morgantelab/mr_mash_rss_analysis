###Load libraries needed
library(optparse)
library(bigsnpr)
library(bigreadr)
library(dplyr)
library(vctrs)
library(glue)
library(parallel)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--regions_location"), type="character")
parser <- add_option(parser, c("--ukb_geno_location"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--selected_inds"), type="character")
outparse <- parse_args(parser)

output <- outparse$output
ncores <- outparse$ncores
regions_location <- outparse$regions_location
ukb_geno_location <- outparse$ukb_geno_location
selected_inds <- outparse$selected_inds


set.seed(1)

ukb_data_location <- "/data2/morgante_lab/data/ukbiobank/genotypes/imputed/"

###Load data
filenames <- list.files(regions_location, pattern="*.txt")

it <- 0
for(nam in filenames){
  it <- it+1
  
  dat <- fread(paste0(regions_location, "/", nam), showProgress=FALSE, 
               colClasses = "character", header=TRUE, data.table = FALSE)
  dat$chr <- readr::parse_number(unlist(strsplit(unlist(strsplit(nam, "_"))[2], ".", fixed=TRUE))[1])
  colnames(dat)[1] <- "rsid"
  if(it == 1){
    dat_full <- dat
  } else {
    dat_full <- bind_rows(dat_full, dat)
  }
}

ids <- bigreadr::fread2(paste0(ukb_geno_location, "ukb22828_c1_b0_v3_s487271.sample"))
sel_ids <- bigreadr::fread2(selected_inds)

###Get list of HM3 SNPs to use
cl <- makeCluster(22)
clusterExport(cl, c("dat_full", "ukb_geno_location"))
list_snp_id <- parLapply(cl, 1:22, function(chr) {
  mfi <- paste0(ukb_geno_location, "ukb_mfi_chr", chr, "_v3.txt")
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
  bgenfiles   = glue::glue(paste0(ukb_geno_location, "ukb22828_c{chr}_b0_v3.bgen"), chr = 1:22), 
  backingfile = output,
  list_snp_id = list_snp_id,
  ind_row     = inds_to_keep,
  ncores      = 22
)  
