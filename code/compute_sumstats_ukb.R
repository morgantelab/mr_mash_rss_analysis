###Load libraries needed
library(optparse)
library(data.table)
library(bigsnpr)

###Function to compute summary stats
compute_univariate_sumstats_bigsnp <- function(geno_obj, Y, chr=0, standardize=FALSE, 
                                               impute_missing=TRUE,
                                               standardize.response=FALSE, mc.cores=1){
  
  if(impute_missing){
    X <- bigsnpr::snp_fastImputeSimple(geno_obj$genotypes, method = "mean2", ncores = ncores)
  } else {
    X <- geno_obj$genotypes
  }

  r <- ncol(Y)
  p <- geno_obj$genotypes$ncol
  
  variable_names <- geno_obj$map$marker.ID
  response_names <- colnames(Y)
  
  if(standardize.response){
    Y <- mr.mash.alpa::scale_fast2(Y, scale=TRUE)$M
  }
  
  
  if(chr==0){
    inds <- seq_len(p)
  } else {
    inds <- which(geno_obj$map$chromosome==chr)
  }
  
  linreg <- function(i, X, Y, cols){
    
    fit <- bigstatsr::big_univLinReg(X=X, y=Y[, i], ind.col=cols, ncores = 1)
    bhat <- fit$estim
    shat <- fit$std.err
    
    return(list(bhat=bhat, shat=shat))
  }
  
  if(mc.cores>1){
    cl <- parallel::makeCluster(mc.cores)
    out <- parallel::parLapply(cl, 1:r, linreg, X, Y, inds)
    parallel::stopCluster(cl)
  } else {
    out <- lapply(1:r, linreg, X, Y, inds)
  }
  
  Bhat <- sapply(out,"[[","bhat")
  Shat <- sapply(out,"[[","shat")
  colnames(Bhat) <- colnames(Shat) <- response_names
  rownames(Bhat) <- rownames(Shat) <- variable_names[inds]
  
  ###Put coefficients and ses on the standardized X scale   
  if(standardize){
    sds <- sqrt(bigstatsr::big_colstats(X, ind.col=inds)$var)
  	
    Bhat <- Bhat*sds
    Shat <- Shat*sds
  }
  
  return(list(Bhat=Bhat, Shat=Shat))
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--pheno"), type="character")
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--test_ids"), type="character")
parser <- add_option(parser, c("--outprefix"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
parser <- add_option(parser, c("--standardize"), type="logical")
parser <- add_option(parser, c("--data_id"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--impute_missing"), type="logical")
outparse <- parse_args(parser)

pheno_dat <- outparse$pheno
geno_dat <- outparse$geno
test_ids <- outparse$test_ids
outprefix <- outparse$outprefix
data_id <- outparse$data_id
chr <- outparse$chr
standardize <- outparse$standardize
ncores <- outparse$ncores
impute_missing <- outparse$impute_missing

###Set seed
set.seed(data_id)

###Read in data and filter them
options(datatable.fread.datatable=FALSE)
geno_fam <- fread(paste0("../data/genotypes/", geno_dat, ".fam"), showProgress=FALSE)
test_ids <- fread(paste0("../data/phenotypes/", test_ids, "_", data_id, ".txt"), showProgress=FALSE)
training_inds_geno <- which(!(geno_fam[,2] %in% test_ids[,2])) ##Get only training individuals

pheno <- readRDS(paste0("../data/phenotypes/", pheno_dat, "_", data_id, ".rds"))$Y
training_inds_pheno <- which(!(rownames(pheno) %in% test_ids[,2])) ##Get only training individuals
pheno <- pheno[training_inds_pheno, ]

tmp <- tempfile(tmpdir="/data2/morgante_lab/fabiom/tmp")
rds <- snp_readBed2(paste0("../data/genotypes/", geno_dat, ".bed"), ind.row=training_inds_geno, backingfile=tmp, ncores=ncores)
geno <- snp_attach(rds)

rm(list=c("geno_fam", "test_ids", "training_inds_geno", "training_inds_pheno"))

###Compute summary stats
chrscar <- unlist(strsplit(chr, ":"))

if(length(chrscar)==1){
  chrs <- as.integer(chrscar)
} else {
  chrs <- as.integer(chrscar[1]):as.integer(chrscar[2])
}

for(i in chrs){ 
  out <- compute_univariate_sumstats_bigsnp(geno_obj=geno, Y=pheno, chr=i, standardize=standardize,
                                            impute_missing=impute_missing,
                                            standardize.response=FALSE, mc.cores=ncores)

  if(length(chrs)==1 && chrs==0){
    i <- "All"
  }
  saveRDS(out, file=paste0("../output/summary_statistics/", outprefix,"_chr", i, "_sumstats_", data_id, ".rds"))
}

###Remove temprary files
file.remove(paste0(tmp, c(".bk", ".rds")))
