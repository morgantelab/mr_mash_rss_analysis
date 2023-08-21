###Load libraries needed
library(optparse)
library(bigsnpr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--genetic_map"), type="character", default=NULL)
parser <- add_option(parser, c("--samples_to_remove"), type="character", default=NULL)
parser <- add_option(parser, c("--chr"), type="integer")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--window_size"), type="integer", default=0)
parser <- add_option(parser, c("--thr_r2"), type="numeric", default=0)
parser <- add_option(parser, c("--impute_missing"), type="logical", default=FALSE)
parser <- add_option(parser, c("--sparse"), type="logical")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
genetic_map <- outparse$genetic_map
samples <- outparse$samples_to_remove
seed <- outparse$seed
chr <- outparse$chr
ncores <- outparse$ncores
window_size <- outparse$window_size
thr_r2 <- outparse$thr_r2
impute_missing <- outparse$impute_missing
sparse <- outparse$sparse
output <- outparse$output

###Set seed
set.seed(seed)

###
if(chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9")){
  chr <- paste0("0", chr)
}

###Read in data
geno <- snp_attach(geno_dat)

p <- geno$genotypes$ncol

n <- geno$genotypes$nrow

samples_to_remove <- tryCatch(data.table::fread(samples, data.table=FALSE, colClasses = "character", header=FALSE), 
                    error = function(e) {
                      return(NULL)
                    },
                    warning = function(w) {
                      return(NULL)
                    }
)

###Impute missing geno1types
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

###Set window size
if(window_size==0){
  win <- p
} else {
  win <- window_size
}

###Set chromosome
if(chr==0){
  col_inds <- seq_len(p)
} else {
  col_inds <- which(geno$map$chromosome==chr)
}

cat(sprintf("After optional filtering, %d variants remain.\n", length(col_inds)))

###Subset individuals (NB this only works with the UKB data)
if(!is.null(samples_to_remove)){
  if(ncol(samples_to_remove) > 1){
    samples_to_remove <- samples_to_remove[, 2]
  }
  
  ukb_data_location <- "/data2/morgante_lab/data/ukbiobank/genotypes/imputed/"
  ids <- bigreadr::fread2(paste0(ukb_data_location, "ukb22828_c1_b0_v3_s487271.sample"))
  ids <- ids[-1, ]
  
  ids1 <- data.table::fread("../data/misc/ukb_cleaned_bc_ind_ids.txt", data.table=FALSE, 
                            colClasses = "character", header=FALSE)
  
  #Keep only bc selected individuals
  ids_sel <- as.character(ids[as.character(ids$ID_2) %in% ids1[,2], 2])
  
  #Keep only training individuals  
  row_inds <- which(!(ids_sel %in% samples_to_remove))
} else {
  row_inds <- seq_len(n)
}

cat(sprintf("After optional filtering, %d individuals remain.\n", length(row_inds)))

###Compute LD
if(sparse){
  if(!is.null(genetic_map)){
    POS <- snp_asGeneticPos(as.integer(geno$map$chromosome), geno$map$physical.pos, 
                            dir=genetic_map, 
                            ncores=ncores)
  } else {
    POS <- geno$map$physical.pos
  }
  
  LD <- snp_cor(Gna=X, ind.col=col_inds, ind.row=row_inds, size=win, infos.pos=POS[col_inds], thr_r2=thr_r2, ncores=ncores)
} else {
  LD <- snp_cor(Gna=X, ind.col=col_inds, ind.row=row_inds, size=win, thr_r2=thr_r2, ncores=ncores)
}
  
###Save to binary file    
conn <- file(output, "wb")
writeBin(as.numeric(LD), conn)
close(conn)

