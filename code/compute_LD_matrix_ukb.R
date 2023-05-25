###Load libraries needed
library(optparse)
library(data.table)
library(bigsnpr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--inds_to_remove"), type="character")
parser <- add_option(parser, c("--vars_to_keep"), type="character")
parser <- add_option(parser, c("--genetic_map"), type="character", default=NULL)
parser <- add_option(parser, c("--chr"), type="integer")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--window_size"), type="integer", default=0)
parser <- add_option(parser, c("--thr_r2"), type="numeric", default=0)
parser <- add_option(parser, c("--impute_missing"), type="logical", default=FALSE)
parser <- add_option(parser, c("--sparse"), type="logical")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--temp_dir"), type="character")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
inds_to_remove <- outparse$inds_to_remove
vars_to_keep <- outparse$vars_to_keep
genetic_map <- outparse$genetic_map
seed <- outparse$seed
chr <- outparse$chr
ncores <- outparse$ncores
window_size <- outparse$window_size
thr_r2 <- outparse$thr_r2
impute_missing <- outparse$impute_missing
sparse <- outparse$sparse
output <- outparse$output
temp_dir <- outparse$temp_dir

###Set seed
set.seed(seed)

###Read in data
options(datatable.fread.datatable=FALSE)
geno_fam <- fread(paste0(unlist(strsplit(geno_dat, ".", fixed=TRUE))[1], ".fam"), showProgress=FALSE, header=FALSE)
inds_to_remove <- fread(inds_to_remove, header=FALSE, showProgress=FALSE)
inds_to_keep <- which(!(geno_fam[,2] %in% inds_to_remove[,2])) ##Get only individuals not used in simulations

tmp <- tempfile(tmpdir=temp_dir)
rds <- snp_readBed2(geno_dat, ind.row=inds_to_keep, backingfile=tmp, ncores=ncores)
geno <- snp_attach(rds)

vars_to_keep <- fread(vars_to_keep, header=FALSE, showProgress=FALSE)
vars_to_keep <- which(geno$map$marker.ID %in% vars_to_keep[,1]) ##Get only variants used in simulations

tmp1 <- tempfile(tmpdir=temp_dir)
rds1 <- snp_subset(geno, ind.col=vars_to_keep, backingfile=tmp1)
geno1 <- snp_attach(rds1)

rm(list=c("geno", "geno_fam", "inds_to_remove", "vars_to_keep"))

p <- geno1$genotypes$ncol

###Impute missing geno1types
if(impute_missing){
  X <- snp_fastImputeSimple(geno1$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno1$genotypes
}

###Set window size
if(window_size==0){
  win <- p
} else {
  win <- window_size
}

###Set chromosome
if(chr==0){
  inds <- seq_len(p)
} else {
  inds <- which(geno1$map$chromosome==chr)
}

###Compute LD
if(sparse){
  if(!is.null(genetic_map)){
    POS <- snp_asGeneticPos(geno1$map$chromosome, geno1$map$physical.pos, 
                            dir=genetic_map, 
                            ncores=ncores)
  } else {
    POS <- geno1$map$physical.pos
  }
  
  LD <- snp_cor(Gna=X, ind.col=inds, size=win, infos.pos=POS[inds], thr_r2=thr_r2, ncores=ncores)
} else {
  LD <- snp_cor(Gna=X, ind.col=inds, size=win, thr_r2=thr_r2, ncores=ncores)
}
  
###Save to binary file    
conn <- file(output, "wb")
writeBin(as.numeric(LD), conn)
close(conn)

###Remove temprary files
file.remove(paste0(tmp, c(".bk", ".rds")))
file.remove(paste0(tmp1, c(".bk", ".rds")))

