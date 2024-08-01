###Load libraries needed
library(optparse)
library(bigsnpr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--genetic_map"), type="character", default=NULL)
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--window_size"), type="integer", default=0)
parser <- add_option(parser, c("--thr_r2"), type="numeric", default=0)
parser <- add_option(parser, c("--thr_maf"), type="numeric", default=0)
parser <- add_option(parser, c("--impute_missing"), type="logical", default=FALSE)
parser <- add_option(parser, c("--sparse"), type="logical")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--temp_dir"), type="character")
outparse <- parse_args(parser)

geno_dat <- outparse$geno
genetic_map <- outparse$genetic_map
seed <- outparse$seed
ncores <- outparse$ncores
window_size <- outparse$window_size
thr_r2 <- outparse$thr_r2
thr_maf <- outparse$thr_maf
impute_missing <- outparse$impute_missing
sparse <- outparse$sparse
output <- outparse$output
temp_dir <- outparse$temp_dir

###Set seed
set.seed(seed)

###Read in data
tmp <- tempfile(tmpdir=temp_dir)
rds <- snp_readBed2(geno_dat, backingfile=tmp, ncores=ncores)
geno <- snp_attach(rds)

p <- geno$genotypes$ncol

###Impute missing genotypes
if(impute_missing){
  X <- snp_fastImputeSimple(geno$genotypes, method = "mean2", ncores = ncores)
} else {
  X <- geno$genotypes
}

###Compute maf and filter out variants by it
if(thr_maf>0){
  maf <- snp_MAF(G=X, ncores=ncores)
  inds <- which(maf>=thr_maf)
} else {
  inds <- seq_len(p)
}

###Set window size
if(window_size==0){
  win <- p
} else {
  win <- window_size
}

###Compute LD
if(sparse){
  if(!is.null(genetic_map)){
    POS <- snp_asGeneticPos(geno$map$chromosome, geno$map$physical.pos, 
                            dir=genetic_map, 
                            ncores=ncores)
  } else {
    POS <- geno$map$physical.pos
  }
  
  LD <- snp_cor(Gna=X, ind.col=inds, size=win, infos.pos=POS[inds], thr_r2=thr_r2, ncores=ncores)
} else {
  LD <- snp_cor(Gna=X, ind.col=inds, size=win, thr_r2=thr_r2, ncores=ncores)
}

###Filter map, add allele frequency and save it to a file
map_filtered <- geno$map[inds, ]
statz <- bigstatsr::big_colstats(X=X, ind.col=inds, ncores=ncores)
map_filtered$freq <- (statz$sum/nrow(X))/2
write.table(map_filtered, paste0(output, ".map"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  
###Save to binary file    
conn <- file(paste0(output, ".ld.bin"), "wb")
writeBin(as.numeric(LD), conn)
close(conn)

###Save to rds
saveRDS(LD, paste0(output, ".rds"))

###Remove temporary files
file.remove(paste0(tmp, c(".bk", ".rds")))

