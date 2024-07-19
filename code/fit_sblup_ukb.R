###Load libraries needed
library(optparse)
library(qgg)
library(SumTool)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_scores"), type="character")
parser <- add_option(parser, c("--geno"), type="character")
parser <- add_option(parser, c("--impute_missing"), type="logical")
parser <- add_option(parser, c("--verbose"), type="logical", default=FALSE)
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--trait"), type="integer")
parser <- add_option(parser, c("--wind_size"), type="numeric", default=2000000)
parser <- add_option(parser, c("--temp_dir"), type="character")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
LD_scores <- outparse$LD_scores
geno <- outparse$geno
impute_missing <- outparse$impute_missing
verbose <- outparse$verbose
output <- outparse$output
n <- outparse$n
ncores <- outparse$ncores
trait <- outparse$trait
seed <- outparse$seed
wind_size <- outparse$wind_size
temp_dir <- outparse$temp_dir

###Set seed
set.seed(seed)

###Read in data
univ_sumstats <- readRDS(sumstats)
ld_scores <- read.table(gzfile(LD_scores), header=TRUE, sep="\t")
tmp <- tempfile(tmpdir=temp_dir)
geno_dat <- SumTool::read_binary(bfile=geno, impute=impute_missing, additive=TRUE, out=tmp, verbose=verbose, threads=ncores)
geno_bed <- geno_dat$geno
geno_map <- geno_dat$map
rm(geno_dat)

###Compute Z scores
Z <- univ_sumstats$Bhat / univ_sumstats$Shat

###Compute genetic parameters
ldscores <- ld_scores$L2
names(ldscores) <- ld_scores$SNP
h2 <- qgg::ldsc(ldscores=ldscores, z=Z, n=rep(n, ncol(Z)), what="h2")

###Organize sumstats for SBLUP
sumstats <- data.frame(SNP=rownames(univ_sumstats$alleles), CHR=univ_sumstats$alleles$chr, BP=univ_sumstats$alleles$bp,
                        A1=univ_sumstats$alleles$a1, A2=univ_sumstats$alleles$a2, BETA=univ_sumstats$Bhat[, trait], 
                        SE=univ_sumstats$Shat[, trait], NMISS=n)

###Fit SBLUP
lambda <- nrow(sumstats)*(1/h2[trait]-1)

tic <- proc.time()[[3]]

fit_sblup <- SumTool::SBLUP(sumstat=sumstats, geno=geno_bed, map=geno_map, lambda=lambda, w=wind_size, threads=ncores, verbose=verbose)

toc <- proc.time()[[3]]

fit_sblup$elapsed_time <- toc-tic

###Save output
saveRDS(fit_sblup, file=output)

###Remove temprary files
file.remove(paste0(tmp, c(".bin", ".map", ".desc")))
