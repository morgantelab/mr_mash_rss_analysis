###Load libraries needed
library(optparse)
library(matrixStats)
library(bigsnpr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats_ols"), type="character")
parser <- add_option(parser, c("--sumstats_mixed"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--sparse"), type="logical", default=FALSE)
parser <- add_option(parser, c("--h2_init"), type="numeric", default=0.5)
parser <- add_option(parser, c("--shrink_corr"), type="numeric", default=1)
parser <- add_option(parser, c("--allow_jump_sign"), type="logical", default=TRUE)
parser <- add_option(parser, c("--verbose"), type="logical", default=FALSE)
parser <- add_option(parser, c("--burn_in"), type="integer", default=500)
parser <- add_option(parser, c("--num_iter"), type="integer", default=500)
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--trait"), type="integer")
parser <- add_option(parser, c("--temp_dir"), type="character")
outparse <- parse_args(parser)

sumstats_ols <- outparse$sumstats_ols
sumstats_mixed <- outparse$sumstats_mixed
LD_matrix <- outparse$LD_matrix
n <- outparse$n
sparse <- outparse$sparse
h2_init <- outparse$h2_init
shrink_corr <- outparse$shrink_corr
allow_jump_sign <- outparse$allow_jump_sign
verbose <- outparse$verbose
burn_in <- outparse$burn_in
num_iter <- outparse$num_iter
ncores <- outparse$ncores
output <- outparse$output
seed <- outparse$seed
trait <- outparse$trait
temp_dir <- outparse$temp_dir

vec_p_init <- seq_log(1e-4, 1, length.out = 30)

###Set seed
set.seed(seed)

###Read in data
univ_sumstats_ols <- readRDS(sumstats_ols)
univ_sumstats_mixed <- readRDS(sumstats_mixed)
p <- nrow(univ_sumstats_ols$Bhat)

if(tail(unlist(strsplit(LD_matrix, ".", fixed = TRUE)), 1) == "bin"){
  LD <- matrix(readBin(LD_matrix, what="numeric", n=p^2), nrow=p, ncol=p, byrow=TRUE)
  LD <- as(LD, "CsparseMatrix")
} else if(tail(unlist(strsplit(LD_matrix, ".", fixed = TRUE)), 1) == "rds"){
  LD <- readRDS(LD_matrix)
}

###Estimate N_eff
chisq_ols <- (univ_sumstats_ols$Bhat[, trait]/univ_sumstats_ols$Shat[, trait])^2
chisq_mixed <- (univ_sumstats_mixed$Bhat[, trait]/univ_sumstats_mixed$Shat[, trait])^2
to_keep <- which(chisq_ols>30)
chisq_ols <- chisq_ols[to_keep]
chisq_mixed <- chisq_mixed[to_keep]
sample_size_adj <- weightedMedian((chisq_mixed/chisq_ols), w=chisq_mixed)

###Prepare the sumstats and LD matrix
df_beta_ols <- data.frame(beta=univ_sumstats_ols$Bhat[, trait],
                          beta_se=univ_sumstats_ols$Shat[, trait],
                          n_eff=rep(n, times=p))

df_beta_mixed <- data.frame(beta=univ_sumstats_mixed$Bhat[, trait],
                            beta_se=univ_sumstats_mixed$Shat[, trait],
                            n_eff=rep(round(n*sample_size_adj), times=p))

tmp <- tempfile(tmpdir=temp_dir)
corr <- as_SFBM(LD, tmp, compact=TRUE)

###If initial estimate of h2 is not provided, compute it using LDSC
if(h2_init<0){
  ld <- Matrix::colSums(LD^2)
  (ldsc <- with(df_beta_ols, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff, blocks = NULL, 
                                  ncores=ncores)))
  h2_init <- ldsc[["h2"]]
  if(h2_init <= 0){
    h2_init <- 1e-4
  }
}

rm(LD)

###Fit mr.mash.rss
tic <- proc.time()[[3]]
fit_ldpred2_auto <- snp_ldpred2_auto(corr=corr, df_beta=df_beta_mixed, h2_init=h2_init, vec_p_init=vec_p_init,
                                     burn_in=burn_in, num_iter=num_iter, sparse=sparse, verbose=verbose,
                                     allow_jump_sign=allow_jump_sign, shrink_corr=shrink_corr, ncores=ncores)
toc <- proc.time()[[3]]

fit_ldpred2_auto$elapsed_time <- toc-tic

saveRDS(fit_ldpred2_auto, file=output)

###Remove temprary files
file.remove(paste0(tmp, ".sbk"))

