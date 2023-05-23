###Load libraries needed
library(optparse)
library(bigsnpr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats_prefix"), type="character")
parser <- add_option(parser, c("--sumstats_suffix"), type="character")
parser <- add_option(parser, c("--LD_matrix_prefix"), type="character")
parser <- add_option(parser, c("--LD_matrix_suffix"), type="character")
parser <- add_option(parser, c("--chr"), type="character")
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
outparse <- parse_args(parser)

sumstats_prefix <- outparse$sumstats_prefix
sumstats_suffix <- outparse$sumstats_suffix
LD_matrix_prefix <- outparse$LD_matrix_prefix
LD_matrix_suffix <- outparse$LD_matrix_suffix
chr <- outparse$chr
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

vec_p_init <- seq_log(1e-4, 0.2, length.out = 30)

###Set seed
set.seed(seed)

###Read in data
chrscar <- unlist(strsplit(chr, ":"))

if(length(chrscar)==1){
  chrs <- as.integer(chrscar)
} else {
  chrs <- as.integer(chrscar[1]):as.integer(chrscar[2])
}

###Loop through chromosomes
it <- 0
tmp <- tempfile(tmpdir="/data2/morgante_lab/fabiom/tmp")

for(i in chrs){
  it <- it+1
  
  ##Read in summary stats and LD matrix
  dat_ss <- readRDS(paste0(sumstats_prefix, i, sumstats_suffix))
  
  p <- nrow(dat_ss$Bhat)
  LD <- matrix(readBin(paste0(LD_matrix_prefix, i, LD_matrix_suffix), what="numeric", n=p^2), 
               nrow=p, ncol=p, byrow=TRUE)
  
  ##Merge summary stats and LD matrix
  if(it==1){
    Bhat_all <- dat_ss$Bhat
    Shat_all <- dat_ss$Shat
    
    corr <- as_SFBM(as(LD, "CsparseMatrix"), tmp, compact = TRUE)
    
  } else {
    Bhat_all <- rbind(Bhat_all, dat_ss$Bhat)
    Shat_all <- rbind(Shat_all, dat_ss$Shat)
    
    corr$add_columns(as(LD, "CsparseMatrix"), nrow(corr))
  }
}

###Prepare the sumstats and LD matrix
p <- nrow(Bhat_all)
df_beta <- data.frame(beta=Bhat_all[, trait],
                      beta_se=Shat_all[, trait],
                      n_eff=rep(n, times=p))

list=c("LD", "dat_ss", "Bhat_all", "Shat_all")

###Fit mr.mash.rss
fit_ldpred2_auto <- snp_ldpred2_auto(corr=corr, df_beta=df_beta, h2_init=h2_init, vec_p_init=vec_p_init,
                                     burn_in=burn_in, num_iter=num_iter, sparse=sparse, verbose=verbose,
                                     allow_jump_sign=allow_jump_sign, shrink_corr=shrink_corr, ncores=ncores)

saveRDS(fit_ldpred2_auto, file=output)

###Remove temprary files
file.remove(paste0(tmp, ".sbk"))

