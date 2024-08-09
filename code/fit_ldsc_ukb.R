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
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
parser <- add_option(parser, c("--trait"), type="integer")
outparse <- parse_args(parser)

sumstats_prefix <- outparse$sumstats_prefix
sumstats_suffix <- outparse$sumstats_suffix
LD_matrix_prefix <- outparse$LD_matrix_prefix
LD_matrix_suffix <- outparse$LD_matrix_suffix
chr <- eval(parse(text=outparse$chr))
n <- outparse$n
ncores <- outparse$ncores
output <- outparse$output
seed <- outparse$seed
trait <- outparse$trait

###Set seed
set.seed(seed)

###Loop over chromosomes
it <- 0

for(i in chr){
  it <- it+1
  
  ##Read in summary stats and LD matrix
  dat_ss <- readRDS(paste0(sumstats_prefix, i, sumstats_suffix))
  
  p <- nrow(dat_ss$Bhat)
  LD <- readRDS(paste0(LD_matrix_prefix, i, LD_matrix_suffix))

  ##Merge summary stats and LD scores
  if(it==1){
    Bhat_all <- dat_ss$Bhat
    Shat_all <- dat_ss$Shat
    
    ld <- Matrix::colSums(LD^2)
    
  } else {
    Bhat_all <- rbind(Bhat_all, dat_ss$Bhat)
    Shat_all <- rbind(Shat_all, dat_ss$Shat)
    
    ld <- c(ld, Matrix::colSums(LD^2))
  }
}

rm(list=c("LD", "dat_ss"))
gc()

###Prepare the sumstats
p <- nrow(Bhat_all)
df_beta <- data.frame(beta=Bhat_all[, trait],
                      beta_se=Shat_all[, trait],
                      n_eff=rep(n, times=p))

###Run LDSC
fit_ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                   sample_size = n_eff, ncores = ncores))

saveRDS(fit_ldsc, file=output)
