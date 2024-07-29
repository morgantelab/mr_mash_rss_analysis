###Load libraries needed
library(optparse)
library(qgg)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--method"), type="character")
parser <- add_option(parser, c("--h2"), type="numeric", default=NULL)
parser <- add_option(parser, c("--pi"), type="numeric", default=0.0001)
parser <- add_option(parser, c("--nit"), type="integer", default=5000)
parser <- add_option(parser, c("--nburn"), type="integer", default=1000)
parser <- add_option(parser, c("--nthin"), type="integer", default=5)
parser <- add_option(parser, c("--trait"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
LD_matrix <- outparse$LD_matrix
method <- outparse$method
h2 <- outparse$h2
Pi <- outparse$pi
nit <- outparse$nit
nburn <- outparse$nburn
nthin <- outparse$nthin
trait <- outparse$trait
output <- outparse$output
seed <- outparse$seed

###Set seed
set.seed(seed)

###Read in data
univ_sumstats <- readRDS(sumstats)
LD <- qgg:::readLD(LD_matrix)

###Organize summary stats in the format required by qgg
stat <- data.frame(b=univ_sumstats$Bhat[, trait], seb=univ_sumstats$Shat[, trait], 
                   n=univ_sumstats$n[trait])

###Fit sBayesX
tic <- proc.time()[[3]]
fit_sbayes <- qgg:::sblr(stat=stat, LD=LD, method=method, pi=Pi, h2=h2,
                         nit=nit, nburn=nburn, nthin=nthin) 

toc <- proc.time()[[3]]

fit_sbayes$elapsed_time <- toc-tic

saveRDS(fit_sbayes, file=output)


