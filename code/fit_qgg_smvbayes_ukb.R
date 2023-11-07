###Load libraries needed
library(optparse)
library(qgg)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--method"), type="character")
parser <- add_option(parser, c("--models"), type="character", default = NULL)
parser <- add_option(parser, c("--h2"), type="numeric", default=0.5)
parser <- add_option(parser, c("--pi"), type="numeric", default=0.0001)
parser <- add_option(parser, c("--nit"), type="integer", default=5000)
parser <- add_option(parser, c("--nburn"), type="integer", default=1000)
parser <- add_option(parser, c("--nthin"), type="integer", default=5)
parser <- add_option(parser, c("--traits"), type="character", default="-1")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
LD_matrix <- outparse$LD_matrix
n <- outparse$n
method <- outparse$method
models <- outparse$models
h2 <- outparse$h2
Pi <- outparse$pi
nit <- outparse$nit
nburn <- outparse$nburn
nthin <- outparse$nthin
traits <- eval(parse(text=outparse$traits))
output <- outparse$output
seed <- outparse$seed

###Set seed
set.seed(seed)

###Read in data
univ_sumstats <- readRDS(sumstats)
LD <- qgg:::readLD(LD_matrix)

if(length(traits)>1){
  univ_sumstats <- lapply(univ_sumstats, function(x, sel){x[, sel]}, traits)
}

if(models == "NULL"){
  models <- NULL
}

###Organize summary stats in the format required by qgg
names(univ_sumstats)[names(univ_sumstats)=="Bhat"] <- "b" 
names(univ_sumstats)[names(univ_sumstats)=="Shat"] <- "seb" 
univ_sumstats$n <- matrix(n,ncol=ncol(univ_sumstats$b),nrow=nrow(univ_sumstats$b))

###Fit smvBayesX
tic <- proc.time()[[3]]
fit_smvbayes <- qgg:::mtsblr(stat=univ_sumstats, LD=LD, method=method, pi=Pi, h2=h2,
                             nit=nit, nburn=nburn, nthin=nthin, models=models) 

toc <- proc.time()[[3]]

fit_smvbayes$elapsed_time <- toc-tic

saveRDS(fit_smvbayes, file=output)


