###Load libraries needed
library(optparse)

###Function needed
is_strong <- function(x){
  strong <- any(abs(x) > 2)
  return(strong)
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input"), type="character")
parser <- add_option(parser, c("--data_id"), type="integer")
outparse <- parse_args(parser)

input <- outparse$input
data_id <- outparse$data_id

###Extract weak effects
for(i in 1:22){
  dat <- readRDS(paste0("../output/summary_statistics/", input,"_chr", i, "_sumstats_", data_id, ".rds"))
  Z <- dat$Bhat/dat$Shat
  strong <- which(apply(Z, 1, is_strong))
  
  if(i==1){
    Z_weak <- Z[-strong, ]
  } else {
    Z_weak <- rbind(Z_weak, Z[-strong, ])
  }
  
}

###Estimate V
Vhat <- matrix(0, nrow=ncol(Z_weak), ncol=ncol(Z_weak))

for(j in 1:nrow(Z_weak)){
  Vhat <- Vhat + tcrossprod(Z_weak[j,])
}

Vhat <- Vhat/nrow(Z_weak)

###Save file
saveRDS(Vhat, file=paste0("../output/misc/", input, "_residual_cov_", data_id, ".rds"))
