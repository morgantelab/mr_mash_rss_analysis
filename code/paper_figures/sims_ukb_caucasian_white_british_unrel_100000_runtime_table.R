set.seed(1)

###Set up
path <- "../output/"
prefix <- "ukb_caucasian_white_british_unrel_100000"
scenarioz <- "blocks_shared_effects_indep_resid"
univ_methodz <- c("ldpred2_auto")
mv_methodz <- c("mr_mash_rss", "mvbayesC")
chr <- 10
nreps <- 20
traits <- 5

n_col <- 4
n_row <- nreps * (length(univ_methodz)+length(univ_methodz)) * length(scenarioz)

res <- as.data.frame(matrix(NA, ncol=n_col, nrow=n_row))
colnames(res) <- c("rep", "scenario", "method", "time")


sce <- scenarioz[1]
repz <- 1:20

i <- 0

###Multivariate methods
for(met in mv_methodz){
  for(repp in repz){
    dat <- readRDS(paste0(path, met, "_fit", "/", prefix, "_", sce, "_chr",chr, "_", met, "_fit_", repp, ".rds"))
  
    i <- i + 1
      
    res[i, 1] <- repp
    res[i, 2] <- sce
    res[i, 3] <- met
    res[i, 4] <- dat$elapsed_time
  }
}

###Univariate methods
for(met in univ_methodz){
  for(repp in repz){
    elapsed_time <- vector("numeric", traits)
  
    for(trait in 1:traits){
      dat <- readRDS(paste0(path, met, "_fit", "/", prefix, "_", sce, "_chr",chr, "_", met, "_fit_trait", trait, "_", repp, ".rds"))
      elapsed_time[trait] <- dat$elapsed_time
    }
  
    elapsed_time <- sum(elapsed_time)
  
    i <- i + 1

    res[i, 1] <- repp
    res[i, 2] <- sce
    res[i, 3] <- met
    res[i, 4] <- elapsed_time
  }
}

print(do.call("rbind", tapply(res$time, res$method, summary)), digit=4)
