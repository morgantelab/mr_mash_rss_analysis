###Load libraries needed
library(optparse)
library(data.table)
library(dplyr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--region_size"), type="integer")
parser <- add_option(parser, c("--mhc"), type="logical")
parser <- add_option(parser, c("--sig_threshold"), type="numeric")
parser <- add_option(parser, c("--gwas_method"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
region_size <- outparse$region_size
mhc <- outparse$mhc
sig_threshold <- outparse$sig_threshold
gwas_method <- outparse$gwas_method
output <- outparse$output
ncores <- outparse$ncores

###Set seed
set.seed(1)

setDTthreads(ncores)

dat <- fread(sumstats, data.table=FALSE, showProgress = FALSE)
  
if(!mhc){
  dat.na <- dat %>% filter(! (chr == 6 & bp>=25000000 & bp<=36000000))
} else {
  dat.na <- dat
}
  
if(gwas_method=="mtag"){
  dat.na$logp <- -log10(dat.na$mtag_p)
} else if(gwas_method=="plink"){
  dat.na$logp <- -log10(dat.na$plink_p)
}
  
chrbp <- dat.na %>% group_by(chr) %>% summarise(start = min(bp), end = max(bp))

res = matrix(NA, 0, 5)
  
i = 0
while(max(dat.na$logp, na.rm=T) > -log10(sig_threshold)){
  i <- i + 1
  signal <- dat.na %>% filter(logp == max(logp, na.rm=T)) %>% top_n(1, bp) %>% select(chr, bp, logp) %>%
    mutate(start = max(bp - region_size, chrbp$start[which(chrbp$chr == chr)]),
           end = min(bp + region_size, chrbp$end[which(chrbp$chr == chr)]))
  res <- rbind(res, signal)
  dat.na[which(dat.na$chr == signal$chr & dat.na$bp <= signal$end & dat.na$bp >= signal$start),] <- NA
}
  
fwrite(res, output, quote = FALSE, col.names=TRUE, row.names = FALSE, sep = "\t")





