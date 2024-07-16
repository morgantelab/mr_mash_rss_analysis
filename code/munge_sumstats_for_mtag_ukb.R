###Load libraries needed
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--trait"), type="integer")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
output <- outparse$output
n <- outparse$n
trait <- outparse$trait

###Read sumstats
univ_sumstats <- readRDS(sumstats)

dat <- data.frame(snpid=rownames(univ_sumstats$alleles), chr=univ_sumstats$alleles$chr, bpos=univ_sumstats$alleles$bp,
                  a1=univ_sumstats$alleles$a1, a2=univ_sumstats$alleles$a2, freq=univ_sumstats$alleles$freq,
                  beta=univ_sumstats$Bhat[, trait], se=univ_sumstats$Shat[, trait], z=univ_sumstats$Bhat[, trait]/univ_sumstats$Shat[, trait],
                  pval=2*pt(abs(univ_sumstats$Bhat[, trait]/univ_sumstats$Shat[, trait]), n-2, lower.tail=FALSE),
                  n=n)
  
write.table(dat, file=output, sep=" ", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  