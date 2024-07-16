###Load libraries needed
library(optparse)
library(Matrix)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--sumstats"), type="character")
parser <- add_option(parser, c("--LD_matrix"), type="character")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

sumstats <- outparse$sumstats
LD_matrix <- outparse$LD_matrix
output <- outparse$output

###Read in data
univ_sumstats <- readRDS(sumstats)
p <- nrow(univ_sumstats$Bhat)
r <- ncol(univ_sumstats$Bhat)
LD <- matrix(readBin(LD_matrix, what="numeric", n=p^2), nrow=p, ncol=p, byrow=TRUE)
LD <- as(LD, "CsparseMatrix")

saveRDS(LD, file=output)