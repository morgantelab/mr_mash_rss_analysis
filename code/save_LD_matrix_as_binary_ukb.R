###Load libraries needed
library(optparse)
library(data.table)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
outparse <- parse_args(parser)

ncores <- outparse$ncores
input <- outparse$input
output <- outparse$output

###Load LD matrix
LD <- fread(input, data.table = FALSE, nThread=ncores)

###Save LD matrix as binary file
conn <- file(output, "wb")
writeBin(as.numeric(LD), conn)
close(conn)
