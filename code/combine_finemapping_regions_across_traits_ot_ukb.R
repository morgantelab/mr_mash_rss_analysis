###Load libraries needed
library(optparse)
library(data.table)
library(dplyr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--regions_prefix"), type="character")
parser <- add_option(parser, c("--regions_suffix"), type="character")
parser <- add_option(parser, c("--traits"), type="character")
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output_prefix_snp"), type="character")
parser <- add_option(parser, c("--output_boundaries"), type="character")
outparse <- parse_args(parser)

regions_prefix <- outparse$regions_prefix
regions_suffix <- outparse$regions_suffix
pheno_names <- eval(parse(text=outparse$traits))
output_prefix_snp <- outparse$output_prefix_snp
output_boundaries <- outparse$output_boundaries
ncores <- outparse$ncores

###Set seed
set.seed(1)

setDTthreads(ncores)


trait_regions = list()
for(trait in pheno_names){
  region <- fread(paste0(regions_prefix, trait, regions_suffix), data.table=FALSE, showProgress = FALSE)
  region <- region %>% arrange(desc(logp))
  region_r <- c()
  for(i in 1:22){
    region.chr <- region %>% filter(chr == i) %>% arrange(start)
    if(nrow(region.chr) == 0){
      next
    }
    tmp <- region.chr %>% group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>%
      summarise(start = first(start), end = max(end), logp = max(logp),.groups = 'drop') %>% 
      mutate(length = end - start) %>%
      mutate(chr = i) %>% select(chr, start, end, length, logp)
    region_r <- rbind(region_r, tmp)
  }
  trait_regions[[trait]] <- region_r
}

tb <- bind_rows(trait_regions, .id = "column_label")
res.final <- c()
for(i in 1:22){
  tb.chr <- tb %>% filter(chr == i) %>% arrange(start)
  if(nrow(tb.chr) == 0){
    next
  }
  tmp <- tb.chr %>% group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>%
    summarise(start = first(start), end = max(end), logp = max(logp), .groups = 'drop') %>% 
    mutate(length = end - start) %>%
    mutate(chr = i) %>% select(chr, start, end, length, logp)
  res.final <- rbind(res.final, tmp)
}

gwas <- fread(paste0("../output/GWAS_for_regions/ukb_ot_chrAll_biallelic_snps_only_recoded_ids_info_06_gwas_for_regions_rsid_", 
                     pheno_names[1], "_plink_mtag_combined.txt"), data.table=FALSE, showProgress = FALSE)

snpsnum <- c()
for(i in 1:nrow(res.final)){
  gwas_reg <- gwas %>% filter(chr == res.final$chr[i],
                              bp >= res.final$start[i], 
                              bp <= res.final$end[i])
  
  df <- data.frame(ID=gwas_reg$rsid)

  write.table(df, file=paste0(output_prefix_snp, "_chr", res.final$chr[i],
                              ".", res.final$start[i], ".", res.final$end[i], "_snp_names.txt"),
              col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
  
  snpsnum <- c(snpsnum, nrow(gwas_reg))
}

res.final$snpsnum <- snpsnum

res.final$logp <- NULL

write.table(res.final, file=output_boundaries, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")