###Load libraries needed
library(optparse)

###Function to perform RINT
inv_normalise <- function(x) { #this would also tolerate NAs
  return( qnorm( (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--pheno"), type="character")
outparse <- parse_args(parser)

output <- outparse$output
pheno <- outparse$pheno

###Set seed
set.seed(1)

###Read in data
pheno_dat <- readRDS(pheno)

pheno_names <- c("DPa", "SPa", "height", "hip", "BMI", "BMR", "WHR")

res_mat <- matrix(as.numeric(NA), nrow=nrow(pheno_dat), ncol=length(pheno_names))
colnames(res_mat) <- pheno_names
rownames(res_mat) <- pheno_dat$id

for(i in 1:length(pheno_names)){
  trait <- pheno_names[i]
  dat <- data.frame(y=pheno_dat[, trait], pheno_dat[, c("sex", "assessment_centre", "age", "genotype_measurement_batch",
                                                        paste0("pc_genetic", 1:10))])
  fit_norm <- lm(y ~ as.factor(sex) + as.factor(assessment_centre) + age + I(age^2) + as.factor(genotype_measurement_batch) +
                     pc_genetic1 + pc_genetic2 + pc_genetic3 + pc_genetic4 + pc_genetic5 + pc_genetic6 + pc_genetic7 + pc_genetic8 +
                     pc_genetic9 + pc_genetic10, dat)
  res <- resid(fit_norm) + coef(fit_norm)[1]
  res_mat[, i] <- inv_normalise(res)
}

dat_final <- data.frame(FID=pheno_dat$id, IID=pheno_dat$id, res_mat)

write.table(dat_final, file=output, col.names = TRUE, row.names=FALSE, quote=FALSE, sep="\t")
