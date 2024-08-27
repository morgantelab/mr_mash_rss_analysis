###Load libraries needed
library(optparse)
library(data.table)
library(dplyr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--raw_pheno"), type="character")
parser <- add_option(parser, c("--bc_pheno"), type="character")
parser <- add_option(parser, c("--n_fold"), type="integer", default=-1)
parser <- add_option(parser, c("--ncores"), type="integer")
parser <- add_option(parser, c("--output_inds_list"), type="character")
parser <- add_option(parser, c("--output_pheno_data"), type="character")
outparse <- parse_args(parser)

ncores <- outparse$ncores
raw_pheno <- outparse$raw_pheno
bc_pheno <- outparse$bc_pheno
n_fold <- outparse$n_fold
output_pheno_data <- outparse$output_pheno_data
output_inds_list <- outparse$output_inds_list

###Set seed
set.seed(1)

###Set data.table threads
setDTthreads(ncores)

###Load data
##Raw phenotypes
dat <- fread(raw_pheno, data.table = FALSE, showProgress=FALSE, colClasses = "character")

##Cleaned blood cell traits
dat_bc <- readRDS(bc_pheno)
dat_covar <- dat_bc[, c("id", "sex", "assessment_centre", "age", "genotype_measurement_batch", paste0("pc_genetic", 1:10), "fold")]
rm(dat_bc); gc()

###Select columns to keep and filter the raw data
#weight, waist, hip, BMI, TFM, BFP
cols_to_keep <- c("eid", "21002-0.0", "48-0.0", "49-0.0", "21001-0.0", "23128-0.0", "23099-0.0")

dat_filt <- dat[, cols_to_keep]
colnames(dat_filt) <- c("id", "weight", "waist", "hip", "BMI", "TFM", "BFP")

rm(dat); gc()

###Convert numerical columns except the first one (the sample ids) to numeric values, and set all empty strings to NA.
n <- length(dat_filt)
for (i in 2:n) {
  x          <- dat_filt[,i]
  x[x == ""] <- as.character(NA)
  if(grepl('med', colnames(dat_filt)[i])){
    dat_filt[,i]    <- x
  }else{
    dat_filt[,i]    <- as.numeric(x)
  }
}

###Join bc covariate data with BP data
dat_filt1 <- inner_join(dat_covar, dat_filt, by = join_by(id))
rm(dat_filt); gc()

# library(reshape2)
# library(ggplot2)
# 
# corr <- cor(dat_filt1[, c("weight", "waist", "hip", "BMI", "TFM", "BFP")], use="pairwise.complete.obs")
# 
# melted_cormat <- melt(corr)
# 
# ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                        midpoint = 0, limit = c(-1,1), space = "Lab",
#                        name="Pearson\nCorrelation") +
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1,
#                                    size = 12, hjust = 1))+
#   coord_fixed()

###Select phenotypes for the analysis and keep only complete cases
dat_filt1 <- dat_filt1[complete.cases(dat_filt1), ]

###Remove individuals with "abnormal" measurements
pheno_names <- c("weight", "waist", "hip", "BMI", "TFM", "BFP")
pheno_names_rint <- paste(pheno_names, "rint", sep="_")

##RINT
dat_filt2 <- matrix(as.numeric(NA), nrow=nrow(dat_filt1), ncol=length(pheno_names_rint))
colnames(dat_filt2) <- pheno_names_rint
dat_filt_final <- cbind(dat_filt1, dat_filt2)

for(name in pheno_names){
  id_rint <- which(colnames(dat_filt_final) == paste(name, "rint", sep="_"))
  id <- which(colnames(dat_filt_final) == name)
  dat_filt_final[, id_rint] = qnorm((rank(dat_filt_final[,id],na.last="keep")-0.5)/sum(!is.na(dat_filt_final[,id])))
}

##Compute empirical covariance matrix
covy <- dat_filt_final %>% select(all_of(pheno_names_rint)) %>% cov
D2 <- stats::mahalanobis(dat_filt_final %>% select(all_of(pheno_names_rint)), center=0, cov=covy) ## mahalanobis distance
dat_filt_final <- dat_filt_final[D2 < qchisq(0.01, df=length(pheno_names), lower.tail = FALSE),]

###Drop unneeded columns
dat_filt_final <- dat_filt_final[, -which(colnames(dat_filt_final) %in% pheno_names_rint)]

###Sample within fold, if requested
if(n_fold > 0){
  ##Write individuals id BEFORE sampling
  indiv_list <- cbind(dat_filt_final$id, dat_filt_final$id)
  
  base_out <- unlist(strsplit(output_inds_list, ".", fixed=TRUE))
  base_out_full <- paste(base_out[1], "full_data", sep="_")
  output_inds_list_full <- paste(base_out_full, base_out[2], sep=".")
  
  fwrite(x=indiv_list, file=output_inds_list_full, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)
  
  dat_filt_final <- dat_filt_final %>% group_by(fold) %>% slice_sample(n=n_fold) %>% ungroup() %>% 
    arrange(as.numeric(id)) %>% as.data.frame()
}

cat("The final data set has", nrow(dat_filt_final), "individuals.\n")

cat("The training set corresponding to each test set has the following numbers of individuals:\n", 
    (nrow(dat_filt_final) - table(dat_filt_final$fold)))

###Prepare output
indiv_list <- cbind(dat_filt_final$id, dat_filt_final$id)

###Write out output
saveRDS(dat_filt_final, file=output_pheno_data)

fwrite(x=indiv_list, file=output_inds_list, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE, showProgress=FALSE)

