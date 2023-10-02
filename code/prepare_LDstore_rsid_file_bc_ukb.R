library(bigsnpr)

###Load the genotype data
geno_dat <- "/scratch1/fabiom/ukb_bc_geno_imp_HM3.rds"
geno <- snp_attach(geno_dat)

###Extract SNPs info
geno_map <- as.data.frame(geno$map)

chrs <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
          "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
          "21", "22")

###Loop over chromosomes and format the data for LDstore 
for(i in chrs){
  chr_ind <- which(geno_map$chromosome==i)
  
  rsid_chr <- geno_map[chr_ind, 3]
  
  write.table(rsid_chr, paste0("../data/LDstore_input_files/ukb_bc_chr", as.integer(i), "_LDstore_rsid.txt"),
              col.names = FALSE, row.names = FALSE, quote=FALSE, sep=" ")
}