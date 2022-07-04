# create backing files for black cohort in UKB
library(bigsnpr)


# load lists of snps to include
path="/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/"
list_snp_id <- lapply(1:22, function(chr) {
  infos_chr <- bigreadr::fread2(paste0(path, 
                                       "chunks/bgens/black_MAFf_chr",
                                       chr,
                                       "_snpstats.txt"), 
                                showProgress = FALSE, 
                                skip=6, 
                                select=c(3:6),
                                nrow=20,
                                col.names=c("chr","pos","a0", "a1"))
  with(infos_chr, paste(chr, pos, a0, a1, sep = "_"))
})

sum(lengths(list_snp_id)) #5053044



ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
ncores
if (is.na(ncores)) { 
  ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
} 
ncores
if (is.na(ncores)) { 
  return(2) # default
} 
ncores<-ncores-1


# do 1 by 1 to find faulty one

for (j in 1:22) {
paste0(j, " started") 

  rds <- bigsnpr::snp_readBGEN(
    bgenfiles   = paste0(path, "chunks/bgens/black_MAFf_chr",j,".bgen"),
    list_snp_id = list_snp_id[j],
    backingfile = paste0(path, "chunks/backingfile/chr_1_22_",j,"_test"),
    ind_row     = 1:20,
    ncores      = ncores
  )

}

bgens<- glue::glue("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/bgens/black_MAFf_chr{chr}.bgen",
                   chr = 1:22)

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles   = bgens,
    list_snp_id = list_snp_id,
    backingfile = paste0(path, "chunks/backingfile/chr_1_22"),
    ind_row     = ,
    ncores      = ncores
  )
) # 43 min with 23 cores

"end of script"

