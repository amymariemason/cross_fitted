# load bigsnp package
library(bigsnpr)

path="/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/"

# load bgen data 
######### ONLY RUN THIS ONCE##########
snp_stats<-bigreadr::fread2(paste0(path,"chunks/black_chr_21_chunk2_snpstats"),
                                select=c(3:6), 
                                col.names=c("chr","pos","a0", "a1"))


list_snp_id <- with(snp_stats, split(paste(chr, pos, a0, a1, sep = "_"),
                                   factor(chr, levels = 1:22)))





snp_readBGEN(paste0(path,"chunks/black_chr_21_chunk2.bgen"), 
             paste0(path,"chunks/backingfile/black_chr_21_chunk2.bgen"), 
             list_snp_id = list_snp_id[21])

# 


# load bgen data 
######### ONLY RUN THIS ONCE##########
snp_stats<-bigreadr::fread2(paste0(path,"data/snp-stats-21.txt"),skip=6,
                            select=c(3:6), 
                            col.names=c("chr","pos","a0", "a1"))


list_snp_id <- with(snp_stats, split(paste(chr, pos, a0, a1, sep = "_"),
                                     factor(chr, levels = 1:22)))





snp_readBGEN(paste0(path,"data/black_chr21.bgen"), 
             paste0(path,"data/backingfile/black_chr21.bgen"), 
             list_snp_id = list_snp_id[21])

# failed, but only because wrong output format

