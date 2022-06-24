# outputs: chunk list of snps for use in snptest GWAS in format output_file/MAF_chr_i_chunk_j
# 

# CHANGE THIS TO CHANGE CHUNK SIZE
chunk_size<-100000
## CHANGE THESE LINES TO CHANGE MAF CUTOFF
MAF_cut<-0.01
## CHANGE THESE LINES TO CHANGE INFO CUTOFF
info_cut<-0.5
## HW CUTOFF
HW_cut<-1*10^{-6}

## CHANGE TO YOUR PREFERED OUTPUT DIRECTORY. WARNING: this makes a lot of small files.
Output_file<-"/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/"

########################################################################################
library(bigreadr)

for (i in 3:22){
  print(paste0("Reading chromosome ",i))
  stats = bigreadr::fread2(
    paste0("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black_chr",i,"_snpstats.txt"),
                     stringsAsFactors=FALSE, header=TRUE, skip=7)
  
  snps = stats$rsid[stats$minor_allele_frequency>MAF_cut & 
                      stats$info>info_cut & stats$HW_exact_p_value>HW_cut]
  snppos= stats$position[stats$rsid %in% snps]
  snpneg= stats$rsid[!stats$rsid %in% snps]
  print(paste0("Filtering chr ",i," for MAF >", MAF_cut))
  length<-floor(length(snps)/chunk_size)
  # wipe control file
  sink(paste0(Output_file,"MAF_chr_",i,"control.txt"), append=FALSE)
  sink()
  # create new exclusion file
  write.table(snpneg,paste0(Output_file,"MAF_chr_",i,"exclude.txt"),quote = F, row.names = F, col.names = F)
  write.table(snps,paste0(Output_file,"MAF_chr_",i,"include.txt"),quote = F, row.names = F, col.names = F)
  

      
      for (j in 0:length){
        test<-snps[(j*chunk_size+1):min((j*chunk_size+chunk_size),length(snps))]
        test<-na.omit(test)
        # print variables
        printchr<-ifelse(i<10,paste0("0",i),paste0(i))
        printstart<-snppos[(j*chunk_size+1)]-10
        printend<- snppos[min((j*chunk_size+chunk_size),length(snps))]+10
        # print this chunk  
        sink(paste0(Output_file,"MAF_chr_",i,"_chunk_",j+1, ".txt"), append=FALSE)
        cat(test, sep=" ")
        sink()
        # add range of this chunk to control document
        sink(paste0(Output_file,"MAF_chr_",i,"control.txt"), append=TRUE)
        cat(paste0(printchr,":",printstart,"-",printend, "\n"))
        sink()
        if(j%%50==0){
          print(paste0("Creating files ", round(j/length*100), "% complete"))
        }
      }
    }  
  
print("Spliting snplist complete")  

