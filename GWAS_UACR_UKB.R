# Run GWAS on black cohort of UKB

# set parameters
params<- list(exposure = "BUN",
  transform = "log",  # log, scale
  exposure_family = "linear", # setting not implemented, change if binary exposure
  sig_cutoffs=  5e-7,
  covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
  covar.factor = c("sex"),
  path.input = "",
  path.genome = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/backingfile/chr_1_22_all.rds",
  path.sample = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample",
  path.sample_stats = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/sample-stats.txt",
  save.all.output = TRUE, # saves all GWAS for all subsamples
  save.topGWAS.output = TRUE, # saves top 10000 hits on GWAS for each sample
  path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
  path.names_out_prefix =  "Black_GWAS_", # will add outcome & split#
  # phewas/covariate association options
  PHEWAS.check = TRUE # checks snps selected with phewas, remove if reach covar.cutoff 
  )
  

# load packages
library(bigsnpr)
library(dplyr)
library(ieugwasr)


#inputs

var.name<-params$exposure
covar.num=params$covar.num
# covariates will be treated as numeric unless in this list
covar.factors = params$covar.factor

exposure_family<-params$exposure_family

## data input

path.genome<- paste0(params$path.input,params$path.genome)
path.sample<- paste0(params$path.input,params$path.sample)
path.sample_stats<-paste0(params$path.input,params$path.sample_stats)

# data output

path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_alldata.csv")
path.GWAS_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_topGWAS.csv")
path.GWAS_manhattan_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_manhattan.png")
path.GWAS_qq_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_qq.png")


# plieotropy checks

PHEWAS.check= params$PHEWAS.check # checks snps selected with phewas, prints table for each 



###################################
# setup data
# load  stored genome data

bigSNP <-snp_attach(params$path.genome)


# Get aliases for useful slots
G   <- bigSNP$genotypes
CHR <- as.integer(bigSNP$map$chromosome) 
POS <- bigSNP$map$physical.pos
rsids<- bigSNP$map$rsid
# Check some counts for the 10 first SNPs
big_counts(G, ind.col = 1:10)

# create snp map

map <- transmute(bigSNP$map,
         chr = as.integer(chromosome), 
         pos = physical.pos,
         a0 = allele1, a1 = allele2) 

# Load sample data
sample <- read.csv(params$path.sample, sep="")
sample<-sample[-1,]

# QC on samples
sample.stats <- read.table(params$path.sample_stats, header=TRUE, quote="\"")
assertthat::assert_that(max(sample.stats$missing_call_proportion)<0.03)


# make sure it is numeric!
covar.list=union(covar.num, covar.factors)

sample<- sample %>%
  mutate(across(all_of(covar.num), as.double)) %>%
  mutate(across(all_of(covar.factors), as.factor))  %>%
  mutate(across(all_of(var.name), as.double)) 


# remove people with missing data

print(paste0("sample size is ", nrow(sample)))

missing<-which(is.na(sample[,var.name]))

print("removing people with missing values")
print(paste0("outcome:", var.name, ": ",length(missing)))

sample2<-sample[!is.na(sample[,var.name]),]


##Clumping on snps (already somewhat filtered  -see step 4)

# Half of the cores you have on your computer
NCORES <- nb_cores()
# Exclude Long-Range Linkage Disequilibrium Regions of the human genome
# based on an online table. https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
excl1 <- snp_indLRLDR(infos.chr = CHR, infos.pos = POS)
# rsid does duplicate some snps so filter on MAF again at same level at step 5 
maf <- snp_MAF(bigSNP$genotypes)
excl2<- which(maf<0.01)
ind.excl = union(excl1, excl2)
# Use clumping (on the MAF) to keep SNPs weakly correlated with each other.
ind.keep <- snp_clumping(G, infos.chr = CHR,
                         thr.r2 = 0.2,
                         exclude = ind.excl,
                         ncores = NCORES)

# create subset of clumped variants
df_beta0<-bigSNP$map %>%
    transmute(rsid=rsid,chr=as.integer(chromosome), pos=as.integer(physical.pos),a0=allele1, a1=allele2, freq=freq)
df_beta2<-df_beta0[ind.keep,]




#inputs

var.name<-params$exposure
covar.num=params$covar.num
# covariates will be treated as numeric unless in this list
covar.factors = params$covar.factor

exposure_family<-params$exposure_family

## data input

path.genome<- paste0(params$path.input,params$path.genome)
path.sample<- paste0(params$path.input,params$path.sample)
path.sample_stats<-paste0(params$path.input,params$path.sample_stats)

# data output

path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_alldata.csv")
path.GWAS_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_topGWAS.csv")
path.GWAS_manhattan_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_manhattan.png")
path.GWAS_qq_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_qq.png")



###################################
# setup data

# make sure it is numeric!
covar.list=union(covar.num, covar.factors)



sample<- sample %>%
  mutate(across(all_of(covar.num), as.double)) %>%
  mutate(across(all_of(covar.factors), as.factor))  %>%
  mutate(across(all_of(var.name), as.double)) 
  
  if (params$transform=="log"){
    sample<- sample %>%
       mutate(across(all_of(var.name), log)) %>%
       mutate(across(all_of(var.name), scale)) 
    
path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, "lg_", var.name, "_GWAS_alldata.csv")
path.GWAS_data_out<- paste0(params$path.output, params$path.names_out_prefix,"lg_", var.name, "_GWAS_topGWAS.csv")
path.GWAS_manhattan_out<- paste0(params$path.output, params$path.names_out_prefix, "lg_", var.name, "_GWAS_manhattan.png")
path.GWAS_qq_out<- paste0(params$path.output, params$path.names_out_prefix, "lg_", var.name, "_GWAS_qq.png")
    }
  
    if (params$transform=="scale"){
    sample<- sample %>%
       mutate(across(all_of(var.name), scale)) 
    }
  



# remove people with missing data

print(paste0("sample size is ", nrow(sample)))

missing<-which(is.na(sample[,var.name]))

print("removing people with missing values")
print(paste0("outcome:", var.name, ": ",length(missing)))

sample2<-sample[!is.na(sample[,var.name]),]



# DO THE GWAS
  new_var= var.name 
  y_index =  which(!is.na(sample[,var.name]))

  if(exposure_family=="linear"){

  GWAS_output <- big_univLinReg(G,
    y.train = sample[y_index, new_var],
    ind.train = y_index,
    ind.col = ind.keep,
    covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
    thr.eigval = 1e-04,
    ncores = NCORES
  )
}
  
  if(exposure_family=="logistic"){

  GWAS_output <- big_univLogReg(G,
    y01.train = sample[y_index, new_var],
    ind.train = y_index,
    ind.col = ind.keep,
    covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
    ncores = NCORES
  )
  

  
}

  


  
  # get beta values from GWAS
  GWAS_output$p<-predict(GWAS_output,log10=FALSE)
  GWAS_output$man_p<- -predict(GWAS_output,log10=TRUE)
  GWAS_output$chr<-CHR[ind.keep]
  GWAS_output$pos<-POS[ind.keep]
  GWAS_output$rsid<-rsids[ind.keep]
  GWAS_output$eaf<- bigSNP$map$freq[ind.keep]
  GWAS_output$effect<-bigSNP$map$allele2[ind.keep]
  GWAS_output$other<-bigSNP$map$allele1[ind.keep]
  
    
# save output

  if(params$save.all.output==TRUE) {

write.table(GWAS_output, file = path.all_data_out)
  }
  

  if(params$save.topGWAS.output==TRUE) {
      
       GWAS_mini<-GWAS_output %>%
          arrange(p) %>%
          slice(1:10000)

write.table(GWAS_mini, file = path.GWAS_data_out)
    }

# save graphs
print(path.GWAS_manhattan_out)
  png(file=path.GWAS_manhattan_out)
    snp_manhattan(GWAS_output, CHR[ind.keep], POS[ind.keep], npoints=10000)
    dev.off()
# qq graph    
print(path.GWAS_qq_out)
  png(file=path.GWAS_qq_out) 
    snp_qq(GWAS_output)
    dev.off()


#############################################################
# edit to do additional GWAS in same run

exit

#############################################################


# Run GWAS on black cohort of UKB

# set parameters
params<- list(exposure = "eGFR",
  transform = "centre",  # log, centre
  exposure_family = "linear", 
  sig_cutoffs=  5e-7,
  covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
  covar.factor = c("sex"),
  path.input = "",
  path.genome = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/backingfile/chr_1_22_all.rds",
  path.sample = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample",
  path.sample_stats = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/sample-stats.txt",
  save.all.output = TRUE, # saves all GWAS for all subsamples
  save.topGWAS.output = TRUE, # saves top 10000 hits on GWAS for each sample
  path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
  path.names_out_prefix =  "Black_GWAS_", # will add outcome & split#
  # phewas/covariate association options
  PHEWAS.check = TRUE # checks snps selected with phewas, remove if reach covar.cutoff 
  )



#inputs

var.name<-params$exposure
covar.num=params$covar.num
# covariates will be treated as numeric unless in this list
covar.factors = params$covar.factor

exposure_family<-params$exposure_family

## data input

path.genome<- paste0(params$path.input,params$path.genome)
path.sample<- paste0(params$path.input,params$path.sample)
path.sample_stats<-paste0(params$path.input,params$path.sample_stats)

# data output

path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_alldata.csv")
path.GWAS_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_topGWAS.csv")
path.GWAS_manhattan_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_manhattan.png")
path.GWAS_qq_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_qq.png")



###################################
# setup data

# make sure it is numeric!
covar.list=union(covar.num, covar.factors)



sample<- sample %>%
  mutate(across(all_of(covar.num), as.double)) %>%
  mutate(across(all_of(covar.factors), as.factor))  %>%
  mutate(across(all_of(var.name), as.double)) 
  
  if (params$transform=="log"){
    sample<- sample %>%
       mutate(across(all_of(var.name), log)) %>%
       mutate(across(all_of(var.name), scale)) 
    }
  
    if (params$transform=="scale"){
    sample<- sample %>%
       mutate(across(all_of(var.name), scale)) 
    }
  



# remove people with missing data

print(paste0("sample size is ", nrow(sample)))

missing<-which(is.na(sample[,var.name]))

print("removing people with missing values")
print(paste0("outcome:", var.name, ": ",length(missing)))

sample2<-sample[!is.na(sample[,var.name]),]



# DO THE GWAS
  new_var= var.name 
  y_index =  which(!is.na(sample[,var.name]))

  if(exposure_family=="linear"){

  GWAS_output <- big_univLinReg(G,
    y.train = sample[y_index, new_var],
    ind.train = y_index,
    ind.col = ind.keep,
    covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
    thr.eigval = 1e-04,
    ncores = NCORES
  )
}
  
  if(exposure_family=="logistic"){

  GWAS_output <- big_univLogReg(G,
    y01.train = sample[y_index, new_var],
    ind.train = y_index,
    ind.col = ind.keep,
    covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
    ncores = NCORES
  )
}

  


  
  # get beta values from GWAS
  GWAS_output$p<-predict(GWAS_output,log10=FALSE)
  GWAS_output$man_p<- -predict(GWAS_output,log10=TRUE)
  GWAS_output$chr<-CHR[ind.keep]
  GWAS_output$pos<-POS[ind.keep]
  GWAS_output$rsid<-rsids[ind.keep]
  GWAS_output$eaf<- bigSNP$map$freq[ind.keep]
  GWAS_output$effect<-bigSNP$map$allele2[ind.keep]
  GWAS_output$other<-bigSNP$map$allele1[ind.keep]
  
    
# save output

  if(params$save.all.output==TRUE) {

write.table(GWAS_output, file = path.all_data_out)
  }
  

  if(params$save.topGWAS.output==TRUE) {
      
       GWAS_mini<-GWAS_output %>%
          arrange(p) %>%
          slice(1:10000)

write.table(GWAS_mini, file = path.GWAS_data_out)
    }

# save graphs
print(path.GWAS_manhattan_out)
  png(file=path.GWAS_manhattan_out)
    snp_manhattan(GWAS_output, CHR[ind.keep], POS[ind.keep], npoints=10000)
    dev.off()
# qq graph    
print(path.GWAS_qq_out)
  png(file=path.GWAS_qq_out) 
    snp_qq(GWAS_output)
    dev.off()
