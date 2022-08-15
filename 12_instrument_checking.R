# Check instrument for BUN & eGFR in UKB
# set parameters
params<- list(exposure = "urea",
              exposure_family = "logistic", # setting not implemented, change if binary exposure
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
              path.names_out_prefix =  "Black_instrm_", # will add outcome & split#
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
sig_cutoff<- params$sig_cutoff



## data input

path.genome<- paste0(params$path.input,params$path.genome)
path.sample<- paste0(params$path.input,params$path.sample)
path.sample_stats<-paste0(params$path.input,params$path.sample_stats)

# data output

path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_GWAS_alldata.csv")



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



#input BUN instrument
library(readr)
BUN_CKDGen_formated_clumped <- read_table("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/BUN_CKDGen_formated_clumped")
BUN_CKDGEN_in_UGR <- read_table("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/BUN_CKDGEN_in_UGR")


# match up to genome
found<- which(rsids %in% BUN_CKDGen_formated_clumped$SNP)
details<-tibble(snp=rsids[found], chr=CHR[found], pos=POS[found])

BUN_CKDGEN_in_UGR<- merge(x=BUN_CKDGEN_in_UGR, y=details, by.x=c("chr","ps"), by.y=c("chr", "pos"))
BUN_CKDGEN_in_UGR <- BUN_CKDGEN_in_UGR %>%
  transmute(rsid=snp,a0=allele0, a1=allele1, 
            beta=beta, beta.UGR=beta, chr=chr, 
            pos=ps, eaf.UGR=af, se.UGR=se)

found<- which(rsids %in% BUN_CKDGEN_in_UGR$snp)
details<-tibble(snp=rsids[found], chr=CHR[found], pos=POS[found])

BUN_CKDGen_formated_clumped<- merge(x=BUN_CKDGen_formated_clumped, y=details, by.x="SNP", by.y="snp", keep.y=TRUE, keep.x=FALSE)
BUN_CKDGen_formated_clumped <- BUN_CKDGen_formated_clumped %>%
  transmute(rsid=SNP,a0=other_allele.exposure, a1=effect_allele.exposure, 
              beta=beta.exposure, beta.CKDgen=beta.exposure, 
            chr=chr, pos=pos, eaf.CKDgen=eaf.exposure, se.CKDgen=se.exposure)





sum_stats<-bigSNP$map %>%
  transmute(rsid=rsid,chr=as.integer(chromosome), pos=as.integer(physical.pos),a0=allele1, a1=allele2, freq=freq)


# check strand flipping before finishing this

match_snps<- snp_match(sumstats=BUN_CKDGen_formated_clumped, info_snp=sum_stats, strand_flip=FALSE, return_flip_and_rev = TRUE)
match_snps$palidrome<- (match_snps$a0 %in% c("C","G") & match_snps$a1 %in% c("C","G"))|(match_snps$a0 %in% c("A", "T") & match_snps$a1 %in% c("A", "T"))
match_snps$ambigeous_freq<-((match_snps$eaf.CKDgen>0.42 & match_snps$eaf.CKDgen<0.58)|(match_snps$freq>0.42 & match_snps$freq<0.58))
match_snps$remove<-(match_snps$palidrome & match_snps$ambigeous_freq)
includ.snps<-match_snps[match_snps$remove==FALSE,]$rsid

match_snps2<- snp_match(sumstats=BUN_CKDGEN_in_UGR, info_snp=sum_stats, strand_flip=FALSE, return_flip_and_rev = TRUE)
match_snps2$palidrome<- (match_snps2$a0 %in% c("C","G") & match_snps2$a1 %in% c("C","G"))|(match_snps2$a0 %in% c("A", "T") & match_snps2$a1 %in% c("A", "T"))
match_snps2$ambigeous_freq<-((match_snps2$eaf.UGR>0.42 & match_snps2$eaf.UGR<0.58)|(match_snps2$freq>0.42 & match_snps2$freq<0.58))
match_snps2$remove<-(match_snps2$palidrome & match_snps2$ambigeous_freq)
includ.snps2<-match_snps2[match_snps2$remove==FALSE,]$rsid

keep.snps<-intersect(includ.snps, include.snps2)

View(match_snps[(match_snps$a0 %in% c("C","G") & match_snps$a1 %in% c("C","G"))|(match_snps$a0 %in% c("A", "T") & match_snps$a1 %in% c("A", "T")),])

# make instrument
sample$BUN_PRS_value <- big_prodMat(X=G,as.matrix(match_snps$beta), 
                                ind.col = match_snps[["_NUM_ID_"]],
                                ncores = nb_cores())

sample2<-sample[!is.na(sample[,var.name]),]


formula<- as.formula(paste("urea~ BUN_PRS_value+", paste(covar.list, collapse="+")))

formula2<- as.formula(paste("urea~", paste(covar.list, collapse="+")))

BUN_PRS_model<- lm(formula, data = sample2)
summary(BUN_PRS_model)
anova(BUN_PRS_model)


## get estimates for cad_int, haem, ischtia

# set parameters :cad-int
params<- list(exposure = "cad_int",
              exposure_family = "logistic", # setting not implemented, change if binary exposure
              covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
              covar.factor = c("sex"),
              path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
              path.names_out_prefix =  "Black_assoc_" # will add outcome & split#
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

outcome.family<-params$outcome_family


# data output

path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_alldata.csv")

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



# create subset of clumped variants
ind.keep= match_snps$`_NUM_ID_`

df_beta0<-bigSNP$map %>%
  transmute(rsid=rsid,chr=as.integer(chromosome), pos=as.integer(physical.pos),a0=allele1, a1=allele2, freq=freq)
df_beta2<-df_beta0[ind.keep,]



# DO THE GWAS
new_var= params$exposure
y_index =  which(!is.na(sample[,var.name]))

if(params$exposure_family=="linear"){
  
  GWAS_output <- big_univLinReg(G,
                                y.train = sample[y_index, new_var],
                                ind.train = y_index,
                                ind.col = ind.keep,
                                covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
                                thr.eigval = 1e-04,
                                ncores = nb_cores()
  )
}

if(params$exposure_family=="logistic"){
  
  GWAS_output <- big_univLogReg(G,
                                y01.train = sample[y_index, new_var],
                                ind.train = y_index,
                                ind.col = ind.keep,
                                covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
                                ncores = nb_cores()
  )
}


# get beta values from GWAS
GWAS_output$p<-predict(GWAS_output,log10=FALSE)
GWAS_output$man_p<- -predict(GWAS_output,log10=TRUE)
GWAS_output$chr<-CHR[ind.keep]
GWAS_output$pos<-POS[ind.keep]
GWAS_output$rsid<-rsids[ind.keep]
GWAS_output$eaf<- bigSNP$map$freq[ind.keep]
GWAS_output$effect<-bigSNP$map$allele1[ind.keep]
GWAS_output$other<-bigSNP$map$allele2[ind.keep]



# save output  
  write.table(GWAS_output, file = path.all_data_out)


  ## get estimates for cad_int, haem, ischtia
  
  # set parameters :cad-int
  params<- list(exposure = "haem",
                exposure_family = "logistic", # setting not implemented, change if binary exposure
                covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
                covar.factor = c("sex"),
                path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
                path.names_out_prefix =  "Black_assoc_" # will add outcome & split#
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
  outcome.family<-params$exposure_family
  
  
  # data output
  
  path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_alldata.csv")
  
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
  
  
  
  # create subset of clumped variants
  ind.keep= match_snps$`_NUM_ID_`
  
  df_beta0<-bigSNP$map %>%
    transmute(rsid=rsid,chr=as.integer(chromosome), pos=as.integer(physical.pos),a0=allele1, a1=allele2, freq=freq)
  df_beta2<-df_beta0[ind.keep,]
  
  
  
  # DO THE GWAS
  new_var= params$exposure
  y_index =  which(!is.na(sample[,var.name]))
  
  if(params$exposure_family=="linear"){
    
    GWAS_output <- big_univLinReg(G,
                                  y.train = sample[y_index, new_var],
                                  ind.train = y_index,
                                  ind.col = ind.keep,
                                  covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
                                  thr.eigval = 1e-04,
                                  ncores = nb_cores()
    )
  }
  
  if(params$exposure_family=="logistic"){
    
    GWAS_output <- big_univLogReg(G,
                                  y01.train = sample[y_index, new_var],
                                  ind.train = y_index,
                                  ind.col = ind.keep,
                                  covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
                                  ncores = nb_cores()
    )
  }
  
  
  # get beta values from GWAS
  GWAS_output$p<-predict(GWAS_output,log10=FALSE)
  GWAS_output$man_p<- -predict(GWAS_output,log10=TRUE)
  GWAS_output$chr<-CHR[ind.keep]
  GWAS_output$pos<-POS[ind.keep]
  GWAS_output$rsid<-rsids[ind.keep]
  GWAS_output$eaf<- bigSNP$map$freq[ind.keep]
  GWAS_output$effect<-bigSNP$map$allele1[ind.keep]
  GWAS_output$other<-bigSNP$map$allele2[ind.keep]
  
  
  
  # save output  
  write.table(GWAS_output, file = path.all_data_out)
  
  ## get estimates for cad_int, haem, ischtia
  
  # set parameters :cad-int
  params<- list(exposure = "ischtia",
                exposure_family = "logistic", # setting not implemented, change if binary exposure
                covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
                covar.factor = c("sex"),
                path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
                path.names_out_prefix =  "Black_assoc_", # will add outcome & split#
                sample # sample file
                  )
  df_beta0=bigSNP$map %>%
    transmute(rsid=rsid,chr=as.integer(chromosome), pos=as.integer(physical.pos),a0=allele1, a1=allele2, freq=freq)
  
  
find_est<-function(exposure = "ischtia",
                   exposure_family = "logistic", # setting not implemented, change if binary exposure
                   covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
                   covar.factor = c("sex"),
                   path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
                   path.names_out_prefix =  "Black_assoc_", # will add outcome & split#
                   sample=sample, #parameters as above
                   G=G, # genome
                   ind.keep= match_snps$`_NUM_ID_`, # numbers of snps
                   ind.train = y_index, # numbers of samples
                   df_beta0= df_beta0
){
  #inputs
  
  var.name<-params$exposure
  covar.num=params$covar.num
  # covariates will be treated as numeric unless in this list
  covar.factors = params$covar.factor
  outcome.family<-params$exposure_family
  
  # data output
  
  path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_alldata.csv")
  
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
  
  
  
  # create subset of clumped variants

 df_beta2<-df_beta0[ind.keep,]
  
  
  
  # DO THE GWAS
  new_var= params$exposure
  y_index =  which(!is.na(sample[,var.name]))
  
  if(params$exposure_family=="linear"){
    
    GWAS_output <- big_univLinReg(G,
                                  y.train = sample[y_index, new_var],
                                  ind.train = y_index,
                                  ind.col = ind.keep,
                                  covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
                                  thr.eigval = 1e-04,
                                  ncores = nb_cores()
    )
  }
  
  if(params$exposure_family=="logistic"){
    
    GWAS_output <- big_univLogReg(G,
                                  y01.train = sample[y_index, new_var],
                                  ind.train = y_index,
                                  ind.col = ind.keep,
                                  covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
                                  ncores = nb_cores()
    )
  }
  
  
  # get beta values from GWAS
  GWAS_output$p<-predict(GWAS_output,log10=FALSE)
  GWAS_output$man_p<- -predict(GWAS_output,log10=TRUE)
  GWAS_output$chr<-CHR[ind.keep]
  GWAS_output$pos<-POS[ind.keep]
  GWAS_output$rsid<-rsids[ind.keep]
  GWAS_output$eaf<- bigSNP$map$freq[ind.keep]
  GWAS_output$effect<-bigSNP$map$allele1[ind.keep]
  GWAS_output$other<-bigSNP$map$allele2[ind.keep]
  
  
  # save output  
  write.table(GWAS_output, file = path.all_data_out)
  
}
  
  

  
  
  
  # set parameters : CKD
  params<- list(exposure = "ckd",
                exposure_family = "logistic", # setting not implemented, change if binary exposure
                covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
                covar.factor = c("sex"),
                path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
                path.names_out_prefix =  "Black_assoc_" # will add outcome & split#
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
  outcome.family<-params$exposure_family
  
  
  # data output
  
  path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, var.name, "_alldata.csv")
  
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
  
  
  
  # create subset of clumped variants
  ind.keep= match_snps$`_NUM_ID_`
  
  df_beta0<-bigSNP$map %>%
    transmute(rsid=rsid,chr=as.integer(chromosome), pos=as.integer(physical.pos),a0=allele1, a1=allele2, freq=freq)
  df_beta2<-df_beta0[ind.keep,]
  
  
  
  # DO THE GWAS
  new_var= params$exposure
  y_index =  which(!is.na(sample[,var.name]))
  
  if(params$exposure_family=="linear"){
    
    GWAS_output <- big_univLinReg(G,
                                  y.train = sample[y_index, new_var],
                                  ind.train = y_index,
                                  ind.col = ind.keep,
                                  covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
                                  thr.eigval = 1e-04,
                                  ncores = nb_cores()
    )
  }
  
  if(params$exposure_family=="logistic"){
    
    GWAS_output <- big_univLogReg(G,
                                  y01.train = sample[y_index, new_var],
                                  ind.train = y_index,
                                  ind.col = ind.keep,
                                  covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
                                  ncores = nb_cores()
    )
  }
  
  
  # get beta values from GWAS
  GWAS_output$p<-predict(GWAS_output,log10=FALSE)
  GWAS_output$man_p<- -predict(GWAS_output,log10=TRUE)
  GWAS_output$chr<-CHR[ind.keep]
  GWAS_output$pos<-POS[ind.keep]
  GWAS_output$rsid<-rsids[ind.keep]
  GWAS_output$eaf<- bigSNP$map$freq[ind.keep]
  GWAS_output$effect<-bigSNP$map$allele1[ind.keep]
  GWAS_output$other<-bigSNP$map$allele2[ind.keep]
  
  
  
  # save output  
  write.table(GWAS_output, file = path.all_data_out)
  