# Check instrument for BUN & eGFR in UKB
# set parameters
params<- list(exposure = "BUN",
              exposure_family = "linear", # logistic or linear
              covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
              covar.factor = c("sex"),
              path.input = "",
              path.genome = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/backingfile/chr_1_22_all.rds",
              path.sample = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample",
              path.sample_stats = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/sample-stats.txt",
              path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
              path.names_out_prefix =  "Black_assoc_", # will add outcome & split#
              path.instrument = "BUN",
              instrument.name = "BUN"
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

## data input

path.genome<- paste0(params$path.input,params$path.genome)
path.sample<- paste0(params$path.input,params$path.sample)
path.sample_stats<-paste0(params$path.input,params$path.sample_stats)

# data output

path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix,"_inst", params$instrument.name, "_exp",  var.name,".csv")



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


####### to run a seocnd time without reloading


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







# make sure it is numeric!
covar.list=union(covar.num, covar.factors)

sample<- sample %>%
  mutate(across(all_of(covar.num), as.double)) %>%
  mutate(across(all_of(covar.factors), as.factor))  %>%
  mutate(across(all_of(var.name), as.double)) 






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

