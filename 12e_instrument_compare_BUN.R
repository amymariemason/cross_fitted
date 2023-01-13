# Check instrument for BUN & eGFR in UKB




# set parameters
params<- list(exposure = "ischtia",
              exposure_family = "logistic", # setting not implemented, change if binary exposure
              covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
              covar.factor = c("sex"),
              path.input = "",
              path.genome = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/backingfile/chr_1_22_all.rds",
              path.sample = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample",
              path.sample_stats = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/sample-stats.txt",
              path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
              path.names_out_prefix =  "Black_assoc_", # will add outcome & split#
              path.instrument = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/UGR_at_0.00005",
              name.instrument ="UGR_BUN_v2"
              # checks snps selected with phewas, remove if reach covar.cutoff 
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

path.all_data_out<- paste0(params$path.output, params$path.names_out_prefix, "inst_",params$name.instrument, "exp_", var.name, ".csv")



## data input

path.genome<- paste0(params$path.input,params$path.genome)
path.sample<- paste0(params$path.input,params$path.sample)
path.sample_stats<-paste0(params$path.input,params$path.sample_stats)


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
BUN_UGR<- read.csv(params$path.instrument, sep="")%>%
  transmute(a0=allele0, a1=allele1, 
            beta=beta,  
            chr=chr, pos=ps, eaf=af, se=se)
# match up to genome
found<- which(POS %in% BUN_UGR$pos)
details<-tibble(snp=rsids[found], chr=CHR[found], pos=POS[found], 
                allele1=bigSNP$map$allele1[found], allele2=bigSNP$map$allele2[found])

BUN_UGR2<- merge(x=BUN_UGR, y=details, by=c("chr", "pos"), keep.y=TRUE, keep.x=FALSE)


sum_stats<-bigSNP$map %>%
  transmute(rsid=rsid,chr=as.integer(chromosome), pos=as.integer(physical.pos),a0=allele1, a1=allele2, freq=freq)


# check strand flipping before finishing this

match2<- snp_match(sumstats=BUN_UGR2, info_snp=sum_stats, strand_flip=FALSE, return_flip_and_rev = TRUE)


# check frequencies for sense
match2$eaf<-ifelse(match2$`_REV_`==TRUE, 1- match2$eaf, match2$eaf)
match2$palidrome<- (match2$a0 %in% c("C","G") & match2$a1 %in% c("C","G"))|(match2$a0 %in% c("A", "T") & match2$a1 %in% c("A", "T"))
match2$ambigeous_freq<-((match2$eaf>0.42 & match2$eaf<0.58)|(match2$freq>0.42 & match2$freq<0.58))
match2$mismatch_freq<-((match2$freq<0.42 & match2$eaf>0.58)|(match2$eaf<0.42 & match2$freq>0.58))


# flip if palidrom and mismatched direction
match2$beta.ss<- ifelse(match2$palidrome & match2$mismatch_freq,
                        -1*match2$beta, 
                        match2$beta)
match2$eaf.ss<- ifelse(match2$palidrome & match2$mismatch_freq,
                       1-match2$eaf, 
                       match2$eaf)


# remove if ambigous and palindrom
print("removing as ambigeous palindromes")
print(match2[(match2$palidrome & match2$ambigeous_freq),]$rsid )
match2<- match2[!(match2$palidrome & match2$ambigeous_freq),] 
print(nrow(match2))


# make instrument
sample$BUN_PRS_value <- big_prodMat(X=G,as.matrix(match2$beta), 
                                    ind.col = match2[["_NUM_ID_"]],
                                    ncores = nb_cores())

sample2<-sample[!is.na(sample[,var.name]),]


formula<- as.formula(paste("BUN~ BUN_PRS_value+", paste(covar.list, collapse="+")))

formula2<- as.formula(paste("BUN~", paste(covar.list, collapse="+")))


formula3<- as.formula("BUN~BUN_PRS_value")

BUN_PRS_model<- lm(formula, data = sample2)
BUN_PRS_model3<- lm(formula3, data = sample2)
summary(BUN_PRS_model)
summary(BUN_PRS_model3)
anova(BUN_PRS_model)

# try individual snps

snps<-G[, match2[,"_NUM_ID_"]]
sample2<-cbind(sample, snps)
colnames(snps)<-rsids[match2[,"_NUM_ID_"]]
sample2<-cbind(sample, snps)
sample3<-sample2[!is.na(sample2[,var.name]),]

formula4<- as.formula(paste("BUN~", paste(covar.list, collapse="+"),"+",paste(rsids[match2[,"_NUM_ID_"]], collapse="+")))

BUN_PRS_model4<- lm(formula4, data = sample3)
