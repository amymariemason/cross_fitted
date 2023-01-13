# Run GWAS on black cohort of UKB

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1) {
 args <- c("--help")
}

## Help section

if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --exposure=   - character, names of outcome to GWAS in the sample file
      --transform=  - character: log, scale or empty : this transforms numeric outcomes as requested
      --exposure_family= - character: linear (continuous outcome) or logistic (binary outcome)
      --help              - print this text
 
      \n\n")
 
  q(save="no")
}

 # Default parameters
params<- list(exposure = "BUN",
  transform = "log",  # log, scale or ""
  exposure_family = "linear", # options = linear (continuous outcome) or logistic (binary outcome)
  covar.num = c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10"),
  covar.factor = c("sex"),
  path.input = "",
  path.genome = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/backingfile/chr_1_22_all.rds",
  path.sample = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample",
  path.sample_stats = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/sample-stats.txt",
  save.all.output = TRUE, # saves all GWAS for all subsamples
  save.topGWAS.output = FALSE, # saves top 10000 hits on GWAS for each sample
  path.output = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/",
  path.names_out_prefix =  "v2Black_GWAS_" # will add outcome 
  )
  
 
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF[,2]))
names(argsL) <- argsDF[,1]
 
## Arg1 default
if(is.null(argsL$exposure)) {
  cat("ERROR:exposure missing")
  q(save="no")
}
 
## Arg2 default
if(is.null(argsL$transform)) {
  argsL$transform<-""
}
 
## Arg3 default
if(is.null(argsL$exposure_family)) {
   argsL$exposure_family<-"linear"
}

   
  # replace defaults with arguments
  params$exposure<-argsL$exposure
  params$transform<-argsL$transform
  params$exposure_family<-argsL$exposure_family
  


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


# Load sample data
sample <- read.csv(params$path.sample, sep="")
sample<-sample[-1,]

# QC on samples
sample.stats <- read.table(params$path.sample_stats, header=TRUE, quote="\"")
assertthat::assert_that(max(sample.stats$missing_call_proportion)<0.03)


# make sure data is numeric or factor as appropiete
covar.list=union(covar.num, covar.factors)

sample<- sample %>%
  mutate(across(all_of(covar.num), as.double)) %>%
  mutate(across(all_of(covar.factors), as.factor))  %>%
  mutate(across(all_of(var.name), as.double)) 
  
  
 # rescale outcome if required
   
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



print("removing people with missing values")
missing<-which(is.na(sample[,var.name]))
print(paste0("outcome:", var.name, ": ",length(missing)))

print("removing people with missing covariates")
missing2<-which(!complete.cases(sample[,covar.list]) & !is.na(sample[,var.name])) 
print(paste0("outcome:", var.name, ": ",length(missing2)))

ind.gwas <- which(!is.na(sample[,var.name]) & complete.cases(sample[,covar.list]))


# Half of the cores you have on your computer
NCORES <- nb_cores()


# Exclude Long-Range Linkage Disequilibrium Regions of the human genome
# based on an online table. https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
excl1 <- snp_indLRLDR(infos.chr = CHR, infos.pos = POS)
# rsid does duplicate some snps so filter on MAF again at same level at step 5 
maf <- snp_MAF(bigSNP$genotypes)
excl2<- which(maf<0.01)
ind.excl = union(excl1, excl2)
# NB snps have already been filtered by MAF > 0.01, hwe > 1*10^-6 and info >0.5
ind.keep <- setdiff(1:ncol(G), ind.excl)               

# DO THE GWAS
  if(exposure_family=="linear"){

  GWAS_output <- big_univLinReg(G,
    y.train = sample[ind.gwas, var.name],
    ind.train = ind.gwas,
    ind.col = ind.keep,
    covar.train = covar_from_df(sample[ind.gwas, covar.list]),  # already regressed on sex, ages, PC1-PC10
    thr.eigval = 1e-04,
    ncores = NCORES
  )
}
  
  if(exposure_family=="logistic"){

  GWAS_output <- big_univLogReg(G,
    y01.train = sample[ind.gwas, var.name],
    ind.train = ind.gwas,
    ind.col = ind.keep,
    covar.train = covar_from_df(sample[ind.gwas, covar.list]),  # already regressed on sex, ages, PC1-PC10
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
  GWAS_output$effect<-bigSNP$map$allele1[ind.keep] #see https://github.com/privefl/bigsnpr/issues/249 & https://zzz.bwh.harvard.edu/plink/anal.shtml
  GWAS_output$other<-bigSNP$map$allele2[ind.keep]
  
    
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

