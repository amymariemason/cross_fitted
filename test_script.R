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

 
 # set default parameters
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
  
  # replace with arguments
  params$exposure<-argsL$exposure
  params$transform<-argsL$transform
  params$exposure_family<-argsL$exposure_family
  
 