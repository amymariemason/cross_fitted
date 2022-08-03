# set path 
path="/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/"

# load packages
library(bigsnpr)
library(dplyr)

#inputs

var.name<-"ln_uacr"
no_of_loops =20
covar.num=c("ages",
             "PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10")
# covariates will be treated as numeric unless in this list
covar.factors = c("sex")
sig_cutoff<- 5*10^(-2)

thr_clump = 0.001
size_clump = 100 /0.001

#### data locations

path.genome <-paste0(path,"chunks/backingfile/chr_1_22_test.rds")
path.sample<- paste0(path,"data/black.sample")
path.sample_stats<-paste0(path,"data/sample-stats.txt")

# data output
path.sample_out <- paste0(path,"black_",var.name,"_instrument.sample")
path.data_out<- paste0(path,"black_",var.name,"_instrument_report.Rdata")

###################################
# load  stored genome data

bigSNP <-snp_attach(path.genome)


# Get aliases for useful slots
G   <- bigSNP$genotypes
CHR <- as.integer(bigSNP$map$chromosome)
POS <- bigSNP$map$physical.pos
# Check some counts for the 10 first SNPs
big_counts(G, ind.col = 1:10)

# create snp map

map <- transmute(bigSNP$map,
         chr = as.integer(chromosome), 
         pos = physical.pos,
         a0 = allele1, a1 = allele2) 

# Load sample data
sample <- read.csv(path.sample, sep="")
sample<-sample[-1,]

# QC on samples
sample.stats <- read.table(path.sample_stats, header=TRUE, quote="\"")
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


# create empty outputs
GWAS_output<-  list(GWAS=vector("list", no_of_loops), 
                    MANHATTEN =  vector("list", no_of_loops),
                    QQ = vector("list", no_of_loops))
PRS_model <- vector("list", no_of_loops) 
instrument <- vector("double", nrow(sample))
report<-vector("list", no_of_loops) 


## create random sample sets

## add sample numbers at random
resample <- function(x, ...) x[sample.int(length(x), ...)]
randomiser<-resample(rep(1:no_of_loops,nrow(sample[-missing,])/no_of_loops+1)[1:nrow(sample[-missing, ])])

sample$sample_no<-NA
sample[-missing,]$sample_no<-randomiser



errormark=0
for (j in seq_len(no_of_loops)) {
  print(j)
  new_var= var.name 
  y_index <-  which(sample$sample_no!=j)
  y_validation <- which(sample$sample_no==j)


  obj.gwas <- big_univLinReg(G,
    y.train = sample[y_index, new_var],
    ind.train = y_index,
    ind.col = ind.keep,
    covar.train = covar_from_df(sample[y_index, covar.list]),  # already regressed on sex, ages, PC1-PC10
    thr.eigval = 1e-04,
    ncores = NCORES
  )

  GWAS_output$GWAS[[j]]<-obj.gwas
  GWAS_output$MANHATTEN[[j]] <- snp_manhattan(GWAS_output$GWAS[[j]], CHR[ind.keep], POS[ind.keep])
  GWAS_output$QQ[[j]]<-snp_qq(GWAS_output$GWAS[[j]])
  
  # get beta values from GWAS
  library(dplyr)
  df_beta1 <- GWAS_output$GWAS[[j]] %>%
    transmute(beta = estim, beta_se = std.err, n_eff = length(y_index))
  df_beta0<-bigSNP$map %>%
    transmute(rsid=rsid,chr=as.integer(chromosome), pos=as.integer(physical.pos),a0=allele1, a1=allele2, freq=freq)
  df_beta2<-df_beta0[ind.keep,]
  sum_stats<-cbind(df_beta2, df_beta1)
  sum_stats$p<-predict(GWAS_output$GWAS[[j]],log10=FALSE)
  sum_stats$man_p<- -predict(GWAS_output$GWAS[[j]],log10=TRUE)
  
  # weight snps by significance
  
  S <- rep(0, ncol(G)); 
  S[ind.keep] <- sum_stats$man_p

  # sig cut off

  signif=(S> -log10(sig_cutoff))
  
  # create not exclusion group - all of snps excluded before, any that are not significant
  
  ind.excl2= union(which(is.na(S)), which(!signif))
  
  ## check appropieteness of cutoff in first loop
  
  # clump again, prioritizing significance?
  
  ind.keep2 <- snp_clumping(G, 
                           infos.chr = map$chr, 
                           infos.pos = map$pos, 
                           S = S,
                           ind.row = y_index, 
                           thr.r2 = 0.001, 
                           size = 100 /0.001, 
                           ncores = NCORES,
                           exclude = ind.excl2 )
  
  if (is.null(ind.keep2)){
    errormark=1
    stop (paste0("no snps found at this p-value threshold in set" , j) )
  }
  df_beta_find <- snp_match(sum_stats, map[ind.keep2,],strand_flip = FALSE)[,c("chr", "pos", "a0", "a1","beta", "beta_se")]
  df_beta <- snp_match(df_beta_find, map,strand_flip = FALSE)
  

  # create genetic prediction for people in the application set i.e. who were NA in discovery set
  
  sample$PRS_value <- big_prodMat(X=G,as.matrix(df_beta$beta), ind.col = df_beta[["_NUM_ID_"]],
                            ncores = NCORES)
  
  
  # fit model on training set and use to predict values in validation set
  formula<- as.formula(paste(new_var, "~ PRS_value+", paste(covar.list, collapse="+")))
  formula2<- as.formula(paste(new_var, "~ PRS_value"))
  PRS_model[[j]]<- lm(formula, data = sample[y_index,])
  
  predict_test = predict(PRS_model[[j]], newdata=sample[y_validation,])
  instrument[y_validation]<- predict_test

  PRS_model[[j]]$error<- sqrt(mean((sample[y_validation,new_var]- predict_test)^2))
  names(PRS_model[[j]]$error)<-"RMSE"
  
  
  # report on each stage in quick summary
  
  report[[j]]<- c("sample" = j ,
                  "train_n"= length(y_index),
                  "predict_n" = length(y_validation),
                  "no_of_snps" = length(ind.keep2),
                  "PRS_score" = summary(PRS_model[[j]])$coef["PRS_value",3],
                  "lm_model_r2" = summary(PRS_model[[j]])$r.squared,
                  "RMSE" = PRS_model[[j]]$error)
 
  
}
print("GWAS and PRS report in subsets. See Rdata output for more details")
as.data.frame(do.call(rbind, report))

if (errormark==0) {

#save(GWAS_output, PRS_model,report, file = path.data_out)
save.image(file=path.data_out)

sample$instrument <- instrument

write.table(sample, path.sample_out, row.names=FALSE, quote=FALSE)

} else {
  print ("no result saved as errors found in GWAS loop")
}




