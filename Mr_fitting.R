# set path 
path="/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/"

# load packages

library(ivreg)
library(ivtools)

#inputs

var.name<-"ln_uacr"
covar.num=c("ages",
            "PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10")
# covariates will be treated as numeric unless in this list
covar.factors = c("sex")
outcome<-"ckd"
outcome.family<-"binomial"
G.est=FALSE
boot=FALSE

## data input

path.sample_out <- paste0(path,"black_",var.name,"_instrument.sample")
path.data_out<- paste0(path,"black_",var.name,"_instrument_report.Rdata")

# Load sample data
sample <- read.csv(path.sample_out, sep="")


# working variables

covar.list=union(covar.num, covar.factors)

# ivreg method
if (outcome.family=="guassian") {
  
  sample<- sample %>%
    mutate(across(all_of(outcome), as.numeric)) 

formula= as.formula(paste(outcome, 
                          "~ ", var.name,
                          "|",  paste(covar.list, collapse="+"),
                          "| instrument"))

fit.ivreg<- ivreg::ivreg(formula, data=sample )
summary(fit.ivreg)

}

# iv tools method 

if (outcome.family=="binomial") {

sample<- sample %>%
  mutate(across(all_of(outcome), as.factor)) 


#two-stage estimation
formulaX.LZ= as.formula(paste(var.name, "~ instrument+", paste(covar.list, collapse="+")))
formulaY.LZ=as.formula(paste(outcome, "~", var.name, " +", paste(covar.list, collapse="+")))
fitX.LZ <- glm(formula=formulaX.LZ, data=sample)
fitY.LX <- glm(formula=formulaY.LZ, data=sample, family="binomial")
fitIV <- ivglm(estmethod="ts", fitX.LZ=fitX.LZ, fitY.LX=fitY.LX, data=sample,
               ctrl=FALSE)

summary(fitIV)

}


if (outcome.family=="binomial" & G.est==TRUE) {
## g-estimation method
formulaZ.L= as.formula(paste("instrument~", paste(covar.list, collapse="+")))
formulaY.LZX =as.formula(paste(outcome, "~ instrument +", var.name, " +", paste(covar.list, collapse="+")))

fitZ.L <- glm(formula=formulaZ.L, data=sample)
fitY.LZX <- glm(formula=formulaY.LZX, data=sample, family="binomial")

fitIV_g <- ivglm(estmethod="g", X="ln_uacr", fitZ.L=fitZ.L, fitY.LZX=fitY.LZX, data=sample, link="logit")


summary(fitIV_g)
}


# bootstrapped errors instead

if (outcome.family=="binomial" & boot==TRUE) {

library(boot)
set.seed(1)

bootfun <- function(data, indicies){ 
  dd <- data[indicies, ] 
  #two-stage estimation 
  fitX.LZ <- glm(formula=formulaX.LZ, data=dd)
  fitY.LX <- glm(formula=formulaY.LZ, data=dd, family="binomial")
  fitIV_ts <- ivglm(estmethod="ts", fitX.LZ=fitX.LZ, fitY.LX=fitY.LX, data=dd,
                    vcov.fit=FALSE, ctrl=FALSE)  
  est_ts <- fitIV_ts$est["ln_uacr"] 
  #G-estimation 
  fitZ.L <- glm(formula=formulaZ.L, data=dd)
  fitY.LZX <- glm(formula=formulaY.LZX, data=dd, family="binomial")
  fitIV_g <- ivglm(estmethod="g", X="ln_uacr", fitZ.L=fitZ.L, 
                   fitY.LZX=fitY.LZX, data=dd, link="logit", vcov.fit=FALSE) 
  est_g <- fitIV_g$est["vitd_std"] 
  return(c(est_ts, est_g)) 
} 

bb <- boot(data=sample, statistic=bootfun, R=100)
apply(bb$t, MARGIN=2, FUN=sd, na.rm=TRUE)

}

