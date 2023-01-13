---
title: "Cross fitting in black cohort report"
author: "Amy Mason"
date: "2022-08-24"
output: html_document
params: 
  exposure: "ln_uacr"
  outcome: "cad_int"
  outcome_family: "binomial"
  sig_cutoffs: !r 5e-7
  covar.num: !r c("ages","PC1","PC2","PC3","PC4", "PC5", "PC6","PC7","PC8","PC9", "PC10")
  covar.factor: !r c("sex")
  splits: 20
  path.input: ""
  path.genome: "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/chunks/backingfile/chr_1_22_all.rds"
  path.sample: "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/black.sample"
  path.sample_stats: "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/data/sample-stats.txt"
  save.all.output: FALSE # saves all GWAS for all subsamples
  save.topGWAS.output: TRUE # saves top 10000 hits on GWAS for each sample
  save.instrument: TRUE # saves subgroup numbers & instrument
  save.sample: FALSE # saves original sample with subgroup numbers and instrument
  path.output: "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/"
  path.names_out_prefix: "Black_CXfit_" # will add outcome & split#
  # phewas/covariate association options
  PHEWAS.check: TRUE # checks snps selected with phewas, remove if reach covar.cutoff 
  PHEWAS.remove: !r NULL # remove snps that match on this trait 
  covariate.check: TRUE # run GWAS on black cohort for following covariates, remove if reach covar.cutoff
  covariate.remove: !r NULL # must be in sample
  covar.cutoff: !r NULL # defaults to p_threshold if null
  snp.remove: !r NULL # to force removal of snps that are correlated with other causal factors
  
---

















