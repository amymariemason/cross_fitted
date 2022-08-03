library(rmarkdown)

render_report = function(exposure, threshold, splits, covariate.remove) {
  rmarkdown::render(
    "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Rcode/bigsnp_crossfit_and_report.Rmd", 
      params = list(
        exposure=exposure,
        sig_cutoffs = threshold,
        splits=splits,
        covariate.remove = covariate.remove, #c("ldl", "trig", "bmi", "lpa", "sbp"),
        covar.cutoff = threshold,
        save.all.output = FALSE,
        save.sample = FALSE  
    ),
    output_file = paste0("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/CFMR-report_", exposure, "_threshold",threshold ,
                         "_splits", splits, "_date", Sys.Date(),".html")
  )
}

render_report("hdl",5e-6,20, c("ldl", "trig", "bmi", "lpa", "sbp"))
render_report("ldl",5e-6,20, c("hdl", "trig", "bmi", "lpa", "sbp"))
render_report("trig",5e-6,20, c("hdl", "ldl", "bmi", "lpa", "sbp"))
render_report("lpa",5e-6,20, c("hdl", "ldl", "bmi", "trig", "sbp"))
render_report("sbp",5e-6,20, c("hdl", "ldl", "bmi", "trig", "lpa"))
render_report("bmi",5e-6,20, c("hdl", "ldl", "sbp", "trig", "lpa"))



render_report_no_covar_rm = function(threshold, splits) {
  rmarkdown::render(
    "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Rcode/bigsnp_crossfit_and_report.Rmd", 
      params = list(
        sig_cutoffs = threshold,
        splits=splits,
        covar.cutoff = threshold,
        save.all.output = FALSE,
        save.sample = FALSE  
    ),
output_file = paste0("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/zz_mr/Amy/black/Output/CFMR-report_Xcv_", exposure, "_threshold",threshold
  )
}


#render_report_no_covar_rm(5e-6,20)

