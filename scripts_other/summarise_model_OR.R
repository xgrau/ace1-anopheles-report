summarise_model_report_string = function(model) {
  # report model
  mod_sum = summary(model)
  mod_sum_out = as.data.frame(mod_sum$coefficients)
  # odds ratios:
  mod_sum_OR = 1/exp(mod_sum_out$Estimate)[-1] # odds ratio for each genotype (exclude intercept)
  mod_sum_ORcu = 1/exp(mod_sum$coefficients[-1,1] + qnorm(0.975) * mod_sum$coefficients[-1,2])
  mod_sum_ORcd = 1/exp(mod_sum$coefficients[-1,1] + qnorm(0.025) * mod_sum$coefficients[-1,2])
  mod_sum_var = rownames(mod_sum_out)[-1] # list of genotypes
  # create a string with odds ratios and CIs
  mod_sum_ORstring = c() ; van = 0
  for (vai in mod_sum_var) {
    van = van + 1
    mod_sum_ORstring[van] = sprintf("%s OR = %.3f ( %.3f - %.3f )", vai, signif(mod_sum_OR[van],5),  signif(mod_sum_ORcu[van],5), signif(mod_sum_ORcd[van],5) )
  }
  mod_sum_ORstring = paste(mod_sum_ORstring, collapse = "\n")
  return(mod_sum_ORstring)
}
