# libraries

library(lme4)
source("../scripts_other/glmodelling.R")
source("../scripts_other/summarise_model_OR.R")


# input
dac = read.table("biobank_data_172subset_CNVs.csv", sep="\t", header = T)
dac = dac[dac$is_either_outlier == 0,]

# populations to test
dac$population = paste(dac$species, dac$Location)
pop_order = sort(unique(dac$population))

# empty tables for model data
table_models = data.frame()
table_modvar = data.frame()


# convert CNV estimate (number of copies relative to Kisumu, where wt copy number = 2) to CNV as in Lucas et al
dac$CNV = dac$CNV_estimate * 2 

#### Genotype-phenotype associations: CNVs ####

pdf(file="Fig3Esup_CNV_GLMs.pdf",height=9,width=4)


#### Pool ####
# GLM with random effects on the entire dataset
# using species and location as random effects

# variables for model
data = dac[, c("phenotype","ratio_FAM_HEX","CNV","Location","species")]

# null model 
# with species and location as random effects
mod_null = glmer(phenotype ~ (1 | Location) + (1 | species), data = data, family = "binomial")

# full model 
# with species and location as random effects
mod_full = glmer(phenotype ~ CNV + ratio_FAM_HEX + (1 | Location) + (1 | species),
                 data = data, family = "binomial")

# significance of the full model
mod_signif = anova(mod_full, mod_null, test = 'Chisq')
mod_pval_chisq_v_null_full = mod_signif$`Pr(>Chisq)`[2]



# BKE model
# with species and location as random effects
# NOT BETTER THAN THE FULL MODEL
mod_BKE = glmodelling(input.table = data, list.of.markers = c("ratio_FAM_HEX","CNV"), control.for = c("Location","species"), 
                      rescolumn = "phenotype", glm.function = "glmer", verbose = F)

# significance of the BKE model relative to null
mod_signif = anova(mod_BKE$final.model, mod_null, test = 'Chisq')
mod_variables = paste(rownames(mod_BKE$final.sig), collapse = ",")
mod_pval_chisq_v_null = mod_signif$`Pr(>Chisq)`[2]

# significance of each variable
mod_cnv = glmer(phenotype ~ CNV + (1 | Location) + (1 | species), data = data, family = "binomial")
mod_rat = glmer(phenotype ~ ratio_FAM_HEX + (1 | Location) + (1 | species), data = data, family = "binomial")
mod_signif_cnv = anova(mod_cnv, mod_null, test = 'Chisq')
mod_signif_rat = anova(mod_rat, mod_null, test = 'Chisq')
mod_signif_cnv_p = mod_signif_cnv$`Pr(>Chisq)`[2]
mod_signif_rat_p = mod_signif_rat$`Pr(>Chisq)`[2]

# some diagnostic plots
layout(c(1,2,3), widths = 2)
plot(data$phenotype, main="all", ylab="ratio", col="slategray2", names=NA,
     sub=sprintf(
       "Alive = %i | Dead = %i\n\nBKE~null: p = %.3E\n%s\n\nfull~null p = %.3E\n%s",
       sum(data$phenotype=="Alive"),
       sum(data$phenotype=="Dead"),
       mod_pval_chisq_v_null,
       summarise_model_report_string(mod_BKE$final.model),
       mod_pval_chisq_v_null_full,
       summarise_model_report_string(mod_full)
     )
)
plot(data$phenotype, data$CNV, main="all", ylab="CNV", ylim=c(0,15), xlab="phenotype", col="slategray2",
     sub=sprintf(
       "av CNV = %.3f +/- %.3f | p = %.3E (1-var)",
       mean(data$CNV), 
       sd(data$CNV)/sqrt(nrow(data)),
       mod_signif_cnv_p
     )
)
plot(data$phenotype, data$ratio_FAM_HEX, main="all", ylab="ratio", ylim=c(0,15), xlab="phenotype", col="slategray2",
     sub=sprintf(
       "av rat = %.3f +/- %.3f | p = %.3E (1-var)",
       mean(data$ratio_FAM_HEX), 
       sd(data$ratio_FAM_HEX)/sqrt(nrow(data)),
       mod_signif_rat_p
     )
)

# save model table
mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = paste("allpops","CNV+ratio"), model_type = "glmer")
write.table(file="Fig3Esup_CNV_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = F)
write.table(file="Fig3Esup_CNV_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)



#### Species ####
# now do the same thing within species, 
# using location as a random effect
for (pop in c("col","gam")){
  
  # variables for model
  data = dac[dac$species == pop, c("phenotype","ratio_FAM_HEX","CNV","Location")]
  
  # null model
  mod_null = glmer(phenotype ~ 1|Location, data = data, family = "binomial")
  
  # # Full model with random effects (UNUSED)
  mod_full = glmer(phenotype ~ CNV + ratio_FAM_HEX + (1|Location), data = data, family = "binomial")
  mod_signif = anova(mod_full, mod_null, test = 'Chisq')
  mod_pval_chisq_v_null_full = mod_signif$`Pr(>Chisq)`[2]
  
  # BKE model with random effects (Location)
  mod_BKE = glmodelling(input.table = data, list.of.markers = c("ratio_FAM_HEX","CNV"), control.for = c("Location"),
                        rescolumn = "phenotype", glm.function = "glmer", verbose = F)
  mod_signif = anova(mod_BKE$final.model, mod_null, test = 'Chisq')
  mod_variables = paste(rownames(mod_BKE$final.sig), collapse = ",")
  mod_pval_chisq_v_null = mod_signif$`Pr(>Chisq)`[2]
  
  # significance of each variable
  mod_cnv = glmer(phenotype ~ CNV + (1 | Location), data = data, family = "binomial")
  mod_rat = glmer(phenotype ~ ratio_FAM_HEX + (1 | Location), data = data, family = "binomial")
  mod_signif_cnv = anova(mod_cnv, mod_null, test = 'Chisq')
  mod_signif_rat = anova(mod_rat, mod_null, test = 'Chisq')
  mod_signif_cnv_p = mod_signif_cnv$`Pr(>Chisq)`[2]
  mod_signif_rat_p = mod_signif_rat$`Pr(>Chisq)`[2]
  
  # some diagnostic plots
  layout(c(1,2,3), widths = 2)
  plot(data$phenotype, main=pop, ylab="ratio", col="slategray2", names=NA,
       sub=sprintf(
         "Alive = %i | Dead = %i\n\nBKE~null: p = %.3E\n%s\n\nfull~null p = %.3E\n%s",
         sum(data$phenotype=="Alive"),
         sum(data$phenotype=="Dead"),
         mod_pval_chisq_v_null,
         summarise_model_report_string(mod_BKE$final.model),
         mod_pval_chisq_v_null_full,
         summarise_model_report_string(mod_full)
       )
  )
  plot(data$phenotype, data$CNV, main=pop, ylab="CNV", ylim=c(0,15), xlab="phenotype", col="slategray2",
       sub=sprintf(
         "av CNV = %.3f +/- %.3f | p = %.3E (1-var)",
         mean(data$CNV), 
         sd(data$CNV)/sqrt(nrow(data)),
         mod_signif_cnv_p
       )
  )
  plot(data$phenotype, data$ratio_FAM_HEX, main=pop, ylab="ratio", ylim=c(0,15), xlab="phenotype", col="slategray2",
       sub=sprintf(
         "av rat = %.3f +/- %.3f | p = %.3E (1-var)",
         mean(data$ratio_FAM_HEX), 
         sd(data$ratio_FAM_HEX)/sqrt(nrow(data)),
         mod_signif_rat_p
       )
  )
  
  # save model table
  mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = paste(pop,"CNV+ratio"), model_type = "glmer")
  write.table(file="Fig3Esup_CNV_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
  write.table(file="Fig3Esup_CNV_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
  
  
  
}



#### Population ####
# now do the same thing for each population independently
# no random effects
for (pop in pop_order) {
  
  # variables for model
  data = dac[dac$population == pop, c("phenotype","ratio_FAM_HEX","CNV")]
  
  # null model
  mod_null = glm(phenotype ~ 1, data = data, family = "binomial")

  # full model (UNUSED)
  mod_full = glm(phenotype ~ CNV+ratio_FAM_HEX, data = data, family = "binomial")
  mod_signif = anova(mod_full, mod_null, test = 'Chisq')
  mod_pval_chisq_v_null_full = mod_signif$`Pr(>Chi)`[2]
  
  # BKE model
  mod_BKE = glmodelling(input.table = data, list.of.markers = c("ratio_FAM_HEX","CNV"), 
                        rescolumn = "phenotype", glm.function = "glm", verbose = F)
  
  # full v null
  mod_signif = anova(mod_BKE$final.model, mod_null, test = 'Chisq')
  mod_pval_chisq_v_null = mod_signif$`Pr(>Chi)`[2]
  mod_variables = paste(rownames(mod_BKE$final.sig), collapse = ",")
  
  # significance of each variable
  mod_cnv = glm(phenotype ~ CNV, data = data, family = "binomial")
  mod_rat = glm(phenotype ~ ratio_FAM_HEX, data = data, family = "binomial")
  mod_signif_cnv = anova(mod_cnv, mod_null, test = 'Chisq')
  mod_signif_rat = anova(mod_rat, mod_null, test = 'Chisq')
  mod_signif_cnv_p = mod_signif_cnv$`Pr(>Chi)`[2]
  mod_signif_rat_p = mod_signif_rat$`Pr(>Chi)`[2]
  
  
  # some diagnostic plots
  layout(c(1,2,3), widths = 2)
  plot(data$phenotype, main=pop, ylab="ratio", col="slategray2", names=NA,
       sub=sprintf(
         "Alive = %i | Dead = %i\n\nBKE~null: p = %.3E\n%s\n\nfull~null p = %.3E\n%s",
         sum(data$phenotype=="Alive"),
         sum(data$phenotype=="Dead"),
         mod_pval_chisq_v_null,
         summarise_model_report_string(mod_BKE$final.model),
         mod_pval_chisq_v_null_full,
         summarise_model_report_string(mod_full)
       )
  )
  plot(data$phenotype, data$CNV, main=pop, ylab="CNV", ylim=c(0,15), xlab="phenotype", col="slategray2",
       sub=sprintf(
         "av CNV = %.3f +/- %.3f | p = %.3E (1-var)",
         mean(data$CNV), 
         sd(data$CNV)/sqrt(nrow(data)),
         mod_signif_cnv_p
       )
  )
  plot(data$phenotype, data$ratio_FAM_HEX, main=pop, ylab="ratio", ylim=c(0,15), xlab="phenotype", col="slategray2",
       sub=sprintf(
         "av rat = %.3f +/- %.3f | p = %.3E (1-var)",
         mean(data$ratio_FAM_HEX), 
         sd(data$ratio_FAM_HEX)/sqrt(nrow(data)),
         mod_signif_rat_p
       )
  )
  
  # save model table
  mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = paste(pop,"CNV+ratio"), model_type = "glm")
  write.table(file="Fig3Esup_CNV_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
  write.table(file="Fig3Esup_CNV_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
  
}

dev.off()


message("FI")
