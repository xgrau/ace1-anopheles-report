# libraries
source("../scripts_other/summarise_model_OR.R")
source("../scripts_other/glmodelling.R")

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



# standard error function
serr <- function(x) sd(x)/sqrt(length(x))



pdf(file="Fig_WApops_mean_CNV_variables_pop.pdf",height=4,width=3)
# plot number of samples
dac_aggregate_sum = aggregate(!is.na(dac$population), by = list(dac$population), FUN = sum)
ploti = barplot(height = dac_aggregate_sum$x,names.arg = dac_aggregate_sum$Group.1, horiz = T, las=1, xlab="num samples", col="slategray2",  main="num samples")
text(x=0,ploti,labels = dac_aggregate_sum$x, col="red",pos=4)
# plot means of ratios
plot(as.factor(dac$population), dac$ratio_FAM_HEX, las=2, ylab="280S allele ratio", col="slategray2", main="ratios")
dac_aggregate_mean = aggregate(dac$ratio_FAM_HEX, by = list(dac$population), FUN = mean)
dac_aggregate_serr = aggregate(dac$ratio_FAM_HEX, by = list(dac$population), FUN = serr)
ploti = barplot(height = dac_aggregate_mean$x,names.arg = dac_aggregate_mean$Group.1, horiz = T, las=1, xlab="280S allele ratio", col="slategray2",  main="ratios")
text(x=0,ploti,labels = paste(signif(dac_aggregate_mean$x, 3),"+/-", signif(dac_aggregate_serr$x,3), "SE"), col="red",pos=4)
# same with cnvs
plot(as.factor(dac$population), dac$CNV, las=2, ylab="copies", col="slategray2", main="Ace1 copies")
dac_aggregate_mean = aggregate(dac$CNV, by = list(dac$population), FUN = mean)
dac_aggregate_serr = aggregate(dac$CNV, by = list(dac$population), FUN = serr)
ploti = barplot(height = dac_aggregate_mean$x,names.arg = dac_aggregate_mean$Group.1, horiz = T, las=1, xlab="copies", col="slategray2",  main="Ace1 copies")
text(x=0,ploti,labels = signif(dac_aggregate_mean$x, 3), col="red",pos=4)
text(x=0,ploti,labels = paste(signif(dac_aggregate_mean$x, 3),"+/-", signif(dac_aggregate_serr$x,3), "SE"), col="red",pos=4)
dev.off()


#### Genotype-phenotype associations: CNVs ####

#### Population ####
# test association of FAMHEX ratio and CNVs with resistance independently
# and then use a joint model + BKE reduction

table_models_min = data.frame()
table_modvar_min = data.frame()
table_models_cnv = data.frame()
table_modvar_cnv = data.frame()
table_models_rat = data.frame()
table_modvar_rat = data.frame()

for (pop in pop_order) {
  
  # variables for model
  data = dac[dac$population == pop, c("phenotype","ratio_FAM_HEX","CNV")]
  
  # null model
  mod_nul = glm(phenotype ~ 1, data = data, family = "binomial")

  # full model
  mod_tot = glm(phenotype ~ CNV+ratio_FAM_HEX, data = data, family = "binomial")
  
  # reduce full model using BIC
  mod_min = step(mod_tot, direction = "both", steps = 1e6, trace = F, k = log(nrow(data))) # k=log(num_obs) for BIC
  
  # reduce full model using BKE
  mod_bke = glmodelling(input.table = data,list.of.markers = c("CNV","ratio_FAM_HEX"), 
                        rescolumn = "phenotype",glm.function = "glm")
                        
  # GLMs for each variable
  mod_cnv = glm(phenotype ~ CNV, data = data, family = "binomial")
  mod_rat = glm(phenotype ~ ratio_FAM_HEX, data = data, family = "binomial")
  
  # significance of all models
  mod_min_signif = glm_tables(model = mod_min, null = mod_nul, model_name = paste(pop,"BIC minimal model, binomial GLM"))
  mod_bke_signif = glm_tables(model = mod_bke$final.model, null = mod_nul, model_name = paste(pop,"BKE minimal model, binomial GLM"))
  mod_cnv_signif = glm_tables(model = mod_cnv, null = mod_nul, model_name = paste(pop,"number of Ace1 copies, binomial GLM"))
  mod_rat_signif = glm_tables(model = mod_rat, null = mod_nul, model_name = paste(pop,"280S allele ratio, binomial GLM"))
  mod_tot_signif = glm_tables(model = mod_tot, null = mod_nul, model_name = paste(pop,"full model, binomial GLM"))
  
  
  # write tables
  write.table(file="Fig_WApops_CNV_models_BIC.csv", t(mod_min_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_BIC.csv", mod_min_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_BIC.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
  
  # write tables
  write.table(file="Fig_WApops_CNV_models_BKE.csv", t(mod_bke_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_BKE.csv", mod_bke_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_BKE.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
  
  # write tables
  write.table(file="Fig_WApops_CNV_models_CNV.csv", t(mod_cnv_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_CNV.csv", mod_cnv_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_CNV.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
  
  # write tables
  write.table(file="Fig_WApops_CNV_models_RAT.csv", t(mod_rat_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_RAT.csv", mod_rat_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_RAT.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
  
  # write tables
  write.table(file="Fig_WApops_CNV_models_TOT.csv", t(mod_tot_signif$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_TOT.csv", mod_tot_signif$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
  write.table(file="Fig_WApops_CNV_models_TOT.csv", data.frame(), quote=FALSE, sep="\t", row.names=FALSE, append = T)
  
}



message("FI")
