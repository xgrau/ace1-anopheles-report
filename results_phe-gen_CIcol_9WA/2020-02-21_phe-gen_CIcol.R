library(gmodels)
library(pheatmap)
library(stringr)
library(pairwiseCI)
source("../scripts_other/summarise_model_OR.R")

col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))

graphics.off()

# read input (run 2020-02-17_heatmaps.R to produce this table)
gtd = read.table("../results_tables/freq_CIcol_CNV-ALTallele.csv", header = T)


gtd$genotype = stringr::str_replace(string = gtd$genotype, pattern = "wt/wt", replacement = "GG")
gtd$genotype = stringr::str_replace(string = gtd$genotype, pattern = "119S/wt", replacement = "GS")
gtd$genotype = stringr::str_replace(string = gtd$genotype, pattern = "119S/119S", replacement = "SS")

data_A65S = read.table("Genotype_3489405_A65S.csv", header = T)
data_A65S$genotype = as.factor(as.character(data_A65S$genotype))


#### G280S genotype ####

pdf(file="Fig3A_CIcol_phe-gty.pdf",height=4,width=4)

# Fisher test
ta = CrossTable(gtd[gtd$population == "CIcol",]$genotype, gtd[gtd$population == "CIcol",]$phenotype, fisher = T, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
ta_corr = pairwiseCI::Prop.or(ta$t[1,], ta$t[2,], conf.level = 0.95, CImethod = "Woolf")


# GLM model
# used to calculate ORs for the resistance genotypes
data = gtd[gtd$population == "CIcol",][, c("phenotype","genotype")]
colnames(data) = c("phenotype","Ace1")
mod_null = glm(phenotype ~ 1, data = data, family = "binomial")
mod_full = glm(phenotype ~ ., data = data, family = "binomial")
mod_signif = anova(mod_full, mod_null, test = 'Chisq')
# pval string
mod_pval_chisq_v_null = mod_signif$`Pr(>Chi)`[2]

# report
pheatmap(t(ta$t), color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~G280S\nFisher's exact test p = %.3E\nGLM p = %.3E\n%s\nWoolf approx OR: %.3f ( %.3f - %.3f)",
                      ta$fisher.ts$p.value,mod_pval_chisq_v_null, summarise_model_report_string(mod_full), 
                      1/ta_corr$estimate, 1/ta_corr$conf.int[2], 1/ta_corr$conf.int[1]))

# save model table
mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = "CIcol G280S~resistance")
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", data.frame(), quote = F, append = T)

# same but ONLY THOSE WITH multiple wt and 1 280S
# Fisher test
ta = CrossTable(gtd[gtd$population == "CIcol" & gtd$estimated_n_ALT <= 1 ,]$genotype, gtd[gtd$population == "CIcol" & gtd$estimated_n_ALT <= 1,]$phenotype, fisher = T, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
ta_corr = pairwiseCI::Prop.or(ta$t[1,], ta$t[2,], conf.level = 0.95, CImethod = "Woolf")


# GLM model
# used to calculate ORs for the resistance genotypes
data = gtd[gtd$population == "CIcol"  & gtd$estimated_n_ALT <= 1,][, c("phenotype","genotype")]
colnames(data) = c("phenotype","Ace1")
mod_null = glm(phenotype ~ 1, data = data, family = "binomial")
mod_full = glm(phenotype ~ ., data = data, family = "binomial")
mod_signif = anova(mod_full, mod_null, test = 'Chisq')
# pval string
mod_pval_chisq_v_null = mod_signif$`Pr(>Chi)`[2]

# save model table
mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = "CIcol G280S 1nALT~resistance")
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", data.frame(), quote = F, append = T)


# report
pheatmap(t(ta$t), color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~G280S (1 ALT allele only)\nFisher's exact test p = %.3E\nGLM p = %.3E\n%s\nWoolf approx OR: %.3f ( %.3f - %.3f)",
                      ta$fisher.ts$p.value,mod_pval_chisq_v_null, summarise_model_report_string(mod_full), 
                      1/ta_corr$estimate, 1/ta_corr$conf.int[2], 1/ta_corr$conf.int[1]))

# same with S65A
# Fisher
ta = CrossTable(data_A65S[data_A65S$population == "CIcol",]$genotype, data_A65S[data_A65S$population == "CIcol",]$phenotype, fisher = T, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
ta_corr = pairwiseCI::Prop.or(ta$t[1,], ta$t[2,], conf.level = 0.95, CImethod = "Woolf")


# GLM model
# used to calculate ORs for the resistance genotypes
data = data_A65S[data_A65S$population == "CIcol",][, c("phenotype","genotype")]
colnames(data) = c("phenotype","A65S")
mod_null = glm(phenotype ~ 1, data = data, family = "binomial")
mod_full = glm(phenotype ~ ., data = data, family = "binomial")
mod_signif = anova(mod_full, mod_null, test = 'Chisq')
# pval string
mod_pval_chisq_v_null = mod_signif$`Pr(>Chi)`[2]

# save model table
mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = "CIcol S65A~resistance")
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", data.frame(), quote = F, append = T)


# GLM
pheatmap(t(ta$t), color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~A65S\nFisher's exact test p = %.3E\nGLM p = %.3E\n%s\nWoolf approx OR: %.3f ( %.3f - %.3f)",
                      ta$fisher.ts$p.value,mod_pval_chisq_v_null, summarise_model_report_string(mod_full), 
                      1/ta_corr$estimate, 1/ta_corr$conf.int[2], 1/ta_corr$conf.int[1]))
dev.off()





#### CNV ####
# table genotype-phenotype associations: CNVs
pdf(file="Fig3BC_CIcol_phe-CNV_phe-nALT.pdf",height=4,width=4)

# Fisher test
ta = CrossTable(gtd[gtd$population == "CIcol",]$CNV, gtd[gtd$population == "CIcol",]$phenotype, fisher = T, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)

# GLM model
# used to calculate ORs for the resistance genotypes
data = gtd[gtd$population == "CIcol",][, c("phenotype","CNV")]
colnames(data) = c("phenotype","Ace1")
mod_null = glm(phenotype ~ 1, data = data, family = "binomial")
mod_full = glm(phenotype ~ ., data = data, family = "binomial")
mod_signif = anova(mod_full, mod_null, test = 'Chisq')
# pval string
mod_pval_chisq_v_null = mod_signif$`Pr(>Chi)`[2]

# save model table
mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = "CIcol CNV~resistance")
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", data.frame(), quote = F, append = T)

# report
pheatmap(t(ta$t), color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~CNV\nFisher's exact test p = %.3E\nGLM p = %.3E\n%s",
                      ta$fisher.ts$p.value,mod_pval_chisq_v_null, summarise_model_report_string(mod_full)))
barplot(t(ta$t), col= c("springgreen3","magenta3"), xlab = "CNV", ylim = c(0,30), las=1)





#### nALT ####
# Fisher test
ta = CrossTable(gtd[gtd$population == "CIcol",]$estimated_n_ALT, gtd[gtd$population == "CIcol",]$phenotype, fisher = T, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)

# GLM model
# used to calculate ORs for the resistance genotypes
data = gtd[gtd$population == "CIcol",][, c("phenotype","estimated_n_ALT")]
colnames(data) = c("phenotype","Ace1")
mod_null = glm(phenotype ~ 1, data = data, family = "binomial")
mod_full = glm(phenotype ~ ., data = data, family = "binomial")
mod_signif = anova(mod_full, mod_null, test = 'Chisq')
# pval string
mod_pval_chisq_v_null = mod_signif$`Pr(>Chi)`[2]

# save model table
mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = "CIcol nALT~resistance")
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", data.frame(), quote = F, append = T)


# report
pheatmap(t(ta$t), color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~nALT\nFisher's exact test p = %.3E\nGLM p = %.3E\n%s",
                      ta$fisher.ts$p.value,mod_pval_chisq_v_null, summarise_model_report_string(mod_full)))
barplot(t(ta$t), col= c("springgreen3","magenta3"), xlab = "nALT", ylim = c(0,30), las=1)


## both

plot(jitter(gtd[gtd$population == "CIcol",]$CNV, factor=0.6), 
     jitter(gtd[gtd$population == "CIcol",]$estimated_n_ALT, factor=0.6), 
     xlab="Ace1 copies",  ylab="280S copies", las=1, xlim=c(2,5),
     main="280S and CNV",
     col=c("springgreen3","magenta3")[gtd[gtd$population == "CIcol",]$phenotype])
abline(h=1.5, lty=2, col="red")
abline(v=3.5, lty=2, col="red")

plot(jitter(gtd[gtd$population == "CIcol",]$CNV, factor=0.1), 
     gtd[gtd$population == "CIcol",]$proportion, factor=0.6, 
     main="ratio and CNV",
     xlab="Ace1 copies",  ylab="Fraction 280S alleles", las=1, ylim=c(0,1),xlim=c(2,5),
     col=c("springgreen3","magenta3")[gtd[gtd$population == "CIcol",]$phenotype])
abline(v=3.5, lty=2, col="red")

dev.off()





#### GLM model with A65S and G280S data ####

source("../scripts_other/summarise_model_OR.R")

# manual genotypes:
# load data: genotypes in each sample
data = read.table("Genotype_3492074_G280S.csv", header = T)
colnames(data) = c("ox_code", "population", "G280S", "phenotype")

# add a65s
data_A65S = read.table("Genotype_3489405_A65S.csv", header = T)
data$A65S = data_A65S[,"genotype"]

# add dup data
data_gtd = gtd[,c("id", "CNV", "estimated_n_ALT")]
data = merge(data, data_gtd, by.x = "ox_code", by.y = "id") 

# clean
rownames(data) = data$ox_code
data = data[,c("A65S","G280S", "CNV", "estimated_n_ALT", "phenotype")]

# convert data to factors
data$A65S = as.factor(data$A65S)
data$G280S = as.factor(data$G280S)


# null model
mod_null = glm(phenotype ~ 1, data = data, family = "binomial")

### BIC PROCEDURE
# initial model (with all variables)
mod_tot = glm(phenotype ~ ., data = data, family = "binomial")
mod_cnv = glm(phenotype ~ CNV, data = data, family = "binomial")
mod_rat = glm(phenotype ~ estimated_n_ALT, data = data, family = "binomial")
# stepwise removal using BIC criterion
mod_BIC = step(mod_tot, direction = "both", steps = 1e6, trace = F, k = log(nrow(data))) # k=log(num_obs) for BIC

# save model table
mod_tau = glm_tables(model=mod_BIC, null=mod_null, model_name = "CIcol all~resistance (BIC)")
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", data.frame(), quote = F, append = T)

mod_tau = glm_tables(model=mod_tot, null=mod_null, model_name = "CIcol all~resistance (full)")
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
write.table(file="Fig3_CIcol_phe-gty_all_models.csv", data.frame(), quote = F, append = T)



# ### BACKWARD ELIMINATION
# mod_BKE = glmodelling(input.table = data, list.of.markers = c("A65S","G280S", "CNV", "estimated_n_ALT"),
#                       rescolumn = "phenotype", glm.function = "glm", verbose = F)
# 
# # save model table
# mod_tau = glm_tables(model=mod_BKE$final.model, null=mod_null, model_name = "CIcol all~resistance (BKE)")
# write.table(file="Fig3ABC_CIcol_phe_minGLM_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
# write.table(file="Fig3ABC_CIcol_phe_minGLM_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)


message("FI CIcol")

# # pearson test (linkage between mutations)
# cor.test(as.numeric(as.character(data$G280S)), as.numeric(as.character(data$A65S)), method = "pearson")
# cor.test(as.numeric(as.character(data$G280S)), as.numeric(data$CNV>2), method = "pearson")
