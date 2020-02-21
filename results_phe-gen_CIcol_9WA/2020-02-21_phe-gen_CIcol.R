library(gmodels)
library(pheatmap)

col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))

# read input (run 2020-02-17_heatmaps.R to produce this table)
gtd = read.table("../results_tables/Fig3_CIcol_CNV-ALTallele.csv", header = T)


#### genotype-phenotype associations in CIcol ####
# Simple Fisher's tests of individual variants

# table genotype-phenotype associations: 119S
pdf(file="Fig3_CIcol_phe-gty119S.pdf",height=12,width=12)
ta = CrossTable(gtd[gtd$population == "CIcol",]$genotype, gtd[gtd$population == "CIcol",]$phenotype, fisher = T, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
pheatmap(t(ta$t), color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~119S in CIcol\nFisher's exact test p=%.3E", ta$fisher.ts$p.value))
dev.off()

# table genotype-phenotype associations: CNVs
pdf(file="Fig3_CIcol_phe-CNV_phe-nALT.pdf",height=12,width=12)
ta = CrossTable(gtd[gtd$population == "CIcol",]$CNV, gtd[gtd$population == "CIcol",]$phenotype, fisher = T, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
pheatmap(t(ta$t), color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~CNV in CIcol\nFisher's exact test p=%.3E", ta$fisher.ts$p.value))

ta = CrossTable(gtd[gtd$population == "CIcol",]$estimated_n_ALT, gtd[gtd$population == "CIcol",]$phenotype, fisher = T, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
pheatmap(t(ta$t), color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~nALT in CIcol\nFisher's exact test p=%.3E", ta$fisher.ts$p.value))
dev.off()





#### GLM model with A65S and G280S data ####
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

### BIC PROCEDURE
# initial model (with all variables)
mod_tot = glm(phenotype ~ ., data = data, family = "binomial")
# stepwise removal using BIC criterion
mod_BIC = step(mod_tot, direction = "both", steps = 1e6, trace = F, k = log(nrow(data))) # k=log(num_obs) for BIC
mod_BIC_sum = summary(mod_BIC)
mod_tot_sum = summary(mod_tot)

# summary table BIC
mod_BIC_out = as.data.frame(mod_BIC_sum$coefficients)
colnames(mod_BIC_out) = c("BIC_coef","BIC_se", "BIC_z", "BIC_p_chisq")
mod_BIC_out$variant   = rownames(mod_BIC_out)
mod_BIC_out$BIC_OR    = 1/exp(mod_BIC_out$BIC_coef)

write.table(mod_BIC_out, "Fig3X_model_GLM_stepBIC_CIcol.csv" , quote = F, sep="\t")


# ### BACKWARD ELIMINATION
# source("../scripts_other/glmodelling.R")
# mod_BKE = glmodelling(input.table = data, list.of.markers = c("A65S","G280S", "CNV", "estimated_n_ALT"),
#                       rescolumn = "phenotype", glm.function = "glm", verbose = F)
# 
# mod_BKE_sum = summary(mod_BKE$final.model)
# mod_BKE_out = as.data.frame(mod_BKE_sum$coefficients)
# colnames(mod_BKE_out) = c("coef","se", "zscore", "p_chisq")
# mod_BKE_out$variant   = rownames(mod_BKE_out)
# mod_BKE_out$OR    = 1/exp(mod_BKE_out$coef)
# mod_BKE$final.sig

