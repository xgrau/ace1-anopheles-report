### Define input ####

# libraries
library(gmodels)
library(pheatmap)
source("../scripts_other/summarise_model_OR.R")

# input: genotyped/phenotyped samples from Biobank 
# genotypes:Dimitra Pipini, Emily Rippon, 
# phenotypes: Samuel Dadzie, Alexander Egyir-Yawson, John Essandoh, Joseph Chabi, Luc Djogbénou
dat_fn = "biobank_data_1080samples_9WA.csv"
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))


#### Genotype-phenotype associations: 119S ####

dat = read.table(dat_fn, header = T, sep="\t")
dat = dat[!is.na(dat$Species),]
dat = dat[!dat$Species == "gamcol",]
dat$population = paste(dat$Species, dat$Location, dat$Country)

pdf(file="Fig3E_gty_per_pop.pdf",height=4,width=4)
tat = CrossTable(dat$population, dat$Ace1_G119S)
pheatmap(tat$t, color = col.fun(20), breaks = seq(0,50,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("genotypes per pop"))
dev.off()

pop_order = rownames(tat$t)


# table genotype-phenotype associations
# loop with all tests
pdf(file="Fig3Esup_G280S.pdf",height=4,width=4)

for (pop in pop_order) {
  dai = dat[dat$population == pop,]
  tei = CrossTable(dai$phenotype, dai$Ace1_G119S, fisher=T, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
  
  # glm model
  if(ncol(tei$t) > 1 && nrow(tei$t) > 1 && !is.null(tei$fisher.ts$p.value)) {
    
    # GLM model
    # used to calculate ORs for the resistance genotypes
    data = dai[, c("phenotype","Ace1_G119S")]
    colnames(data) = c("phenotype","Ace1")
    mod_null = glm(phenotype ~ 1, data = data, family = "binomial")
    mod_full = glm(phenotype ~ ., data = data, family = "binomial")
    mod_signif = anova(mod_full, mod_null, test = 'Chisq')
    # pval string
    mod_pval_chisq_v_null = mod_signif$`Pr(>Chi)`[2]
    
    # report
    pheatmap(tei$t, color = col.fun(20), breaks = seq(0,30,length.out = 20),
             cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
             border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
             main=sprintf("phe~G280S\n%s\nFisher's exact test p = %.3E\nGLM p = %.3E\n%s",
                          pop, tei$fisher.ts$p.value,mod_pval_chisq_v_null, 
                          summarise_model_report_string(mod_full)))
    
    # save model table
    mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = paste(pop,"G280S"))
    write.table(file="Fig3Esup_G280S_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
    write.table(file="Fig3Esup_G280S_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
    
  }
}

for (pop in c("gam","col")) {
  dai = dat[dat$Species == pop,]
  tei = CrossTable(dai$phenotype, dai$Ace1_G119S, fisher=T, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
  
  # glm model
  if(ncol(tei$t) > 1 && nrow(tei$t) > 1 && !is.null(tei$fisher.ts$p.value)) {
    
    # GLM model
    # used to calculate ORs for the resistance genotypes
    data = dai[, c("phenotype","Ace1_G119S")]
    colnames(data) = c("phenotype","Ace1")
    mod_null = glm(phenotype ~ 1, data = data, family = "binomial")
    mod_full = glm(phenotype ~ ., data = data, family = "binomial")
    mod_signif = anova(mod_full, mod_null, test = 'Chisq')
    # pval string
    mod_pval_chisq_v_null = mod_signif$`Pr(>Chi)`[2]
    
    # report
    pheatmap(tei$t, color = col.fun(20), breaks = seq(0,30,length.out = 20),
             cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
             border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
             main=sprintf("phe~G280S\n%s\nFisher's exact test p = %.3E\nGLM p = %.3E\n%s",
                          pop, tei$fisher.ts$p.value,mod_pval_chisq_v_null, 
                          summarise_model_report_string(mod_full)))
    
    # save model table
    mod_tau = glm_tables(model=mod_full, null=mod_null, model_name = paste(pop,"G280S"))
    write.table(file="Fig3Esup_G280S_models.csv", t(mod_tau$model_table), quote=FALSE, sep="\t", col.names=FALSE, append = T)
    write.table(file="Fig3Esup_G280S_models.csv", mod_tau$variable_table, quote=FALSE, sep="\t", row.names=FALSE, append = T)
    
  }
}

dev.off()


#### Plot frequency of each species in each location ####

pdf(file="Fig3Esup_populations.pdf",height=12,width=12)

tei = CrossTable(dat$Ace1_G119S, dat$Species, fisher=F, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
pheatmap(tei$t, color = col.fun(20), breaks = seq(0,300,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("genotype~species"))

tei = CrossTable(dat$Location, dat$Species, fisher=F, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
pheatmap(tei$t, color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("location~species"))

tei = CrossTable(dat$Country, dat$Species, fisher=F, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
pheatmap(tei$t, color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("country~species"))

tei = CrossTable(dat$Country, dat$Location, fisher=F, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
pheatmap(tei$t, color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("country~location"))

dev.off()

message("FI 119S")

