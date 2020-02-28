### Define input ####

# libraries
library(gmodels)
library(pheatmap)


# input: genotyped/phenotyped samples from Biobank 
# genotypes:Dimitra Pipini, Emily Rippon, 
# phenotypes: Samuel Dadzie, Alexander Egyir-Yawson, John Essandoh, Joseph Chabi, Luc Djogbénou
dat_fn = "biobank_samples_gen-phe.csv"
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))


#### Genotype-phenotype associations: 119S ####

dat = read.table(dat_fn, header = T, sep="\t")
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
pdf(file="Fig3-sup_phe-gty_119S.pdf",height=4,width=4)

for (pop in pop_order) {
  dai = dat[dat$population == pop,]
  tei = CrossTable(dai$phenotype, dai$Ace1_G119S, fisher=T, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
  
  if(ncol(tei$t) > 1 && nrow(tei$t) > 1 && !is.null(tei$fisher.ts$p.value)) {
    pheatmap(tei$t, color = col.fun(20), breaks = seq(0,30,length.out = 20),
             cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
             border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
             main=sprintf("phe~119S in %s\nFisher's exact test p=%.3E",
                          pop, tei$fisher.ts$p.value))
  }
}
dev.off()


#### Plot frequency of each species in each location ####

pdf(file="Fig3-sup_populations.pdf",height=12,width=12)

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

#### Genotype-phenotype associations: CNVs ####

# load CNV data for a subset of samples
dac = read.table("biobank_samples_gen-phe_subset_CNV.csv", sep="\t", header = T)

# convert CNV estimate (number of copies relative to Kisumu, where wt copy number = 2) to CNV as in Lucas et al
dac$CNV = dac$CNV_estimate * 2 
dac = merge(dat,dac,by.x="Serial", by.y = "serial.no")


# populations to test
pop_order = sort(unique(dac$population))

for (pop in pop_order) {
  
  # variables for model
  data = dac[dac$population == pop, c("phenotype","ratio_FAM_HEX","CNV")]
  
  # full model
  mod_pop = glm(phenotype ~ CNV+ratio_FAM_HEX, data = data, family = "binomial")
  mod_BKE = glmodelling(input.table = data, list.of.markers = c("ratio_FAM_HEX","CNV"), 
                        rescolumn = "phenotype", glm.function = "glm", verbose = F)
  mod_null = glm(phenotype ~ 1, data = data, family = "binomial")
  
  # full v null
  mod_pop_signif = anova(mod_BKE$final.model, mod_null, test = 'Chisq')
  mod_pop_pval_chisq_v_null = mod_pop_signif$`Pr(>Chi)`[2]
  mod_pop_variables = paste(rownames(mod_BKE$final.sig), collapse = ",")
  
  # report
  print(paste(pop,nrow(data), mod_pop_pval_chisq_v_null, mod_pop_variables, 
              signif(mean(data$CNV),3), signif(sd(data$CNV)/sqrt(nrow(data)),3) ) )
  print(table(data$phenotype))
  
  
  layout(c(1,2), widths = 2)
  plot(data$phenotype, data$CNV, main=pop, ylab="CNV")
  plot(data$phenotype, data$ratio_FAM_HEX, main=pop, ylab="ratio")
  
  # # spearman correlation
  # vec_cnv = dac[dac$population == pop,"CNV"]
  # vec_phe = as.numeric(as.vector(dac[dac$population == pop,"phenotype"]) == "Alive")
  # cor_cnv_phe = cor.test(vec_cnv , vec_phe, method = "spearman")
  # plot(vec_cnv, vec_phe, main=pop, col="blue")
  # print(paste(pop,length(vec_cnv), cor_cnv_phe$p.value, cor_cnv_phe$estimate))
  
}

# do the same thing but with a larger model where species is a 
# random effect, with glmer function (lme4)

library(lme4)
source("../scripts_other/glmodelling.R")

# variables for model
data = dac[,c("phenotype","ratio_FAM_HEX","CNV","Species","Location")]

# null model
mod_null = glmer(phenotype ~ 1|Species, data = data, family = "binomial")

# initial model (with all variables)
mod_tot = glmer(phenotype ~ CNV+ratio_FAM_HEX|Species, data = data, family = "binomial")

mod_tot_signif = anova(mod_tot, mod_null, test = 'Chisq')
mod_tot_pval_chisq_v_null = mod_tot_signif$`Pr(>Chisq)`[2]

# now do the same thing within species, using locatio as a random effect

for (spi in c("col","gam")){
  
  # variables for model
  data = dac[dac$Species == spi, c("phenotype","ratio_FAM_HEX","CNV","Location")]
  
  # RANDOM EFFECTS
  mod_spi = glmer(phenotype ~ CNV+ratio_FAM_HEX|Location, data = data, family = "binomial")
  mod_BKE = glmodelling(input.table = data, list.of.markers = c("ratio_FAM_HEX","CNV"), control.for = c("Location"),
                        rescolumn = "phenotype", glm.function = "glmer", verbose = F)
  mod_null = glmer(phenotype ~ 1|Location, data = data, family = "binomial")
  mod_spi_signif = anova(mod_BKE$final.model, mod_null, test = 'Chisq')
  mod_spi_variables = paste(rownames(mod_BKE$final.sig), collapse = ",")
  mod_spi_pval_chisq_v_null = mod_spi_signif$`Pr(>Chisq)`[2]
  
  # # NO RANDOM 
  # mod_spi = glm(phenotype ~ CNV+ratio_FAM_HEX, data = data, family = "binomial")
  # mod_null = glm(phenotype ~ 1, data = data, family = "binomial")
  # mod_spi_signif = anova(mod_spi, mod_null, test = 'Chisq')
  # mod_spi_pval_chisq_v_null = mod_spi_signif$`Pr(>Chi)`[2]
  
  
  # report
  print(paste(spi,nrow(data), mod_spi_pval_chisq_v_null, mod_spi_variables,
              signif(mean(data$CNV),3), signif(sd(data$CNV)/sqrt(nrow(data)),3) ) )
  print(table(data$phenotype))
  
}


message("FI")
