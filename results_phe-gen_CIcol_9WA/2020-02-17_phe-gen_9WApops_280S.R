### Define input ####

# libraries
library(gmodels)
library(pheatmap)


#### Genotype-phenotype associations: 119S ####

# input: genotyped/phenotyped samples from Biobank (D. Pipini)
dat = read.table("input_genotype-phenotype_2020-02-17.all_data.csv", header = T, sep="\t")
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))


dat$population = paste(dat$Species, dat$Location, dat$Country)

daf = dat[dat$Species %in% c("col","gam") && dat$phenotype,]

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



#### Genotype-phenotype associations: 119S ####
#### TODO: SAME WITH DUPLICATIONS


# pdf(file="FigY_phe-gty_CNVs.pdf",height=12,width=12)
# dev.off()

