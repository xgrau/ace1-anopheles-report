### Define input ####

# libraries
library(gmodels)


#### Genotype-phenotype associations: 119S ####

# input: genotyped/phenotyped samples from Biobank (D. Pipini)
dat = read.table("input_genotype-phenotype_2020-02-17.all_data.csv", header = T, sep="\t")
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))


# table genotype-phenotype associations: CNVs
pdf(file="FigBiobank_phe-gty.pdf",height=12,width=12)
# loop with all tests
for (spi in c("col","gam")) {
  for (loi in unique(dat$Location)) {
    dai = dat[dat$Location == loi & dat$Species == spi,]
    country = as.vector(unique(dai$Country))
    if (nrow(dai)>0 && sum(table(dai$phenotype)>1) > 1)  {
      tei = CrossTable(dai$phenotype, dai$Ace1_G119S, fisher=T, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
      pheatmap(tei$t, color = col.fun(20), breaks = seq(0,10,length.out = 20), 
               cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
               border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
               main=sprintf("phe~119S in %s %s %s\nFisher's exact test p=%.3E",
                            spi, country, loi, tei$fisher.ts$p.value))
    }
  }
}
dev.off()



#### Plot frequency of each species in each location ####

pdf(file="FigBiobank_populations.pdf",height=12,width=12)

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

