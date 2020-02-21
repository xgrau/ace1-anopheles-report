### Define input ####

# libraries
library(gmodels)
library(pheatmap)


#### Genotype-phenotype associations: 119S ####

# input: genotyped/phenotyped samples from Biobank (D. Pipini)
dat = read.table("input_genotype-phenotype_2020-02-17.all_data.csv", header = T, sep="\t")
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))


# table genotype-phenotype associations: CNVs
pdf(file="FigY_phe-gty_119S.pdf",height=4,width=4)

# loop with all tests
for (spi in c("col","gam")) {
  for (loi in unique(dat$Location)) {
    dai = dat[dat$Location == loi & dat$Species == spi,]
    country = as.vector(unique(dai$Country))
    if (nrow(dai)>0 && sum(table(dai$phenotype)>1) > 1)  {
      tei = CrossTable(dai$phenotype, dai$Ace1_G119S, fisher=T, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
      pheatmap(tei$t, color = col.fun(20), breaks = seq(0,30,length.out = 20), 
               cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
               border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
               main=sprintf("phe~119S in %s %s %s\nFisher's exact test p=%.3E\nResistant GS = %.3f",
                            spi, country, loi, tei$fisher.ts$p.value, tei$prop.col[,"GS"]["Alive"]))
    }
  }
}

# add table with ALL col
tei = CrossTable(dat[dat$Species == "col",]$phenotype, dat[dat$Species == "col",]$Ace1_G119S, fisher=T, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
pheatmap(tei$t, color = col.fun(20), breaks = seq(0,200,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~119S in ALL col\nFisher's exact test p=%.3E\nResistant GS = %.3f\nResistant SS = %.3f",
                      tei$fisher.ts$p.value, tei$prop.col[,"GS"]["Alive"], tei$prop.col[,"SS"]["Alive"]))
# add table with ALL gam
tei = CrossTable(dat[dat$Species == "gam",]$phenotype, dat[dat$Species == "gam",]$Ace1_G119S, fisher=T, prop.r = F, prop.c = F, prop.t = F, expected = T, missing.include = T)
pheatmap(tei$t, color = col.fun(20), breaks = seq(0,200,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=sprintf("phe~119S in ALL gam\nFisher's exact test p=%.3E\nResistant GS = %.3f\nResistant SS = %.3f",
                      tei$fisher.ts$p.value, tei$prop.col[,"GS"]["Alive"], tei$prop.col[,"SS"]["Alive"]))

dev.off()



#### Plot frequency of each species in each location ####

pdf(file="FigY_populations.pdf",height=12,width=12)

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

