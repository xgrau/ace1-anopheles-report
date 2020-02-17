### Define input ####

# input files
si          = "Anogam"                                  # Species name (requires Spi_long.annot.gff and Spi_gDNA.fasta)
dupsf       = "ace1_duplications.csv"
outcode     = "cida_dupl"
sampf       = "../metadata/samples.meta_phenotypes.txt"

# libraries
library(gmodels)
getmode = function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# define central loci
loci_p   = 3492074
loci_sta = 3484107
loci_end = 3495790


# load files
sam = read.csv(sampf,sep="\t",na.strings=c("","NA"))
dup = read.csv(dupsf,sep="\t",na.strings=c("","NA"))

# load genotypes for G280S
getyf  = "Genotype_3492074_G280S.csv"
gty    = read.csv(getyf,sep="\t",na.strings=c("","NA"))

# CNV state at loci
# dup coordinate that's closest to loci of interest:
dup_ix = base::which(abs(dup$Position-loci_p)==min(abs(dup$Position-loci_p)), useNames = F)
dup_i  = dup[dup_ix,]
# dup coordinates overlapping gene of interest
dup_co = base::which(dup$Position > loci_sta & dup$Position < loci_end)
dup_a  = dup[dup_co,]
dup_a  = aggregate(.~id,data=dup_i, getmode)
gtd = merge(gty, dup_a,by.x="ox_code", by.y="id")


print(paste(">>>>",getyf,"<<<<"))
# cross tables CI
print(CrossTable(gtd[gtd$population=="CIcol",]$genotype, 
                 gtd[gtd$population=="CIcol",]$phenotype,
                 fisher = T, prop.r = F, prop.c = F, prop.t = F ,expected = T))

print(CrossTable(gtd[gtd$population=="CIcol",]$CNV, 
                 gtd[gtd$population=="CIcol",]$genotype,
                 fisher = T, prop.r = F, prop.c = F, prop.t = F ,expected = T))



# load genotypes for A65S
getyf  = "Genotype_3489405_A65S.csv"
gty    = read.csv(getyf,sep="\t",na.strings=c("","NA"))

# CNV state at loci
# dup coordinate that's closest to loci of interest:
dup_ix = base::which(abs(dup$Position-loci_p)==min(abs(dup$Position-loci_p)), useNames = F)
dup_i  = dup[dup_ix,]
# dup coordinates overlapping gene of interest
dup_co = base::which(dup$Position > loci_sta & dup$Position < loci_end)
dup_a  = dup[dup_co,]
dup_a  = aggregate(.~id,data=dup_i, getmode)
gtd = merge(gty, dup_a,by.x="ox_code", by.y="id")


print(paste(">>>>",getyf,"<<<<"))
# cross tables CI
print(CrossTable(gtd[gtd$population=="CIcol",]$genotype, 
                 gtd[gtd$population=="CIcol",]$phenotype,
                 fisher = T, prop.r = F, prop.c = F, prop.t = F ,expected = T))

print(CrossTable(gtd[gtd$population=="CIcol",]$CNV, 
                 gtd[gtd$population=="CIcol",]$genotype,
                 fisher = T, prop.r = F, prop.c = F, prop.t = F ,expected = T))

message("FI")