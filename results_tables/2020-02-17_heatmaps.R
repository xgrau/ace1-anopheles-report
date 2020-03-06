### Define input ####

# input files
getyf       = "ace1_genotypes.csv"
dupsf       = "ace1_duplications.csv"
pop1        = "alive"
pop2        = "dead"
sampf       = "../metadata/samples.meta_phenotypes.txt"

graphics.off()

# libraries
library(gmodels)
library(pheatmap)
library(Rfast)
library(stringr)
getmode = function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



popl = c("AOcol","BFcol","CIcol","GHcol","GNcol","BFgam","CMgam","FRgam","GAgam","GHgam","GNgam","GQgam","UGgam","GM","GW","KE")
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))


##### Frequency of non-syn variants in Ace1 #####
# per population
# load allele freqs
allel   = read.table("../results_hap_analysis/AlleleFq_tab.csv",header = T,sep="\t")
# allel   = read.table("input_Ace1due.AlleleFq_tab.csv",header = T,sep="\t")
popl_fq = paste(popl,"fqmin",sep="_")

# is allele nonsyn?
aa_ALT = sapply(stringr::str_split(allel$PEP_eff, pattern="[0-9]+"), "[", 2)
aa_REF = sapply(stringr::str_split(allel$PEP_eff, pattern="[0-9]+"), "[", 1)
aa_REF = sapply(stringr::str_split(aa_REF, pattern="\\."), "[", 2)
allel$REF_aa = aa_REF
allel$ALT_aa = aa_ALT
allel$is_nonsyn = allel$REF_aa != allel$ALT_aa

# find nonsyn in Ace1
allel_f = allel[allel$gene_eff == "AGAP001356" & allel$is_nonsyn & !is.na(allel$is_nonsyn),]
allel_f[,popl_fq][is.na(allel_f[,popl_fq])] = 0
allel_f = allel_f[apply(allel_f[,popl_fq], 1, FUN=max) > 0.05,]
rownames(allel_f) = paste(allel_f$chr,":",allel_f$POS," ",gsub("n.","",allel_f$CDS_eff)," ",gsub("p.","",allel_f$PEP_eff)," | reverse",rowMeans(allel_f[,popl_fq])>0.5,sep="")

# canvia major a minor a un allel concret
allel_f[rowMeans(allel_f[,popl_fq])>0.5,popl_fq] = 1-allel_f[rowMeans(allel_f[,popl_fq])>0.5,popl_fq]


##### Calculate number of ALT alleles per specimen #####
# define central loci
loci_p   = 3492074
loci_sta = 3484107
loci_end = 3495790

# load files
sam = read.csv(sampf,sep="\t",na.strings=c("","NA"))
dup = read.csv(dupsf,sep="\t",na.strings=c("","NA"))
gty = read.csv(getyf,sep="\t",na.strings=c("","NA"))
colnames(gty)[colnames(gty) == "CNstate"] = "CNstate_OLD"

# define populations
pop1_l = as.character(sam[sam$population=="CIcol" & sam$phenotype == pop1, "ox_code"])
pop2_l = as.character(sam[sam$population=="CIcol" & sam$phenotype == pop2, "ox_code"])
popo_l = as.character(sam[sam$population!="CIcol" , "ox_code"])

# CNV state at loci
# dup coordinate that's closest to loci of interest:
dup_ix = base::which(abs(dup$Position-loci_p)==min(abs(dup$Position-loci_p)), useNames = F)
dup_i  = dup[dup_ix,]
# dup coordinates overlapping gene of interest
dup_co = base::which(dup$Position > loci_sta & dup$Position < loci_end)
dup_a  = dup[dup_co,]
dup_a  = aggregate(.~id,data=dup_i, getmode)


# merge genotype with CNV state
gty$id = rownames(gty)
gtd = merge(gty, dup_a,by="id")
gtd = merge(gtd, sam, by.x="id", by.y="ox_code")

# estimated number of alternate alleles
gtd$estimated_n_ALT = apply(gtd,1, FUN = function(x) {
  cnv   = as.numeric(x["CNV"])
  fra   = as.numeric(x["proportion"])
  cnv_v = seq(0,cnv, length.out = cnv+1)/cnv
  cnv_i = which(abs(cnv_v-fra)==min(abs(cnv_v-fra)))
  est   = cnv_v[cnv_i]*cnv
})


##### heatmaps with genotype and duplication frequencies #####

popl = c("AOcol","BFcol","BFgam","CIcol","CMgam","FRgam","GAgam","GHcol","GHgam","GM","GNcol","GNgam","GQgam","GW","KE","UGgam")
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))


# plot duplications Fig1A
pdf(file="Fig1A_freqs_nonsyn.pdf",height=12,width=12)
pheatmap(t(allel_f[,popl_fq]), color = col.fun(20), breaks = seq(0,1,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "red",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,
         main=paste("nonsyn freqs"))
dev.off()

# 119S genotypes per population Fig1B
pdf(file="Fig1B_freqs_119Sgty_per_pop.pdf",height=12,width=12)
ta = CrossTable(gtd$population, gtd$genotype,fisher = F, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
ta$t = cbind(ta$t,rowSums(ta$t))
pheatmap(ta$t, color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=paste("genotype~pop"))
dev.off()

# number of Ace1 copies per population Fig2A
pdf(file="Fig2A_freqs_CNV_per_pop.pdf",height=12,width=12)
ta = CrossTable(gtd$population, gtd$CNV, fisher = F, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
ta$t = cbind(ta$t,rowSums(ta$t))
pheatmap(ta$t, color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=paste("CNV~pop"))
dev.off()

# number of Ace1 copies per 119S genotype Fig2B
pdf(file="Fig2B_freqs_CNV_per_gty.pdf",height=12,width=12)
ta = CrossTable(gtd$genotype, gtd$CNV, fisher = F, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
ta$t = cbind(ta$t,rowSums(ta$t))
pheatmap(ta$t, color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=paste("genotype~CNV"))
dev.off()


# fig10 summary duplications
pdf(file="Fig2C_summarydups.pdf",height=12,width=12)
ta = CrossTable(gtd$estimated_n_ALT, gtd$CNV, fisher = F, prop.r = F, prop.c = F, prop.t = F ,prop.chisq = F)
ta$t = cbind(ta$t,rowSums(ta$t))
pheatmap(ta$t, color = col.fun(20), breaks = seq(0,10,length.out = 20), 
         cellwidth = 18, cellheight = 12, na_col = "dodgerblue4",number_color = "aliceblue",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,number_format = "%i",
         main=paste("CNV~estimated"))
dev.off()



#### frequency of ALT alleles per sample (fig2CD) ####
write.table(gtd, "Fig3_CIcol_CNV-ALTallele.csv", sep="\t",quote = F, row.names = F)
pdf(file="Fig2CD_ALTallele.pdf",height=8,width=12)
layout(matrix(1:6,nrow=2))


# plot duplication data: fraction
plot(gtd[gtd$id %in% popo_l & gtd$genotype == "wt/wt","proportion"],gtd[gtd$id %in% popo_l  & gtd$genotype == "wt/wt","CNV"], col="slategray2",
     xlab="fraction ALT alleles",ylab="# copies",main="ALT allele fraction",xlim=c(0,1),ylim=c(2,10))
points(gtd[gtd$id %in% popo_l & gtd$genotype != "wt/wt","proportion"],gtd[gtd$id %in% popo_l  & gtd$genotype != "wt/wt","CNV"],col="slategray4")
points(gtd[gtd$id %in% pop2_l,"proportion"],gtd[gtd$id %in% pop2_l,"CNV"],col="magenta3")
points(gtd[gtd$id %in% pop1_l,"proportion"],gtd[gtd$id %in% pop1_l,"CNV"],col="springgreen3")
legend("topleft",legend=c("alive","dead","other"),col=c("springgreen3","magenta3","slategray"),pch=1,bty="n")

# plot duplication data: counts
plot(gtd[gtd$id %in% popo_l,"A"],gtd[gtd$id %in% popo_l,"CNV"], col="slategray",
     xlab="# ALT alleles",ylab="# copies",main="ALT allele count")
points(gtd[gtd$id %in% pop2_l,"A"],gtd[gtd$id %in% pop2_l,"CNV"],col="magenta3")
points(gtd[gtd$id %in% pop1_l,"A"],gtd[gtd$id %in% pop1_l,"CNV"],col="springgreen3")
legend("topleft",legend=c("alive","dead","other"),col=c("springgreen3","magenta3","slategray"),pch=1,bty="n")

# plot alternate fraction vs phenotype 
popx_f_ecdf = seq(0,max(gtd$proportion),length.out = 1000)
pop1_f_ecdf = ecdf(gtd[gtd$id %in% pop1_l,"proportion"])
pop2_f_ecdf = ecdf(gtd[gtd$id %in% pop2_l,"proportion"])
popo_f_ecdf = ecdf(gtd[gtd$id %in% popo_l,"proportion"])
kst         = ks.test(gtd[gtd$id %in% pop1_l,"proportion"],gtd[gtd$id %in% pop2_l,"proportion"])
plot(popx_f_ecdf,popo_f_ecdf(popx_f_ecdf),col="slategray",type="l",ylim=c(0,1),
     xlab="fraction ALT alleles",ylab="Fraction",main="CDF ALT allele fraction",
     sub=paste(pop1,"~",pop2,"KS-2s p =",signif(kst$p.value,3)))
lines(popx_f_ecdf,pop1_f_ecdf(popx_f_ecdf),col="springgreen3",type="l")
lines(popx_f_ecdf,pop2_f_ecdf(popx_f_ecdf),col="magenta3",type="l")
legend("bottomright",legend=c("alive","dead","other"),col=c("springgreen3","magenta3","slategray"),lty=1,bty="n")

# plot alternate count vs phenotype 
popx_f_ecdf = seq(0,max(gtd$A),length.out = 1000)
pop1_f_ecdf = ecdf(gtd[gtd$id %in% pop1_l,"A"])
pop2_f_ecdf = ecdf(gtd[gtd$id %in% pop2_l,"A"])
popo_f_ecdf = ecdf(gtd[gtd$id %in% popo_l,"A"])
kst         = ks.test(gtd[gtd$id %in% pop1_l,"A"],gtd[gtd$id %in% pop2_l,"A"])
plot(popx_f_ecdf,popo_f_ecdf(popx_f_ecdf),col="slategray",type="l",ylim=c(0,1),
     xlab="# ALT alleles",ylab="Fraction",main="CDF ALT allele count",
     sub=paste(pop1,"~",pop2,"KS-2s p =",signif(kst$p.value,3)))
lines(popx_f_ecdf,pop1_f_ecdf(popx_f_ecdf),col="springgreen3",type="l")
lines(popx_f_ecdf,pop2_f_ecdf(popx_f_ecdf),col="magenta3",type="l")
legend("bottomright",legend=c("alive","dead","other"),col=c("springgreen3","magenta3","slategray"),lty=1,bty="n")

# plot density alt fraction
kst = ks.test(gtd[gtd$id %in% pop1_l,"proportion"],gtd[gtd$id %in% pop2_l,"proportion"])
plot(density(gtd[gtd$id %in% popo_l,"proportion"]),col="slategray",type="l",
     xlab="fraction ALT alleles",ylab="Fraction",main="CDF ALT allele fraction",xlim=c(0,1),
     sub=paste(pop1,"~",pop2,"KS-2s p =",signif(kst$p.value,3)))
lines(density(gtd[gtd$id %in% pop1_l,"proportion"]),col="springgreen3",type="l")
lines(density(gtd[gtd$id %in% pop2_l,"proportion"]),col="magenta3",type="l")
legend("topright",legend=c("alive","dead","other"),col=c("springgreen3","magenta3","slategray"),lty=1,bty="n")

# plot density alt count
kst = ks.test(gtd[gtd$id %in% pop1_l,"A"],gtd[gtd$id %in% pop2_l,"A"])
plot(density(gtd[gtd$id %in% popo_l,"A"]),col="slategray",type="l",
     xlab="# ALT alleles",ylab="Fraction",main="CDF ALT allele count",
     sub=paste(pop1,"~",pop2,"KS-2s p =",signif(kst$p.value,3)))
lines(density(gtd[gtd$id %in% pop1_l,"A"]),col="springgreen3",type="l")
lines(density(gtd[gtd$id %in% pop2_l,"A"]),col="magenta3",type="l")
legend("topright",legend=c("alive","dead","other"),col=c("springgreen3","magenta3","slategray"),lty=1,bty="n")

dev.off()



##### Frequency of minor duplications ####
# Unused; data from Lucas et al 2019 

# load duplications
dupf    = read.table("eric_dups_ace1.csv",header = T)
dupf_f  = dupf[,popl]
dupf_id = paste(dupf$chr,":",dupf$start,"-",dupf$end," ",dupf$CNV_id,sep="")
rownames(dupf_f) = dupf_id

# plot duplications
pdf(file=paste("FigXX_minor_duplications.pdf",sep=""),height=12,width=12)
pheatmap(t(dupf_f), color = col.fun(20), breaks = seq(0,1,length.out = 20), 
         cellwidth = 18, cellheight = 12,
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = T,
         main=paste("dup freqs"))
dev.off()


print("FI!")
