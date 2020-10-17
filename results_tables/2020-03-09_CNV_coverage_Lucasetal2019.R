### Define input ####

library("scales")

# input files
getyf       = "ace1_genotypes.csv"
dupsf       = "ace1_duplications.csv"
pop1        = "alive"
pop2        = "dead"
sampf       = "freq_CIcol_CNV-ALTallele.csv"

graphics.off()

# ace location
ace_start = 3489213 / 1e6
ace_end   = 3493788 / 1e6
ace_mut = 3492074 / 1e6

# duplication breakpoints 
dup_start_major = 3436800 / 1e6
dup_end_major = 3639600 / 1e6
dup_start_minor = 3447900 / 1e6
dup_end_minor = 3518700 / 1e6

# add tagging variants
tags = c(3465693,3469441,3481632,3504796)

dup = read.table(dupsf, sep="\t", header = T)
met = read.csv(sampf,sep="\t",na.strings=c("","NA"))

# axis limts
pos_min_dataset = min(dup$Position) / 1e6
pos_max_dataset = max(dup$Position) / 1e6

# plot CN states in individual samples
# list of samples
list_samples = met$id

pdf(file="cnv_per_sample_coverage.pdf",height=4,width=8)

for (i in 1:length(list_samples)) {
  
  # find coverage for sample i
  dup_i = dup[dup$id == list_samples[i],]
  
  # where does its dup start and stop?
  dup_i_positions = dup_i[dup_i$CNV > 2, "Position"]
  if (length(dup_i_positions) == 0) {
    dup_i_pos_start = NA
    dup_i_pos_end = NA
  } else {
    dup_i_pos_start = min(dup_i_positions)
    dup_i_pos_end = max(dup_i_positions)
  }
  
  # plot
  plot(
    dup_i$Position/1e6,
    dup_i$Normalised_coverage,col="slategray", 
    xlim = c(pos_min_dataset, pos_max_dataset),  ylim=c(0,12),  xlab="Position", ylab="Normalised coverage", cex=0.7, las=1,
    main=sprintf("%s %s\n%i-%i", list_samples[i], met[i,"population"], dup_i_pos_start, dup_i_pos_end))
  lines(dup_i$Position/1e6,dup_i$CNV, col="blue")
  lines(c(ace_start,ace_end), c(0,0), lwd= 10, lend='butt', col="magenta")
  abline(h=2, lty=2, col="black")
  abline(v=c(dup_start_major, dup_end_major), lty=2, col="red")
  abline(v=c(dup_start_minor, dup_end_minor), lty=2, col="orange")
  
  abline(v=tags/1e6, lty=3, col="springgreen4")
  
  # for (tagi in tags) {
  #   
  #   pos_tagi_bool = tagi > dup_i$Position - 50 & tagi <= dup_i$Position + 300
  #   points(dup_i[pos_tagi_bool,]$Position/1e6, dup_i[pos_tagi_bool,]$Normalised_coverage, col="springgreen3", cex=0.7, las=1, pch=19)
  #   
  # }
  
}

dev.off()




# lot coverage per population, separating dup and nodup
list_populations = levels(met$population)

pdf(file="cnv_per_population_coverage.pdf",height=6,width=4)

popl = c("AOcol","BFcol","CIcol","GHcol","GNcol","BFgam","CMgam","FRgam","GAgam","GHgam","GNgam","GQgam","UGgam","GM","GW","KE")

sample_with_minor_dup = "AV0014-C"

for (pop in popl) {
  
  samples_in_pop_dup = as.character(met[met$population == pop & met$hasdup,"id"])
  samples_in_pop_nodup = as.character(met[met$population == pop & !met$hasdup,"id"])
  
  layout(c(1,2))
  
  # plot no dups
  plot(
    0,0, xlim = c(pos_min_dataset, pos_max_dataset),  ylim=c(0,12), las=1,
    xlab="Position", ylab="Normalised coverage", main=sprintf("%s no duplication", pop))
  lines(c(ace_start,ace_end), c(0,0), lwd= 10, lend='butt', col="magenta")
  abline(h=2, lty=2, col="black")
  abline(v=c(dup_start_major, dup_end_major), lty=2, col="red")
  abline(v=c(dup_start_minor, dup_end_minor), lty=2, col="orange")
  
  if (length(samples_in_pop_nodup) > 0) {
    for (i in 1:length(samples_in_pop_nodup)) {
      dup_i = dup[dup$id == samples_in_pop_nodup[i],]
      lines(dup_i$Position/1e6,dup_i$Normalised_coverage,col=alpha("blue", 0.2))
    }  
  }
  
  # plot with dups
  plot(
    0,0, xlim = c(pos_min_dataset, pos_max_dataset),  ylim=c(0,12), las=1,
    xlab="Position", ylab="Normalised coverage", main=sprintf("%s with duplication", pop))
  lines(c(ace_start,ace_end), c(0,0), lwd= 10, lend='butt', col="magenta")
  abline(h=2, lty=2, col="black")
  abline(v=c(dup_start_major, dup_end_major), lty=2, col="red")
  abline(v=c(dup_start_minor, dup_end_minor), lty=2, col="orange")
  
  if (length(samples_in_pop_dup) > 0) {
    for (i in 1:length(samples_in_pop_dup)) {
      dup_i = dup[dup$id == samples_in_pop_dup[i],]
      lines(dup_i$Position/1e6,dup_i$Normalised_coverage,col=alpha("blue", 0.2))
      # for (tagi in tags) {
      #   
      #   pos_tagi_bool = tagi > dup_i$Position - 100 & tagi <= dup_i$Position + 300
      #   points(dup_i[pos_tagi_bool,]$Position/1e6, dup_i[pos_tagi_bool,]$Normalised_coverage, col=alpha("springgreen3",0.6), cex=0.7, las=1, pch=19)
      #   
      # }
      
    } 
    
    # if (sample_with_minor_dup %in% samples_in_pop_dup) { 
    #   dup_i = dup[dup$id == sample_with_minor_dup,]
    #   lines(dup_i$Position/1e6,dup_i$Normalised_coverage,col="darkorange")
    # 
    # }
    
  }
  
  abline(v=tags/1e6, lty=3, col="springgreen4")
  
  
}

dev.off()

# same, with rolled averages
list_populations = levels(met$population)
library(zoo)

pdf(file="cnv_per_population_coverage_rolled.pdf",height=6,width=4)

popl = c("AOcol","BFcol","CIcol","GHcol","GNcol","BFgam","CMgam","FRgam","GAgam","GHgam","GNgam","GQgam","UGgam","GM","GW","KE")

sample_with_minor_dup = "AV0014-C"

for (pop in popl) {
  
  samples_in_pop_dup = as.character(met[met$population == pop & met$hasdup,"id"])
  samples_in_pop_nodup = as.character(met[met$population == pop & !met$hasdup,"id"])
  
  layout(c(1,2))
  
  # plot no dups
  plot(
    0,0, xlim = c(pos_min_dataset, pos_max_dataset),  ylim=c(0,12), las=1,
    xlab="Position", ylab="Normalised coverage", main=sprintf("%s no duplication", pop))
  lines(c(ace_start,ace_end), c(0,0), lwd= 10, lend='butt', col="magenta")
  abline(h=2, lty=2, col="black")
  abline(v=c(dup_start_major, dup_end_major), lty=2, col="red")
  abline(v=c(dup_start_minor, dup_end_minor), lty=2, col="orange")
  
  if (length(samples_in_pop_nodup) > 0) {
    for (i in 1:length(samples_in_pop_nodup)) {
      dup_i = dup[dup$id == samples_in_pop_nodup[i],]
      i_pos = zoo::rollapply(dup_i$Position/1e6, width=10, by=5, FUN=mean)
      i_nco = zoo::rollapply(dup_i$Normalised_coverage, width=10, by=5, FUN=mean)
      lines(i_pos,i_nco,col=alpha("blue", 0.2))
    }  
  }
  
  # plot with dups
  plot(
    0,0, xlim = c(pos_min_dataset, pos_max_dataset),  ylim=c(0,12), las=1,
    xlab="Position", ylab="Normalised coverage", main=sprintf("%s with duplication", pop))
  lines(c(ace_start,ace_end), c(0,0), lwd= 10, lend='butt', col="magenta")
  abline(h=2, lty=2, col="black")
  abline(v=c(dup_start_major, dup_end_major), lty=2, col="red")
  abline(v=c(dup_start_minor, dup_end_minor), lty=2, col="orange")
  
  if (length(samples_in_pop_dup) > 0) {
    for (i in 1:length(samples_in_pop_dup)) {
      dup_i = dup[dup$id == samples_in_pop_dup[i],]
      i_pos = zoo::rollapply(dup_i$Position/1e6, width=10, by=5, FUN=mean)
      i_nco = zoo::rollapply(dup_i$Normalised_coverage, width=10, by=5, FUN=mean)
      lines(i_pos,i_nco,col=alpha("blue", 0.2))
      # for (tagi in tags) {
      #   
      #   pos_tagi_bool = tagi > dup_i$Position - 100 & tagi <= dup_i$Position + 300
      #   points(dup_i[pos_tagi_bool,]$Position/1e6, dup_i[pos_tagi_bool,]$Normalised_coverage, col=alpha("springgreen3",0.6), cex=0.7, las=1, pch=19)
      #   
      # }
      
    } 
    
    # if (sample_with_minor_dup %in% samples_in_pop_dup) { 
    #   dup_i = dup[dup$id == sample_with_minor_dup,]
    #   lines(dup_i$Position/1e6,dup_i$Normalised_coverage,col="darkorange")
    # 
    # }
    
  }
  
  abline(v=tags/1e6, lty=3, col="springgreen4")
  
  
}

dev.off()

# coverage in duplicated specimens
pdf(file="cnv_per_duplication_coverage_boxplots.pdf",height=8,width=4)

samples_in_dup = as.character(met[met$hasdup,"id"])
samples_in_nodup = as.character(met[!met$hasdup,"id"])

layout(c(1,2))

# plot dups: normalised coverage
dup_i = dup[dup$id %in% samples_in_dup,c("id","Normalised_coverage","Position")]
dup_w = reshape(dup_i, idvar = "Position", direction = "wide", timevar = "id")
rownames(dup_w) = dup_w[,1]
dup_w = as.matrix(dup_w[,-1])

ace_window = which.min(abs(as.numeric(rownames(dup_w)) - ace_start*1e6))
ace_coverage = as.numeric(dup_w[ace_window,])
plot(ecdf(ace_coverage), col="magenta", xlim=c(0,8), verticals = T, pch = NA,ylab="CDF", xlab="Normalised coverage", main="Duplicated specimens")
abline(v=2, lty=2)
dat_coverage = matrix(nrow = length(ace_coverage), ncol = 5)
dat_coverage[,1]=ace_coverage
n=1
for (tag in tags) {
  n=n+1
  tag_window = which.min(abs(as.numeric(rownames(dup_w)) - tag))
  tag_coverage = dup_w[tag_window,]
  lines(ecdf(tag_coverage), col="blue", verticals = T, pch = NA)
  dat_coverage[,n] = tag_coverage
}
boxplot(dat_coverage, names = c("Ace1", tags), col="gray", ylab="Normalised coverage", ylim=c(0,8))
abline(h=2, lty=2)

# plot dups: raw coverage
dup_i = dup[dup$id %in% samples_in_dup,c("id","Counts","Position")]
dup_w = reshape(dup_i, idvar = "Position", direction = "wide", timevar = "id")
rownames(dup_w) = dup_w[,1]
dup_w = as.matrix(dup_w[,-1])

ace_window = which.min(abs(as.numeric(rownames(dup_w)) - ace_start*1e6))
ace_coverage = as.numeric(dup_w[ace_window,])
plot(ecdf(ace_coverage), col="magenta", xlim=c(0,600), verticals = T, pch = NA,ylab="CDF", xlab="Raw coverage", main="Duplicated specimens")
dat_coverage = matrix(nrow = length(ace_coverage), ncol = 5)
dat_coverage[,1]=ace_coverage
n=1
for (tag in tags) {
  n=n+1
  tag_window = which.min(abs(as.numeric(rownames(dup_w)) - tag))
  tag_coverage = dup_w[tag_window,]
  lines(ecdf(tag_coverage), col="blue", verticals = T, pch = NA)
  dat_coverage[,n] = tag_coverage
}
boxplot(dat_coverage, names = c("Ace1", tags), col="gray", ylab="Raw coverage", ylim=c(0,600))

#### INTERMISSION: PLOTS AND ECDFS FOR HAPLOTYPE SCORES ####

hsd = read.table("dupcoverage_stats.HaplotypeScore.csv", sep="\t", header = T)
ace_coverage = hsd[hsd$pos>ace_start*1e6 - 300 & hsd$pos < ace_end*1e6 + 300, "HaplotypeScore" ]
plot(ecdf(ace_coverage), col="magenta", xlim=c(0,20), ylim=c(0,1), verticals = T, pch = NA,ylab="CDF", xlab="Raw coverage", main="HS, entire dataset")
dat_coverage = list("Ace1" = ace_coverage)
n=1
for (tag in tags) {
  n=n+1
  tag_coverage = hsd[hsd$pos>tag - 300 & hsd$pos < tag + 300, "HaplotypeScore" ]
  lines(ecdf(tag_coverage), col="blue", verticals = T, pch = NA)
  dat_coverage[[n]] = tag_coverage
}
abline(v=2, lty=2, col="red")

boxplot(dat_coverage, names = c("Ace1", tags), col="gray", ylab="HaplotypeScore",
        sub="distribution of values in positions within region of interest")
abline(h=2, lty=2, col="red")

dev.off()


# coverage in duplicated specimens
pdf(file="cnv_per_duplication_coverage.pdf",height=10,width=8)
layout(c(1,2))


samples_in_dup = as.character(met[met$hasdup,"id"])
samples_in_nodup = as.character(met[!met$hasdup,"id"])

# plot with dups: norm coverage
plot(
  0,0, xlim = c(pos_min_dataset, pos_max_dataset),  ylim=c(0,12), las=1,
  xlab="Position", ylab="Normalised coverage", main=sprintf("%s with duplication", "all"))
lines(c(ace_start,ace_end), c(0,0), lwd= 10, lend='butt', col="magenta")
abline(h=2, lty=2)
abline(v=c(dup_start_major, dup_end_major), lty=2, col="red")
abline(v=c(dup_start_minor, dup_end_minor), lty=2, col="orange")


for (i in 1:length(samples_in_dup)) {
  dup_i = dup[dup$id == samples_in_dup[i],]
  i_pos = zoo::rollapply(dup_i$Position/1e6, width=10, by=2, FUN=mean)
  i_nco = zoo::rollapply(dup_i$Normalised_coverage, width=10, by=2, FUN=mean)
  lines(i_pos,i_nco,col=alpha("blue", 0.2))
} 
abline(v=tags/1e6, lty=3, col="springgreen4")


# plot with dups: raw coverage
plot(
  0,0, xlim = c(pos_min_dataset, pos_max_dataset),  ylim=c(0,600), las=1,
  xlab="Position", ylab="Raw coverage", main=sprintf("%s with duplication", "all"))
lines(c(ace_start,ace_end), c(0,0), lwd= 10, lend='butt', col="magenta")
abline(h=0, lty=2)
abline(v=c(dup_start_major, dup_end_major), lty=2, col="red")
abline(v=c(dup_start_minor, dup_end_minor), lty=2, col="orange")


for (i in 1:length(samples_in_dup)) {
  dup_i = dup[dup$id == samples_in_dup[i],]
  i_pos = zoo::rollapply(dup_i$Position/1e6, width=10, by=2, FUN=mean)
  i_nco = zoo::rollapply(dup_i$Counts, width=10, by=2, FUN=mean)
  lines(i_pos,i_nco,col=alpha("blue", 0.2))
} 
abline(v=tags/1e6, lty=3, col="springgreen4")


dev.off()


# library(vioplot)
# 
# for (pop in "CIcol") {
#   
#   samples_in_pop_dup = as.character(met[met$population == pop & met$hasdup,"id"])
#   samples_in_pop_nodup = as.character(met[met$population == pop & !met$hasdup,"id"])
#   
# 
#   for (tagi in tags)  {
#   vioplot(dup[dup$id %in% samples_in_pop_dup,]$Normalised_coverage,
#           dup[dup$id %in% samples_in_pop_dup & dup$Position >= tagi & dup$Position <= tagi+300,]$Normalised_coverage,
#           dup[dup$id %in% samples_in_pop_nodup,]$Normalised_coverage,
#           names = c("duplicated regions","tagging variants","non-duplicated regions"),las=2, col = "slategray3",)
#   }
# }

