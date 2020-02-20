library(Biostrings)
library(stringr)

options(scipen = 999)

cat('Loading kmer table.\n')
sig.table <- read.table('full_sig_kmer_table.csv', header = T, row.names = 1)
sample.names <- colnames(sig.table)[3:ncol(sig.table)]

cat('Loading Ag1000G metadata\n')
meta <- read.table('../metadata/samples.meta_phenotypes.txt', sep = '\t', header = T, row.names = 1, quote = '')
CI.meta <- subset(meta, population == 'CIcol')
CI.alive <- grepl('aPM', CI.meta$src_code)
names(CI.alive) <- sub('-', '.', row.names(CI.meta))
alive <- names(CI.alive)[CI.alive]
dead <- names(CI.alive)[!CI.alive]

# So now we have associated the kmer names with their distribution patterns. We can now use the kmer 
# dictionary to work out which kmers were assembled and then aligned to non-Ace1 regions. 
kmer.dict <- read.table('full_sig_kmer_dictionary.txt', row.names = 1, stringsAsFactors = F, col.names = c('Sequence', 'kmers'))
kmer.dict$num.kmer <- unlist(lapply(strsplit(kmer.dict$kmers, '\\.'), length))
# Identify the sequences composed of at least 2 kmers
filtered.sequences <- rownames(kmer.dict)[kmer.dict$num.kmer > 1]

# And load the results of the alignment
alignment.results <- read.table('full_sig_kmers.csv', sep = '\t', header = T, stringsAsFactors = F)
# Filter out all the alignments that we are not interested in
filtered.alignments <- subset(alignment.results, Kmer %in% filtered.sequences)
filtered.alignments <- subset(filtered.alignments, !is.na(Chrom))
filtered.alignments <- subset(filtered.alignments, !(Chrom %in% c('Y_unplaced', 'UNKN')))
filtered.alignments <- subset(filtered.alignments, Supplementary == 'False')


# Find the alignments that were not to Ace-1, and remove those that aligned nowhere or to the unknown
# chromosome or Y chromosomes, or that are supplementary alignments.
non.Ace1.alignments <- subset(filtered.alignments, In_Dup1 != 'True')
# Use the dictionary to identify the full list of kmers that went into assembling these kmers.
non.Ace1.kmers.unsplit <- kmer.dict[non.Ace1.alignments$Kmer, , drop = F]
rownames(non.Ace1.kmers.unsplit) <- non.Ace1.alignments$Kmer
# Split up each of those strings into their constituent kmer names
non.Ace1.kmers <- strsplit(non.Ace1.kmers.unsplit[[1]], '\\.')
names(non.Ace1.kmers) <- non.Ace1.alignments$Kmer
# Create a table where each row is an assembled kmer, and contains the mean normalised kmer counts for
# each kmer associated with that assembled kmer. 
non.Ace1.kmers.mean.counts <- do.call(rbind, lapply(non.Ace1.kmers, function(x) apply(sig.table[x, 3:ncol(sig.table)], 2, mean)))

# Now do the same for kmers not in the Dup1 region
Ace1.alignments <- subset(filtered.alignments, In_Dup1 == 'True')
Ace1.kmers.unsplit <- kmer.dict[Ace1.alignments$Kmer, , drop = F]
rownames(Ace1.kmers.unsplit) <- Ace1.alignments$Kmer
Ace1.kmers <- strsplit(Ace1.kmers.unsplit[[1]], '\\.')
names(Ace1.kmers) <- Ace1.alignments$Kmer
Ace1.kmers.mean.counts <- do.call(rbind, lapply(Ace1.kmers, function(x) apply(sig.table[x, 3:ncol(sig.table)], 2, mean)))

# Load Xavi's Ace-1 coverage data. 
Ace1.coverage <- read.table('../results_tables/Fig3_CIcol_CNV-ALTallele.csv', header = T, row.names =1)
rownames(Ace1.coverage) <- sub('-', '.', rownames(Ace1.coverage))
CI.Ace1.CNV <- Ace1.coverage[sample.names, 'CNV']

# Write a function that will calculate the residuals of kmer frequency within each of the dead and alive categories
# then get a correlation coefficient within each of those. USing tapply to calculate these residuals will change 
# the order of the values, but will do so in the same way for the CNV and the kmer frequency.
Ace1.CNV.residuals <- unlist(tapply(CI.Ace1.CNV, CI.alive, function(x) x-mean(x)))
residual.correlation <- function(x){
  if (any(is.na(x)))
    return(NA)
  x.residual <- unlist(tapply(x, CI.alive, function(x) x-mean(x)))
  cor.test(Ace1.CNV.residuals, x.residual, method = 'pearson')$p.value
}
residual.correlation.estimate <- function(x){
  if (any(is.na(x)))
    return(NA)
  x.residual <- unlist(tapply(x, CI.alive, function(x) x-mean(x)))
  cor.test(Ace1.CNV.residuals, x.residual, method = 'pearson')$estimate
}

# Calculate the correlation of each assembled sequence with Ace-1 copy number
non.Ace1.CNV.correlation <- apply(non.Ace1.kmers.mean.counts, 1, residual.correlation)
Ace1.CNV.correlation <- apply(Ace1.kmers.mean.counts, 1, residual.correlation)
non.Ace1.CNV.estimate <- apply(non.Ace1.kmers.mean.counts, 1, residual.correlation.estimate)
Ace1.CNV.estimate <- apply(Ace1.kmers.mean.counts, 1, residual.correlation.estimate)

cat('\n', sum(Ace1.CNV.correlation < 0.05), ' out of ', length(Ace1.CNV.correlation), ' kmers mapping to the Ace-1 ',
    'CNV were correlated with CNV copy number.\n', sep = '')
cat(sum(non.Ace1.CNV.correlation < 0.05), ' out of ', length(non.Ace1.CNV.correlation), ' kmers NOT mapping to the ',
    'Ace1 CNV were correlated with CNV copy number.\n\n', sep = '')

# There are two assembled sequences that are not correlated with Ace1 copy number
non.Ace1.kmers <- names(non.Ace1.CNV.correlation)[non.Ace1.CNV.correlation > 0.05]
non.Ace1.table <- subset(filtered.alignments, Kmer %in% non.Ace1.kmers)

print(non.Ace1.table)


# output: dataframe with Pearson correlation and one plot
out_df = data.frame(
  row.names = names(c(non.Ace1.CNV.estimate, Ace1.CNV.estimate)),
  pearson_estimate = c(non.Ace1.CNV.estimate, Ace1.CNV.estimate),
  pearson_pval = c(non.Ace1.CNV.correlation, Ace1.CNV.correlation),
  is_in_Ace1 = names(c(non.Ace1.CNV.estimate, Ace1.CNV.estimate)) %in% names(Ace1.CNV.estimate)
)

write.table(out_df, file = "full_sig_kmer_table_Ace1_Pearson.csv", quote = F, sep="\t", row.names = T)
pdf("full_sig_kmer_table_Ace1_Pearson.pdf", width = 4, height = 4)
order_estimate = order(out_df$pearson_estimate)
plot(out_df$pearson_estimate[order_estimate], ylim=c(-1,1), col=as.factor(out_df$is_in_Ace1)[order_estimate], xlab="k-mers", ylab="Pearson correlation", cex=0.8)
legend("bottomright", legend = c("not in Ace1","in Ace1"), col = c("black", "red"), pch = 1, bty="n")
dev.off()


# more output: frequenies of kmers per sample in a heatmap
kmer_order = rownames(out_df)
kmer_order_table = sig.table[match(kmer_order, rownames(sig.table)),]
kmer_order_table <- subset(kmer_order_table, select = -c(V1,V2))

sample_order = rownames(CI.meta[order(CI.meta$phenotype),])
sample_order = gsub("-",".", sample_order)
kmer_order_table = kmer_order_table[,match(sample_order, colnames(kmer_order_table))]


library(pheatmap)
col.fun = colorRampPalette(interpolate="l",c("aliceblue","deepskyblue","dodgerblue4"))
cod.fun = colorRampPalette(c("firebrick4","orangered", "floralwhite", "deepskyblue","dodgerblue4"))
pdf(file="full_sig_kmer_table_Ace1_heatmap.pdf",height=4,width=2)
pheatmap(kmer_order_table, color = col.fun(20), breaks = seq(0,1,length.out = 20), 
         na_col = "dodgerblue4",
         border_color = "white", cluster_cols=F, cluster_rows=F,display_numbers = F, labels_row = "k-mers", labels_col = "Samples",
         main=paste("kmer freq"))
dev.off()



