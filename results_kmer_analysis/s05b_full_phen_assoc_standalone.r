library(rhdf5)
library(data.table)
library(stringr)
library(parallel)
library(fdrtool)
library(bit64)
options(scipen = 999)

start.time <- proc.time()
cat('Process time at start:\n\n')
print(start.time)

# Load the total number of kmers per sample
cat('Loading sample-wise kmer totals.\n')
total.kmers <- read.table('../total_kmers_per_sample.csv', sep = '\t')

n.cores = 15

# Get the phenotype data
cat('Loading metadata and Ace1 data\n')
meta <- read.table('/home/elucas/Liverpool/phase2.AR1/samples/samples.meta.txt', sep = '\t', header = T, row.names = 1, quote = '')
CI.meta <- subset(meta, population == 'CIcol')
CI.alive <- grepl('aPM', CI.meta$src_code)
names(CI.alive) <- row.names(CI.meta)

# Function to run spearman's rank correlation of a given kmer against the resistance phenotype. 
analysis.function <- function(a)
	cor.test(as.numeric(a), as.numeric(CI.alive), method = 'spearman')$p.value

# Find all the Rdata files. 
all.files <- list.files('..', full.names = T)
all.R.environments <- all.files[grep('.*norecode.*Rdata', all.files)]

# Open each of the kmer counts .Rdata environments in turn and run spearman's rank correlation against the 
# phenotype for every single kmer
p.values <- list()
for (e in all.R.environments[1:length(all.R.environments)]){
	cat('Loading environment ', e, '\n', sep = '')
	load(e)
	# Double check that the order of samples is the same for kmers and phenotypes
	sample.names <- colnames(d)[3:ncol(d)]
	if (any(sample.names != names(CI.alive)))
		stop('Sample names different in kmer counts and phenotypes.\n')
	# Replace NA values with 0
	d <- d[, (sample.names) := lapply(.SD[, 3:ncol(.SD)], function(x) {x[is.na(x)] <- 0; x})]
	# Normalise the kmer counts in d
	d[, 3:ncol(d)] <- data.frame(mapply('/', d[, 3:ncol(d)], total.kmers$Total.kmers)*1e9)
	# Use parallel processing to perform the spearman's rank correlations on all the kmers.
	table.chunks <- floor(seq(0, nrow(d), length.out = n.cores + 1))
	table.chunk.sizes <- table.chunks[2:length(table.chunks)] - table.chunks[1:(length(table.chunks)-1)]
	groupings <- rep(LETTERS[1:n.cores], table.chunk.sizes)
	spear.test.p <- unlist(mclapply(split(d[, 3:ncol(d)], groupings), function(x) unlist(apply(x, 1, analysis.function)), mc.cores = n.cores))
	p.values[[e]] <- cbind(d[, 1:2], spear.test.p)
	rm(d)
	rm(groupings)
	rm(spear.test.p)
	rm(table.chunks)
	rm(table.chunk.sizes)
	gc()
}

all.p.values <- do.call(rbind, p.values)
all.p.values$index <- 1:nrow(all.p.values)

# Identify the significant kmers
non.na.indices <- which(!is.na(all.p.values$spear.test.p))
non.na.p <- all.p.values$spear.test.p[non.na.indices]
all.fdr <- fdrtool(non.na.p, statistic = 'pvalue')$qval
significant.indices <- non.na.indices[all.fdr < 0.001]
sig.table <- all.p.values[significant.indices, ]

# Convert the numeric kmer representation back to a nucleotide string
toseq <- function(numseq)
	gsub('1', 'A', gsub('2', 'C', gsub('3', 'G', gsub('4', 'T', numseq))))
sig.kmer.num.strings <- paste(as.character.integer64(sig.table[[1]]), as.character.integer64(sig.table[[2]]), sep = '')
sig.kmer.strings <- toseq(sig.kmer.num.strings)
# Give each kmer a name
names(sig.kmer.strings) <- paste('Sig_kmer', str_pad(1:length(sig.kmer.strings), nchar(length(sig.kmer.strings)), 'left', '0'), sep = '')

# Output to a table
cat('Writing kmers to "full_sig_kmers.txt"\n')
write.table(sig.table, 'full_sig_kmers.txt', sep = ',', col.names = F, row.names = F)

# Now write it as a fasta file
kmer.fasta <- function(kmers, fastafile){
	ff <- file(fastafile, 'w')
	for (i in 1:length(kmers)){
		write(paste('>', names(kmers)[i], sep = ''), ff)
		write(kmers[i], ff)
	}
	close(ff)
}
cat('Writing kmers to "full_sig_kmers.fa"\n')
kmer.fasta(sig.kmer.strings, 'full_sig_kmers.fa')

# Now get the 100 most significant kmers
most.sig.table <- all.p.values[order(all.p.values$spear.test.p)[1:100], ]
most.sig.kmer.num.strings <- paste(as.character.integer64(most.sig.table[[1]]), as.character.integer64(most.sig.table[[2]]), sep = '')
most.sig.kmer.strings <- toseq(most.sig.kmer.num.strings)
names(most.sig.kmer.strings) <- paste('Most_sig_kmer', str_pad(1:length(most.sig.kmer.strings), nchar(length(most.sig.kmer.strings)), 'left', '0'), sep = '')

cat('Writing kmers to "most_sig_kmers.txt"\n')
write.table(most.sig.table, 'most_sig_kmers.txt', sep = ',', col.names = F, row.names = F)

cat('Writing kmers to "most_sig_kmers.fa"\n')
kmer.fasta(most.sig.kmer.strings, 'most_sig_kmers.fa')

save.image('full_phen_assoc.Rdata')

cat('Process time at end:\n')
print(proc.time())
cat('\n\n')

