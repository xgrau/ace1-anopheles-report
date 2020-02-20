# This script loads up the full table data, but only keeps the ones previously identified as significant.
# This creates a much smaller table, which is then saved will the full kmer count data preserved. 

library(rhdf5)
library(data.table)
library(bit64)
library(stringr)
options(scipen = 999)

start.time <- proc.time()
cat('Process time at start:\n\n')
print(start.time)

# Load the total number of kmers per sample
cat('Loading sample-wise kmer totals.\n')
total.kmers <- read.table('../total_kmers_per_sample.csv', sep = '\t')

sig.kmers <- read.table('full_sig_kmers.txt', sep = ',')
colnames(sig.kmers) <- c('kmer.string.1', 'kmer.string.2', 'p.value', 'index')
# Give the kmers the same names as those given in the fasta file
rownames(sig.kmers) <- paste('Sig_kmer', str_pad(1:nrow(sig.kmers), nchar(nrow(sig.kmers)), 'left', '0'), sep = '')

all.files <- list.files('..', full.name = T)
all.R.environments <- all.files[grep('.*norecode.*Rdata', all.files)]

# Load each kmer table and pull out the kmers previously identified as significant
sig.tables.list <- list()
kmers.so.far <- 0
for (e in all.R.environments[1:length(all.R.environments)]){
	cat('Loading environment ', e, '\n', sep = '')
	load(e)
	sample.names <- colnames(d)[3:ncol(d)]
	d <- d[, (sample.names) := lapply(.SD[, 3:ncol(.SD)], function(x) {x[is.na(x)] <- 0; x})]
	these.kmers <- subset(sig.kmers, index > kmers.so.far & index <= (kmers.so.far + nrow(d)))
	these.indices <- these.kmers$index - kmers.so.far
	kmers.so.far <- kmers.so.far + nrow(d)
	dd <- d[these.indices, ]
	# Normalise the kmer counts in dd
	dd[, 3:ncol(dd)] <- data.frame(mapply('/', dd[, 3:ncol(dd)], total.kmers$Total.kmers)*1e9)
	dd$Kmer.name <- rownames(these.kmers)
	sig.tables.list[[e]] <- dd[, c('Kmer.name', 'V1', 'V2', sample.names), with = F]
	rm(d)
	rm(dd)
	gc()
}

# Join the outputs into a single table
s <- do.call(rbind, sig.tables.list)
# Order the table by kmer string. The order function can only take one object at a time when dealing 
# with int64, so we need to order in two steps
ss <- s[order(s$V2), ]
sig.table <- ss[order(ss$V1), ]

write.table(sig.table, 'full_sig_kmer_table.csv', sep = '\t', col.names = T, row.names = F, quote = F)

cat('Process time at end:\n')
print(proc.time())
cat('\nTime taken:\n')
print((proc.time() - start.time)['elapsed']/60)
cat('\n\n')

