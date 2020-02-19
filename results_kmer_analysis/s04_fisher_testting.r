# This script runs Fisher's exact test on the presence / absence data. The R workspace smalltable.Rdata is
# created by the script combine_all_tables.r
# This is just a subset of presence_absence_phen_assoc.r, but keeping only the parts of the code relevant
# to phenotype assocaition testing, which is more suitable for inclusion in the Ace-1 paper. 

library(rhdf5)
library(data.table)
library(parallel)
library(fdrtool)
options(scipen = 999)

n.cores = 15

cat('\nLoading R workspace\n')
load('smalltable.Rdata')

cat('Loading metadata\n')
# Load the Ag1000G phase2 metadata
meta <- read.table('../metadata/samples.meta_phenotypes.txt', sep = '\t', header = T, row.names = 1, quote = '')
CI.meta <- subset(meta, population == 'CIcol')
CI.alive <- grepl('aPM', CI.meta$src_code)
names(CI.alive) <- row.names(CI.meta)

# Double check that the order of samples is the same for both cases
any(colnames(small.table)[4:ncol(small.table)] != names(CI.alive))

# Run Fisher's exact tests on each of the k-mer profiles
cat('\nTesting kmers for phenotypic association\n')
sm.table.chunks <- floor(seq(0, nrow(small.table), length.out = n.cores + 1))
sm.table.chunk.sizes <- sm.table.chunks[2:length(sm.table.chunks)] - sm.table.chunks[1:(length(sm.table.chunks)-1)]
sm.groupings <- rep(LETTERS[1:n.cores], sm.table.chunk.sizes)
f <- unlist(mclapply(split(small.table[, 4:ncol(small.table)], sm.groupings), function(x) unlist(apply(x, 1, function(a) fisher.test(table(a, CI.alive))$p.value)), mc.cores = n.cores))

# It seems that there are rounding errors in f, so the highest values are infinitessimally greater than 1, 
# which causes problems downstream. So fix that here. 

ff <- round(f, 15)

# Check whether any kmers are significant after correction for multiple testing. 
cat('\tNumber of significant kmers after Bonferroni correction:', sum(f < (0.05/length(ff))), '\n')
cat('\tNumber of significant kmers after FDR control at 5%:', sum(fdrtool(as.numeric(ff), statistic = 'pvalue', plot = F)$qval < 0.05), '\n')

# Save the output to the h5 file.
h5.filename = 'smalltable.h5'
cat('\nSaving P-values to ', h5.filename, 'object Fisher_test_pvalues\n', sep = '')
h5createDataset(file = h5.filename, dataset = 'Fisher_test_pvalues', dims = length(ff), chunk = 1e6)
h5write(ff, h5.filename, 'Fisher_test_pvalues')
gc()




