# This script loads up the full table data, but only keeps the ones defined in "significant_matching_kmers.h5".
# This creates a much smaller table, which is then saved will the full kmer count data preserved. 

library(rhdf5)
library(data.table)
library(bit64)
options(scipen = 999)

start.time <- proc.time()
cat('Process time at start:\n\n')
print(start.time)

all.files <- list.files('.')
all.R.environments <- all.files[grep('.*norecode.*Rdata', all.files)]

# Loop through all the kmer tables and count the total number of k-mers in each sample for each sub-table.
sums.list <- list()
for (e in all.R.environments[1:length(all.R.environments)]){
	cat('Loading environment ', e, '\n', sep = '')
	load(e)
	sums.list[[e]] <- apply(d[, 3:ncol(d)], 2, sum, na.rm = T)
	rm(d)
	gc()
}

# Now get the total number of kmers in each sample
total.kmers <- apply(do.call(rbind, sums.list), 2, sum)

# Save this object to file
write.table(total.kmers, file = 'total_kmers_per_sample.csv', sep = '\t', col.names = 'Total.kmers', quote = F)

cat('Process time at end:\n')
print(proc.time())
cat('\nTime taken:\n')
print((proc.time() - start.time)['elapsed']/60)
cat('\n\n')

