library(rhdf5)
library(data.table)
library(bit64)
options(scipen = 999)

start.time <- proc.time()
cat('Process time at start:\n\n')
print(start.time)

all.files <- list.files('.')
all.R.environments <- all.files[grep('norecode.*Rdata', all.files)]

# Load each lexogrphical group and convert kmer counts to presence/absence
tables.list <- list()
for (e in all.R.environments[1:length(all.R.environments)]){
	cat('Loading environment ', e, '\n', sep = '')
	load(e)
	sample.names <- colnames(d)[3:ncol(d)]
	tables.list[[e]] <- d[, (sample.names) := lapply(.SD[, 3:ncol(.SD)], function(x) as.integer(!is.na(x)))]
	rm(d)
	gc()
}

# Combine all lexographical groups
full.table <- do.call(rbind, tables.list)
rm(tables.list)
gc()

rdata.filename = 'fulltable.Rdata'
cat('\nSaving table as ', rdata.filename, '\n', sep = '')
save(full.table, file = rdata.filename)

# This first part of saving to h5 needs to be before we remove the first two columns of full.table
h5.filename = 'smalltable.h5'
h5createFile(h5.filename)

cat('\nSaving first half of k-mer strings to ', h5.filename, ' object kmer_string_1\n', sep = '')
h5createDataset(file = h5.filename, dataset = 'kmer_string_1', dims = nrow(full.table), chunk = 1e6, H5type = 'H5T_NATIVE_INT64')
h5write(as.double(full.table[[1]]), h5.filename, 'kmer_string_1')
gc()
#
cat('\nSaving second half of k-mer strings to ', h5.filename, 'object kmer_string_2\n', sep = '')
h5createDataset(file = h5.filename, dataset = 'kmer_string_2', dims = nrow(full.table), chunk = 1e6, H5type = 'H5T_NATIVE_INT64')
h5write(as.double(full.table[[2]]), h5.filename, 'kmer_string_2')
gc()

# Free up some computer memory by deleting the columns that we have just written from the data.table
full.table[, c('V1', 'V2'):=NULL]
gc()

tobin <- function(x){
	sum(x*(extrabinvec[1:length(x)]))
}

extrabinvec <- 2^(0:70)

# To make it easier and faster to group kmers with the same presence/absence profile, we convert the data for
# each kmer into two numbers, which together act as bit-wise flags to code the binary presence/absence data.
# If we try to do the following in one go on full.table, it fails because it the object is too large for it
# to allocate, even though the actual object size is much smaller than the amount of memory it says it needs. 
# We therefore need to do it in chunks instead. 
Out1 <- list()
Out2 <- list()
start.time2 <- proc.time()
table.chunks <- 1e8
for (i in 1:ceiling(nrow(full.table)/table.chunks)){
	print(i)
	start.index <- (i-1)*table.chunks + 1
	end.index <- ifelse((i*table.chunks) > nrow(full.table), nrow(full.table), i*table.chunks)
	Out1[[i]] <- apply(full.table[start.index:end.index, 36:1], 1, tobin)
	Out2[[i]] <- apply(full.table[start.index:end.index, 71:37], 1, tobin)
	print(object.size(Out1))
	cat('Time taken so far: ', (proc.time() - start.time2)['elapsed']/60, '\n', sep = '')
	gc()
}

out <- cbind(do.call(c, Out1), do.call(c, Out2))
rm(Out1)
rm(Out2)
gc()

# Keep only one entry for each presence/absence profile
cat('Removing duplicate rows\n')
start.time3 <- proc.time()
ss <- data.table(cbind(out, 1:nrow(full.table)))
colnames(ss) <- c('out1', 'out2', 'index')
ind <- ss[, c(.SD[1,], .N), by = c('out1', 'out2')]
colnames(ind)[4] <- 'counts'
small.table <- cbind(ind$out1, ind$out2, ind$counts, full.table[ind$index, ])
colnames(small.table)[1:3] <- c('kmer.key1', 'kmer.key2', 'counts')
cat('Time taken for removing duplicate rows: ', (proc.time() - start.time3)['elapsed']/60, '\n', sep = '')


rm(full.table)
rm(ss)
gc()

save.image('smalltable.Rdata')

# Now write the rest as the final object in the h5 file
# Because the hdf5 in R inverts rows and columns (for computational reasons), we would need to give nrow for the
# number of columns and ncol for the number of rows, and then give "native = T" when writing the table. This 
# will give data that can be properly read into python. It can also be properly read into R by including native = T 
# in the read function. However, it seems that there is a limit to the number of rows that can be saved if we do
# it this way, so it looks like we may need to stick to the "wrong" orientation and deal with that downstream. 

cat('\nSaving kmer groupings to ', h5.filename, 'object kmer_groups\n', sep = '')
h5createDataset(file = h5.filename, dataset = 'kmer_groups', dims = c(nrow(out), ncol(out)), chunk = c(1e6, ncol(out)))
h5write(as.matrix(out), h5.filename, 'kmer_groups')
gc()

cat('\nSaving sample names to ', h5.filename, 'object sample_names\n', sep = '')
sample.names <- colnames(small.table)[4:ncol(small.table)]
h5createDataset(file = h5.filename, dataset = 'sample_names', dims = length(sample.names), storage.mode = 'character', size = max(nchar(sample.names)) + 1)
h5write(sample.names, h5.filename, 'sample_names')
gc()

cat('\nSaving kmer summary data to ', h5.filename, 'object kmer_presense_summary\n', sep = '')
h5createDataset(file = h5.filename, dataset = 'kmer_presence_summary', dims = c(nrow(small.table), 3), chunk = c(1e6, 3))
h5write(as.matrix(small.table[,1:3]), h5.filename, 'kmer_presence_summary')#, native = T)
gc()

cat('\nSaving kmer presence/absence data to ', h5.filename, 'object kmer_presence\n', sep = '')
h5createDataset(file = h5.filename, dataset = 'kmer_presence', dims = c(nrow(small.table), ncol(small.table) - 3), chunk = c(1e6, ncol(small.table) - 3), storage.mode = 'integer')
h5write(as.matrix(small.table[,4:ncol(small.table)]), h5.filename, 'kmer_presence')#, native = T)
gc()

cat('Process time at end:\n')
print(proc.time())
cat('\nTime taken:\n')
print((proc.time() - start.time)['elapsed']/60)
cat('\n\n')

