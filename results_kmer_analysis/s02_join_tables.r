library(data.table)
options(scipen = 999)

leading.bases <- commandArgs(trailingOnly = T)[1]

cat('Running script on', leading.bases, '\n\n')

cat('Process time at start:\n\n')
print(proc.time())

folders <- list.dirs('..')[grep('-C', list.dirs('..'))]
names(folders) <- sub('.*/', '', folders)

get.data.path <- function(f, i, bp)
	paste(f[i], '/', names(f)[i], '-table', bp, '-sorted.txt', sep = '')

cat('Loading sample', names(folders)[1], '\n')
all.data <- fread(get.data.path(folders, 1, leading.bases), header = F)
colnames(all.data)[3] <- names(folders)[1]

for (i in 2:length(folders)){
	cat('Adding sample', names(folders)[i], '\n')
	all.data <- merge(all.data, fread(get.data.path(folders, i, leading.bases), header = F), by = c('V1', 'V2'), all = T)
	colnames(all.data)[ncol(all.data)] <- names(folders)[i]
	gc()
}

# When I ran this on the Cote d'Ivoire samples, the final all.data object was 38.5 Gb. But in the process of 
# generating that object, the process sometimes reached the 125Gb capacity (not sure why, even in the worst
# case scenario, I feel like the process would at most require twice the size of the object, so ~ 80Gb). I 
# have added a gc() command at the end of the loop to deal with this.
# But Neptune has a 492Gb capacity, so that should be fine, at least for this number of samples. 

# count the number of na entries in every row
cat('Counting NAs\n')
na.counts <- apply(all.data, 1, function(x) sum(is.na(x)))
gc()
num.samples <- ncol(all.data) - 2
num.singletons <- sum(na.counts == (num.samples - 1))
num.ubiquitous <- sum(na.counts == 0)
num.allbutone <- sum(na.counts == 1)

# Keep the k-mers that are present in more than 2 samples, and fewer than n-2
to.keep <- (na.counts > 2) & (na.counts < (num.samples - 2))

cat('\nThere were a total of ', nrow(all.data), ' kmers, of which ', num.ubiquitous, ' were ubiquitous, ',
    num.singletons, ' were singletons and ', num.allbutone, ' were singleton absent. Keeping only kmers ',
    'that were present in MORE than 2 samples and absent in MORE than 2 samples, we ended up keeping ',
    sum(to.keep), ' kmers.\n\n', sep = '')

# Reduce the size of the table 
d <- all.data[to.keep, ]
rm(all.data)
gc()

cat('Saving image\n\n')
save.image(paste('tempnorecode', leading.bases, '.Rdata', sep = ''))

cat('Process time at end:\n\n')
print(proc.time())

