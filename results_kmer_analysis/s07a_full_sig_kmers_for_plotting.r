library(Biostrings)
options(scipen=999)

all.sig.kmers <- read.table('full_sig_kmers.csv', sep = '\t', header = T, stringsAsFactors = F)
all.sig.kmers <- all.sig.kmers[order(all.sig.kmers$Chrom, all.sig.kmers$Pos), ]

cat('\nThere were ', nrow(all.sig.kmers), ' alignments in total, from a total of ', length(unique(all.sig.kmers$Kmer)), 
    ' assembled sequences.\n', sep = '')

# Get the median length of the assembled kmers
assembled.kmers <- readDNAStringSet('/home/eric/Liverpool/ML/CI_fastq/full_analysis/full_phen_assoc_standalone/full_sig_merged_kmers.fa')
cat('The median size of the assembled kmers was', median(width(assembled.kmers)), 'bp.\n')

# Load the table relating each assembled sequence to the kmers that constitute it
kmer.dict <- read.table('full_sig_kmer_dictionary.txt', row.names = 1, stringsAsFactors = F, col.names = c('Sequence', 'kmers'))
# For each assembled sequence, get the number of kmers that assembled to it
kmer.dict$num.kmer <- unlist(lapply(strsplit(kmer.dict$kmers, '\\.'), length))
# Identify the sequences composed of at least 2 kmers
filtered.sequences <- rownames(kmer.dict)[kmer.dict$num.kmer > 1]
filtered.sig.kmers <- subset(all.sig.kmers, Kmer %in% filtered.sequences)

cat('\nAfter filtering sequences composed of only 1 k-mer, there were ', nrow(filtered.sig.kmers), ' alignments in ',
    'total, from a total of ', length(unique(filtered.sig.kmers$Kmer)), ' assembled sequences.\n', sep = '')

# Get the median length of the assembled kmers
filtered.assembled.kmers <- assembled.kmers[names(assembled.kmers) %in% filtered.sequences]
cat('The median size of the assembled kmers after filtering was', median(width(filtered.assembled.kmers)), 'bp.\n')

# Get rid of the supplementary alignments (ie: the alignment of the soft-clipped reads)
cat(sum(filtered.sig.kmers$Supplementary == 'True'), ' of the alignments were supplementary. We will now consider only primary alignments.\n', sep = '')
sig.kmers <- subset(filtered.sig.kmers, Supplementary == 'False')

# Get rid of the NA alignments
cat(sum(is.na(sig.kmers$Chrom)), ' assembled kmers did not align.\n', sep = '')
sig.kmers <- subset(sig.kmers, !is.na(Chrom))

# Get rid of the alignments to the unknown chromosome
chroms <- c('2L', '2R', '3L', '3R', 'X')
cat(sum(sig.kmers$Chrom == 'UNKN'), ' assembled kmers aligned to the unknown chromosome.\n', sep = '')
cat(sum(sig.kmers$Chrom == 'Y_unplaced'), ' assembled kmers aligned to the Y chromosome.\n', sep = '')
cat(sum(sig.kmers$Chrom %in% chroms), ' assembled kmers aligned to one of the five main chromosome arms.\n', sep = '')
sig.kmers <- subset(sig.kmers, Chrom %in% chroms)

# Get the proportion of the kmers aligned to Dup1
num.Dup1 <- sum(sig.kmers$In_Dup1 == 'True')
prop.Dup1 <- (num.Dup1 / nrow(sig.kmers))
cat(num.Dup1, ' (', round(100*prop.Dup1, 1), '%) of these ', nrow(sig.kmers), ' assembled kmers aligned to the Ace1_Dup1 region.\n', sep = '')

# Let's plot these kmers along the genome. Here we just load the Anopheles gambiae genome.
genome <- readDNAStringSet('/home/eric/Liverpool/AR3/genome/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa')
chrom.size <- width(genome)[c(2,1,4,3,6)]
names(chrom.size) <- c('2R', '2L', '3R', '3L', 'X')

# Work out some genomic windows
window.size <- 1000000
genomic.windows <- sapply(chrom.size, function(x, w = window.size) c(seq(0, x, w), x))

# Now count the number of reads occuring in each window
ff <- function(x) hist(subset(sig.kmers, Chrom == x)$Pos, breaks = genomic.windows[[x]], plot = F)
window.hist <- lapply(names(genomic.windows), ff)
names(window.hist) <- names(genomic.windows)
window.kmer.counts <- lapply(window.hist, function(x) x$counts)
window.mids <- lapply(window.hist, function(x) x$mids)

# Save the counts and the window coordinates as a table
windowed.kmer.counts.table <- data.frame(Chrom = unlist(lapply(names(window.kmer.counts), function(x) rep(x, length(window.kmer.counts[[x]])))),
                                         Window.midpoint = unlist(window.mids),
                                         Kmer.num = unlist(window.kmer.counts))

write.table(windowed.kmer.counts.table, 'windowed_kmer_counts_table.csv', sep = '\t', col.names = T, row.names = F, quote = F)
