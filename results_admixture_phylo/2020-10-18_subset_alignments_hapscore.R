# libraries
library(ape)
library(stringr)

# input files
ali_fn = "hapalignment_duplication.fasta"
alp_fn = "hapalignment_duplication.fasta.csv"
hsd_fn = "../results_tables/dupcoverage_stats.HaplotypeScore.csv"

# input vars
ace_start = 3489213 # start gene
ace_end   = 3493788 # end gene
ace_119S  = 3492074 # resistant variant
ace_dups  = 3436800 # start duplication
ace_dupe  = 3639600 # end duplication

# load
hsd = read.table(hsd_fn, sep="\t", header = T)
alp = read.table(alp_fn, sep="\t", header = T)
ali = read.FASTA(ali_fn)
