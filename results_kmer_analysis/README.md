# $k$-mer analysis

Recipe to identify $k$-mers from `fastq` that are significantly associated with pirimiphos-methyl resistance, using 71 WGS data for *A. coluzzii* from Côte d'Ivoire.

## Methods

0. **Obtain $k$-mer counts** for each sample `fastq` file using `jellyfish` v. 2.2.10 [Marcais and Kingsford 2011]:

```bash
jellyfish count -C -m 31 --out-counter-len 2 -s 300M --bf-size 10G
```

1. Obtain **lexographical $k$-mer groups**. To reduce the computer memory footprint of the k-mer count tables, the k-mer strings were recoded as integers, split into separate files according to their leading nucleotides ($k$-mers beginning with `AAA` were saved into one file, those beginning with `AAC` into another file, and so on) and, within each file, sorted lexographically. We call the set of $k$-mers within each file a "lexographical group".

```bash
# run in a loop using as the first argument each of the
# lexographical groups, e.g. AAA, then AAC, etc.
python s01_recode_split_sort.py <file_list> <path_to_jellyfish> <num_threads>
```

2. The resulting tables were loaded into `R`, and $k$-mers present in fewer than 3 samples, or absent in fewer than 3 samples, were discarded in order to reduce the size of the dataset while keeping the K-mers most likely to show variation in the poulation (called "variant k-mers"). This was done separately for each lexographical group, using the script:

```bash
# run in a loop using as the first argument each of the
# lexographical groups, e.g. AAA, then AAC, etc.
Rscript s02_join_tables.r <lexographical_group>

# reduce sample set by grouping together kmers with the same presence/absence profile
# output: smalltable.Rdata
Rscript s03_combine_all_tables.r
```

4. **Significance of phenotypic association** was determined independently for each marker using Fisher's exact test with false discovery rate correction, implemented in the script:

```bash
# input: smalltable.Rdata
# ALSO: ensure that it points to the correct metadata file.
# output: smalltable.h5, containing all fisher's tests in a compressed Hdf5 format
Rscript s04_fisher_testing.r
```

5. To **test whether the frequencies of any k-mers were associated with PM resistance**, we normalise k-mer frequencies. Then, we calculated Spearman's rank correlation against resistance phenotype independently on all 767,560,108 k-mers, with false discovery rate control at 0.1%. This was done using the scripts:

```bash
# output: total_kmers_per_sample.csv
Rscript s05a_get_total_kmers_per_sample.r

# input: total_kmers_per_sample.csv, and each of the Rdata objects with kmer frequencies
# output: most_sig_kmers.txt, most_sig_kmers.fa, most_sig_kmers.fa, full_phen_assoc.Rdata
Rscript s05b_full_phen_assoc_standalone.r
```

6. **$k$-mer assembly**. We took the k-mers that showed a significant association with PM resistance and assembled them by joining any k-mers that overlapped perfectly over at least 10 bp. This was implemented using the script:

```bash
# input: full_sig_kmers.fa
# output: full_sig_merged_kmers.fa, full_sig_kmer_dictionary.txt
python s06_assemble_full_sig_kmers.py
```

7. **$k$-mer alignment**. The resulting assembled k-mers were aligned against the AgamP3 reference genome using `bwa mem`. Results were then prepared into tables for plotting with the scripts:

```bash
bwa mem -T 0 Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa full_sig_merged_kmers.fa
python s07a_investigate_full_sig_kmers.py
Rscript s07b_full_sig_kmers_for_plotting.r
```

8. **$k$-mer frequency correlation with *Ace1***: check whether the frequencies of each significant $k$-mer in each sample correlate with *Ace* copy numbers.

```bash
Rscript s08a_extract_full_sig_kmers.r
s08b_Ace1_correlation.r
```

9. Produce table with counts of $k$-mer mapping coordinates along the genome:

```bash
python s09_2020-02-20_kmerplot_coordinates_CIcolb.py
```

## Unused

Assembly with `abyss`:

```bash
abyss-pe k=25 name=full_sig_kmers_assembly se="full_sig_kmers.fa"
```

BLAT alignments:

```bash
blat Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa full_sig_merged_kmers.fa full_sig_merged_kmers.blat -out=blast8
```