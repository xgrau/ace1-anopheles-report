# Genetic diversity and evolution of *Ace1* in *Ag1000G*

## What is this

The scripts and data in this repository can be used to reproduce all analyses from the manuscript [**Resistance to pirimiphos-methyl in West African *Anopheles* is spreading via duplication and introgression of the Ace1 locus**](https://www.biorxiv.org/content/10.1101/2019.12.17.879775v1) (Grau-Bové et al., bioRxiv 2020).

Genome variation data for this project has been generated as part of the [***Anopheles gambiae* 1000 Genomes Consortium**](https://www.malariagen.net/projects/ag1000g).

## Analyses

### Genetic diversity & evolution of *Ace1*

All main analyses are organised in `ipython` notebooks, that you can run/examine in the following order:

1. `s01_haplotype_analysis_Ace1_2020-02-14.ipynb` can be used to calculate genotype frequencies, build haplotype networks, perfom positive selection scans along the gene & chromosome, and to obtain haplotype alignments. Output goes to `results_hap_analysis`.

2. `s02_admixture_Ace1_arab_2020-02-14.ipynb`: perform Patterson's D tests of introgression (aka ABBA-BABA test) between various pairs of populations. Output goes to `results_admixture`. There are four scripts, using *arab*, *quad*, *mela* and *meru* as outgroups.

3. `s03_alignments_haplotypes_Ace1_2020-02-14.ipynb`: produce alignments of haplotypes in the *Ace1* duplication and two control regions upstream and downstream of it. Results go to `results_admixture_phylo`. This folder also contains log files from iqtree ML phylogenetic analyses, and a `R` script to create plots for each phylogeny(`00_plot_trees_2020-02-14.R`).

4. `s04_popgen_CIcol_PCA_2020-02-14.ipynb`: calculate genetic differentiation and selection statistics between the PM-resistant and PM-susceptible subpopulations of *A. coluzzii* from Côte d'Ivoire. Results go to `results_differentiation_CIcol`.

These scripts are available as `ipython` notebooks (you can open them here on github, using jupyter notebooks, or VSCode).

### Genotype-phenotype association

Folder `results_phe-gen_CIcol_9WA` contains scripts and data to run genotype-phenotype assocation analyses for the CIcol samples and 9 additional West African populations of *gam* and *col*.

### Other

* `results_gene_phylo`: alignments and phylogenetic analyses of acetylcholinesterase homologs from multiple animal species (list and data sources available as `SM16`), used to establish homology of *Ace1* mutations across species.

* `results_kmer_analysis`: scripts and commands to count k-mers, detect association with phenotypes, and assemble and map the significant k-mers to the genome.

* `results_tables`: additional scripts used to produce heatmaps and other figures (not required for primary analyses).

Other folders:

* `metadata` folder with metadata for the scripts above (sample info, karyotypes, genome annotations,etc.).

* `scripts_hapclust`, `scripts_printtranscripts`, `scripts_other`: some helper functions used by the main scripts.

## Data

Where is the input data?

* All metadata required is in the `metadata` folder
* Population genomic data **has to be downloaded** from the [Ag1000G project archive](https://www.malariagen.net/projects/ag1000g). These are huge files that don't fit in this repository. Download links for Phase1-AR3 and Phase2-AR1:

```bash
ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/
ftp://ngs.sanger.ac.uk/production/ag1000g/phase2/AR1/
```

Notes on data download:

* All genome genome variation ised in the scripts above need to be **specified at the beginning of each python notebook**. Once you've downloaded them, edit the scripts to point to the relevant files.
* Data is available for download in various formats (VCFs, zarr, and HDF5). The scripts above use the zarr arrays and HDF5 files, which are highly compressed and very handy to use compared to VCFs. The python scripts require some special libraries to deal with these formats, mostly implemented in the `scikit-allel`, `zarr` and `h5py` libraries (see dependencies below).
* **phased variants** are available under the `haplotype/main` subfolder:

```bash
ftp://ngs.sanger.ac.uk/production/ag1000g/phase2/AR1/haplotypes/main/
ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/haplotypes/main/
```

* nucleotide **accessibility arrays** in HDF5 format:

```bash
ftp://ngs.sanger.ac.uk/production/ag1000g/phase2/AR1/accessibility/
ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/accessibility/
```

* other **metadata** files:

```bash
ftp://ngs.sanger.ac.uk/production/ag1000g/phase2/AR1/samples/
ftp://ngs.sanger.ac.uk/production/ag1000g/phase1/AR3/samples/
```

## Dependencies

**Python** notebooks work with Python 3.7.4 and the following libraries, which can all be installed using `conda`:

* numpy 1.17.3
* zarr 2.3.2
* pandas 0.25.3
* scikit-allel, allel 1.2.1
* h5py 2.10.0
* scipy 1.3.2
* bcolz 1.2.1
* matplotlib 3.1.2
* seaborn 0.9.0
* itertools 7.2.0

**R scripts** work with R 3.6.1 and require the following libraries:

* seqinr 3.4-5
* ape 5.3
* phytools 0.6-60
* pheatmap 1.0.12

If you use these scripts in your own work, please do not forget to cite the relevant packages as well. It's free and it makes everyone happy :)

For example, in R:

```R
> citation("ape")

To cite ape in a publication use:

  Paradis E. & Schliep K. 2018. ape 5.0: an environment for modern
  phylogenetics and evolutionary analyses in R. Bioinformatics 35:
  526-528.

A BibTeX entry for LaTeX users is

  @Article{,
    title = {ape 5.0: an environment for modern phylogenetics and evolutionary analyses in {R}},
    author = {E. Paradis and K. Schliep},
    journal = {Bioinformatics},
    year = {2018},
    volume = {35},
    pages = {526-528},
  }

As ape is evolving quickly, you may want to cite also its version
number (found with 'library(help = ape)' or
'packageVersion("ape")').

```
