# Additional scripts

Some additional scripts. These are used to produce various panels in figs 1, 2 and 10 and tables 1 & 2; but they are not reqired for the primary analyses.

* `2020-02-17_heatmaps.R`: produce tables and heatmaps of the frequencies of *119S* alleles and *Ace1* duplications in each *Ag1000G* population.

* `2020-02-21_phe-gen_CIcol.R`: runs **genotype-phenotype association** tests for the combinations of CNVs and mutations in the `CIcol` population. First, it runs independent Fisher's exact tests for each mutation. Then, it uses a GLM binomial model and step-wise elimination of redundant variants using the BIC criterion to identify the minimal significant model.

Calaix de sastre.
