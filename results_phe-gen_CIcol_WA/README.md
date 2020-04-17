# Genotype-phenotype analysis

* `2020-02-21_phe-gen_CIcol.R`: runs **genotype-phenotype association** tests for the combinations of CNVs and mutations in the `CIcol` population. First, it runs independent Fisher's exact tests for each mutation. Then, it uses a GLM binomial model and step-wise elimination of redundant variants using the BIC criterion to identify the minimal significant model. These plots are used in Fig3.

* `2020-02-17_phe-gen_9WApops_280S.R`: runs **genotype-phenotype association** tests for G280S alleles and pirimiphos-methyl resistance, for each of the 9 populations of West African *col* and *gam* (Fisher's exact test).

* `XXXX.R`: Same for CNVs?
