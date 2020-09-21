# import libraries
import allel
import zarr
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import itertools
import logging
import scipy.stats

# input files
metasam_fn = "metadata/samples.meta_phenotypes.txt"
callset_fn = "/media/xavi/Saigon/Variation/ag1000g.phase2.ar1"
results_fo = "results_tables"

## LOAD METADATA
# population metadata
samples = pd.read_csv(metasam_fn, sep='\t')

# sample: haplotypes
sampleh = pd.DataFrame(data = {
	"haplotype" : np.array(list(itertools.chain(*[[s + ' a', s + ' b'] for s in samples["ox_code"] ]))),
	"population": np.array(list(itertools.chain(*[[s , s] for s in samples["population"] ]))),
	"phenotype": np.array(list(itertools.chain(*[[s , s] for s in samples["phenotype"] ])))
})

## POPULATION DICTIONARIES
# population dictionary
poplist = ["AOcol","BFcol","BFgam","CIcol","CMgam","FRgam","GAgam","GHcol","GHgam","GM","GNcol","GNgam","GQgam","GW","KE","UGgam"]
popdict = dict()
popdict["all"] = []
for pop in poplist:
	popdict[pop]   = samples[  samples["population"] == pop ].index.tolist()
	popdict["all"] = popdict["all"] + popdict[pop]


callset_fn = "/home/xavi/Documents/VariationAg1k/data/phase2.AR1/variation/main/hdf5/ag1000g.phase2.ar1.2R.h5"

import h5py
callset     = h5py.File(callset_fn,mode="r")
callset["2R"]["variants"]
# callset = zarr.open("%s/%s" % (callset_fn, chrom))
genvars = allel.VariantChunkedTable( callset["2R"]["variants"], names=["POS","REF","ALT"],index="POS" )
genotyp = allel.GenotypeChunkedArray(callset["2R"]["calldata"]["genotype"])
samples = callset["2R"]["samples"][:].astype(str)

# snps I care about
snps_loc = [3489316,            3489319,           3489344,           3489394,           3489405,          3489435,          3492074,           3493750]
snps_nom = ["104T>A Phe35Tyr", "107C>A Ala36Glu", "132G>T Glu44Asp", "182G>A Gly61Asp", "193G>T Ala65Ser","223G>A Ala75Thr","838G>A Gly280Ser","2176G>T Val726Leu"]

# retrieve them
snps_boo = np.isin(element=genvars["POS"], test_elements=snps_loc)

genvars = genvars.compress(snps_boo)
genotyp = genotyp.compress(snps_boo)


# export genotypes
gennalt_df = pd.DataFrame(columns=snps_loc, index=samples)
genotyp_df = pd.DataFrame(columns=snps_loc, index=samples)
for n,snp in enumerate(snps_loc):
	gennalt_df[snp] = genotyp.to_n_alt()[n]
	genotyp_df[snp] = genotyp[n]
	
genotyp_df.to_csv("results_tables/genotypes_per_sample_Ace1nonsyn_gts.csv",sep="\t")
gennalt_df.to_csv("results_tables/genotypes_per_sample_Ace1nonsyn_nalt.csv",sep="\t")


# export variants and alt/ref
genvars_df = pd.DataFrame(columns=["chrom","POS","REF"])
genvars_df["POS"] = genvars["POS"]
genvars_df["REF"] = genvars["REF"][:].astype(str)
genvars_df["ALT1"] = genvars["ALT"][:,0].astype(str)
genvars_df["ALT2"] = genvars["ALT"][:,1].astype(str)
genvars_df["ALT3"] = genvars["ALT"][:,2].astype(str)
genvars_df["chrom"] = "2R"
genvars_df.to_csv("results_tables/genotypes_per_sample_Ace1nonsyn_vars.csv",sep="\t")
