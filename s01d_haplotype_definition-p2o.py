# import libraries
import allel
import zarr
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import logging
import scipy.stats
import h5py

import os
os.chdir("/home/xavi/Documents/ace1-anopheles-report/")

p2_poplist = ["BFcol","CIcol","GHcol","BFgam","GHgam","GNgam"]
# p2_poplist = ["AOcol","BFcol","CIcol","GHcol","GNcol","BFgam","CMgam","FRgam","GAgam","GHgam","GNgam","GQgam","UGgam","GM","GW","KE"]
# p2_poplist = ["CIcol"]


# input files
p2_metasam_fn = "metadata/samples.meta_phenotypes_acegenotype.simple.txt"
p2_callset_fn = "/home/xavi/Documents/VariationAg1k/data/phase2.AR1/haplotypes/zarr2/ag1000g.phase2.ar1.samples/"

# output file
results_fo = "results_hap_analysis"

# snps I care about
# Ace1 G280S location:
chrom = "2R"
ace_start = 3489213 # start gene
ace_end   = 3493788 # end gene
ace_119S  = 3492074 # resistant variant
ace_dups  = 3436800 # start duplication
ace_dupe  = 3639600 # end duplication


## LOAD METADATA
# population metadata
p2_samples = pd.read_csv(p2_metasam_fn, sep='\t')

# samples with mutation presence
sample_119S_list = p2_samples[p2_samples["119Sgen"] == "119S"]["ox_code"].values

# sample: haplotypes
p2_sampleh = pd.DataFrame(data = {
	"ox_code" : np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in p2_samples["ox_code"] ]))),
	"ox_code_sam" : np.array(list(itertools.chain(*[[s , s ] for s in p2_samples["ox_code"] ]))),
	"population": np.array(list(itertools.chain(*[[s , s] for s in p2_samples["population"] ]))),
	"phenotype": np.array(list(itertools.chain(*[[s , s] for s in p2_samples["phenotype"] ]))),
	"genotype": np.array(list(itertools.chain(*[[s , s] for s in p2_samples["119Sgen"] ])))
})

## POPULATION DICTIONARIES
# population dictionary for phase 2
p2_popdich = dict()
p2_popdict = dict()
p2_popdich["all"] = []
p2_popdict["all"] = []
for pop in p2_poplist:
	# populate haplotype dict
	p2_popdich[pop]   = p2_sampleh[  p2_sampleh["population"] == pop ].index.tolist()
	p2_popdich["all"] = p2_popdich["all"] + p2_popdich[pop]

	# populate samples dict
	p2_popdict[pop]   = p2_samples[  p2_samples["population"] == pop ].index.tolist()
	p2_popdict["all"] = p2_popdict["all"] + p2_popdict[pop]

# load phase2 data
p2_callset = zarr.open(p2_callset_fn)
p2_genvars = allel.VariantChunkedTable( p2_callset[chrom]["variants"], names=["POS","REF","ALT"],index="POS" )
p2_genotyp = allel.GenotypeChunkedArray(p2_callset[chrom]["calldata"]["genotype"])
p2_samlist = p2_callset[chrom]["samples"][:].astype(str)
p2_samlilh = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in p2_samlist ])))

# compress to haplotypes around duplication
flanking_bp = 1e4
p2_var_boo = ((p2_genvars["POS"][:] > ace_dups - flanking_bp) & (p2_genvars["POS"][:] < ace_dups + flanking_bp)) | ((p2_genvars["POS"][:] > ace_dupe - flanking_bp) & (p2_genvars["POS"][:] < ace_dupe + flanking_bp)) | (p2_genvars["POS"] == ace_119S)

# # ...or variants in the duplication?
# flanking_bp = 0
# p2_var_boo = ((p2_genvars["POS"][:] > ace_dups - flanking_bp) & (p2_genvars["POS"][:] < ace_dupe + flanking_bp)) | (p2_genvars["POS"] == ace_119S)

# ...or variants in the downstream region?
# flanking_bp = 1e4
# p2_var_boo = ((p2_genvars["POS"][:] > ace_dupe) & (p2_genvars["POS"][:] < ace_dupe + flanking_bp )) | (p2_genvars["POS"] == ace_119S)

# # ...or variants in the upstream region?
# flanking_bp = 1e5
# p2_var_boo = ((p2_genvars["POS"][:] > ace_dups - flanking_bp) & (p2_genvars["POS"][:] < ace_dups )) | (p2_genvars["POS"] == ace_119S)

# compress
p2_genvars_sub = p2_genvars[:].compress(p2_var_boo)

# expand to haplotyes
p2_genotyp_sub = p2_genotyp.compress(p2_var_boo, axis=0)
p2_haploty_sub = p2_genotyp_sub.to_haplotypes()
# restrict to populations of interest
p2_haploty_sub = np.take(a=p2_haploty_sub, indices = p2_popdich["all"], axis=1)
p2_genotyp_sub = np.take(a=p2_genotyp_sub, indices = p2_popdict["all"], axis=1)

# # get segregating haplotypes from phase 2
# p2_hapalco_sub = p2_haploty_sub.count_alleles()
# p2_is_segh = p2_hapalco_sub.is_segregating()
# p2_haploty_sub_seg = p2_haploty_sub.compress(p2_is_segh, axis=0)
# # same for genotypes
# p2_genalco_sub = p2_genotyp_sub.count_alleles()
# p2_is_segg = p2_genalco_sub.is_segregating()
# p2_genotyp_sub_seg = p2_genotyp_sub.compress(p2_is_segg, axis=0)

# find pure wt genotypes
p2_samples_wt_ix =   np.where( p2_samples["119Sgen"][ p2_popdict["all"] ].values == "wt" )[0]
p2_genotyp_purewt = np.take(a=p2_genotyp_sub, indices = p2_samples_wt_ix, axis=1)
p2_genotyp_purewt_nalt = p2_genotyp_purewt.to_n_alt()[:]
p2_genotyp_purewt_vect = (np.sum(p2_genotyp_purewt_nalt, axis=1) > 0) * 1
p2_genotyp_purewt_vect = np.sum(p2_genotyp_purewt_nalt>0, axis=1)
p2_genotyp_purewt_vect = (np.sum(p2_genotyp_purewt_nalt>0, axis=1) > 0) * 1


#### TRY DISTANCE

# out_purewt_dxy = allel.pairwise_dxy(pos=p2_genvars_sub["POS"], gac=p2_genotyp_purewt.to_allele_counts())
# out_purewt_dxy.shape
# p2_samples_nowt_ix = np.where( p2_samples["119Sgen"][p2_popdict["all"]].values != "wt" )[0]
# p2_sampleh_nowt_ix = np.where( p2_sampleh["genotype"][p2_popdich["all"]].values != "wt" )[0]
# p2_sampleh_wt_ix =   np.where( p2_sampleh["genotype"][p2_popdich["all"]].values == "wt" )[0]

import scipy.spatial
p2_genotyp_purewt_nalt_bin = p2_genotyp_purewt_nalt
p2_genotyp_purewt_nalt_bin[p2_genotyp_purewt_nalt_bin>1] = 1

dis_to_purewt_r = np.zeros(shape = p2_haploty_sub.shape[1])
for i in range(p2_haploty_sub.shape[1]):

	prob=p2_haploty_sub[:,i].tolist()
	arr_prob_puwt = np.transpose(np.vstack((np.transpose(p2_genotyp_purewt_nalt_bin), p2_haploty_sub[:,i] )))
	out_pairdist = allel.pairwise_distance(arr_prob_puwt, metric="cityblock")
	out_pairdisq =  scipy.spatial.distance.squareform(out_pairdist)
	dis_to_purewt_r[i] = np.min(out_pairdisq[-1][:-1])


quantile = .05
top_correlated_quantile = np.quantile(dis_to_purewt_r, q=1-quantile)
bot_correlated_quantile = np.quantile(dis_to_purewt_r, q=quantile)

ixs_top_correlated_haplotyes = np.where(dis_to_purewt_r >= top_correlated_quantile)[0]
ixs_bot_correlated_haplotyes = np.where(dis_to_purewt_r <= bot_correlated_quantile)[0]

# save pdf
plt.figure(figsize=(5,5))
plt.plot(np.sort(dis_to_purewt_r), color="blue")
plt.axhline(y=bot_correlated_quantile, color="red", linestyle='dashed')
plt.axhline(y=top_correlated_quantile, color="red", linestyle='dashed')
plt.title("Mean distance to pure wt sample")
plt.savefig("%s/dist_classification.pdf" % (results_fo), bbox_inches='tight')
plt.close()

p2_sampleh_interest = p2_sampleh.iloc[p2_popdich["all"]]
p2_sampleh_interest = p2_sampleh_interest.reset_index()
p2_resistant_haps = p2_sampleh_interest.iloc[ixs_top_correlated_haplotyes]
p2_resistant_haps.to_csv("%s/dist_classification.resistant_top10p.csv" % (results_fo), sep="\t", index=False)



#### TRY CORRELATION

from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr

cor_to_purewt_r = np.zeros(shape = p2_haploty_sub.shape[1])
for i in range(p2_haploty_sub.shape[1]):

	# cor_to_i = pearsonr(p2_haploty_sub[:,i] , p2_genotyp_purewt_vect)
	# cor_to_purewt_r[i] = cor_to_i[0]
	cor_to_purewt_r[i] = np.max(np.corrcoef(p2_genotyp_purewt_nalt, p2_haploty_sub[:,i], rowvar=False)[-1][:-1])

cor_to_purewt_r = np.nan_to_num(cor_to_purewt_r)


quantile = .1
top_correlated_quantile = np.quantile(cor_to_purewt_r, q=1-quantile)
bot_correlated_quantile = np.quantile(cor_to_purewt_r, q=quantile)

ixs_top_correlated_haplotyes = np.where(cor_to_purewt_r >= top_correlated_quantile)[0]
ixs_bot_correlated_haplotyes = np.where(cor_to_purewt_r <= bot_correlated_quantile)[0]

# save pdf
plt.figure(figsize=(5,5))
plt.plot(np.sort(cor_to_purewt_r), color="blue")
plt.axhline(y=bot_correlated_quantile, color="red", linestyle='dashed')
plt.axhline(y=top_correlated_quantile, color="red", linestyle='dashed')
plt.title("Best correlation to pure wt sample")
plt.savefig("%s/corr_classification.pdf" % (results_fo), bbox_inches='tight')
plt.close()

p2_sampleh_interest = p2_sampleh.iloc[p2_popdich["all"]]
p2_sampleh_interest = p2_sampleh_interest.reset_index()
p2_resistant_haps = p2_sampleh_interest.iloc[ixs_bot_correlated_haplotyes]
p2_resistant_haps.to_csv("%s/corr_classification.resistant_top10p.csv" % (results_fo), sep="\t", index=False)


