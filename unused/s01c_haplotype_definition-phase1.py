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

# input files
p2_metasam_fn = "metadata/samples.meta_phenotypes_acegenotype.simple.txt"
p1_metasam_fn = "/home/xavi/Documents/VariationAg1k/data/phase1.AR3.1/haplotypes/haplotypes.meta.txt"
p2_callset_fn = "/home/xavi/Documents/VariationAg1k/data/phase2.AR1/haplotypes/zarr2/ag1000g.phase2.ar1.samples/"
p1_callset_fn = "/home/xavi/Documents/VariationAg1k/data/phase1.AR3.1/haplotypes/main/hdf5/ag1000g.phase1.ar3.1.haplotypes.2R.h5"

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
p1_sampleh = pd.read_csv(p1_metasam_fn, sep='\t')

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
p2_poplist = ["BFcol","CIcol","GHcol","BFgam","GHgam","GNgam"]
# p2_poplist = ["BFgam","GHgam","GNgam"]
p2_popdich = dict()
p2_popdich["all"] = []
for pop in p2_poplist:
	p2_popdich[pop]   = p2_sampleh[  p2_sampleh["population"] == pop ].index.tolist()
	p2_popdich["all"] = p2_popdich["all"] + p2_popdich[pop]

# population dictionary for phase 2
p1_poplist = ["BFS","BFM","GNS"]
# p1_poplist = ["BFS","GNS"]
p1_popdich = dict()
p1_popdich["all"] = []
for pop in p1_poplist:
	p1_popdich[pop]   = p1_sampleh[  p1_sampleh["population"] == pop ].index.tolist()
	p1_popdich["all"] = p1_popdich["all"] + p1_popdich[pop]

# load phase1 data
p1_callset = h5py.File(p1_callset_fn,mode="r")
p1_genvars = allel.VariantChunkedTable( p1_callset[chrom]["variants"], names=["POS","REF","ALT"],index="POS" )
p1_genotyp = allel.GenotypeChunkedArray(p1_callset[chrom]["calldata"]["genotype"])
p1_samlist = p1_callset[chrom]["samples"][:].astype(str)
# expand to haplotypes
# p1_samlilh = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in p1_samlist ])))

# load phase2 data
p2_callset = zarr.open(p2_callset_fn)
p2_genvars = allel.VariantChunkedTable( p2_callset[chrom]["variants"], names=["POS","REF","ALT"],index="POS" )
p2_genotyp = allel.GenotypeChunkedArray(p2_callset[chrom]["calldata"]["genotype"])
p2_samlist = p2_callset[chrom]["samples"][:].astype(str)
p2_samlilh = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in p2_samlist ])))

# compress to haplotypes around duplication
flanking_bp = 1e3
p1_var_boo = ((p1_genvars["POS"][:] > ace_dups - flanking_bp) & (p1_genvars["POS"][:] < ace_dups + flanking_bp)) | ((p1_genvars["POS"][:] > ace_dupe - flanking_bp) & (p1_genvars["POS"][:] < ace_dupe + flanking_bp)) | (p1_genvars["POS"] == ace_119S)
p2_var_boo = ((p2_genvars["POS"][:] > ace_dups - flanking_bp) & (p2_genvars["POS"][:] < ace_dups + flanking_bp)) | ((p2_genvars["POS"][:] > ace_dupe - flanking_bp) & (p2_genvars["POS"][:] < ace_dupe + flanking_bp)) | (p2_genvars["POS"] == ace_119S)

# # ...or variants in the duplication?
# flanking_bp = 1e4
# p1_var_boo = ((p1_genvars["POS"][:] > ace_dups - flanking_bp) & (p1_genvars["POS"][:] < ace_dupe + flanking_bp)) | (p1_genvars["POS"] == ace_119S)
# p2_var_boo = ((p2_genvars["POS"][:] > ace_dups - flanking_bp) & (p2_genvars["POS"][:] < ace_dupe + flanking_bp)) | (p2_genvars["POS"] == ace_119S)

# # ...or variants in the downstream region?
# flanking_bp = 2e5
# p1_var_boo = ((p1_genvars["POS"][:] > ace_dupe) & (p1_genvars["POS"][:] < ace_dupe + flanking_bp )) | (p1_genvars["POS"] == ace_119S)
# p2_var_boo = ((p2_genvars["POS"][:] > ace_dupe) & (p2_genvars["POS"][:] < ace_dupe + flanking_bp )) | (p2_genvars["POS"] == ace_119S)

# # ...or variants in the upstream region?
# flanking_bp = 1e5
# p1_var_boo = ((p1_genvars["POS"][:] > ace_dups - flanking_bp) & (p1_genvars["POS"][:] < ace_dups )) | (p1_genvars["POS"] == ace_119S)
# p2_var_boo = ((p2_genvars["POS"][:] > ace_dups - flanking_bp) & (p2_genvars["POS"][:] < ace_dups )) | (p2_genvars["POS"] == ace_119S)

# compress
p1_genvars_sub = p1_genvars[:].compress(p1_var_boo)
p2_genvars_sub = p2_genvars[:].compress(p2_var_boo)

# expand to haplotyes
p1_genotyp_sub = p1_genotyp.compress(p1_var_boo, axis=0)
p1_haploty_sub = p1_genotyp_sub.to_haplotypes()
p2_genotyp_sub = p2_genotyp.compress(p2_var_boo, axis=0)
p2_haploty_sub = p2_genotyp_sub.to_haplotypes()


# restrict to populations of interest
p1_sam_ixs = p1_popdich["all"]
p2_sam_ixs = p2_popdich["all"]

# restrict to samples of interest
p1_sam_ixs = np.where(np.isin(p1_sampleh["ox_code"], sample_119S_list))[0]
p2_sam_ixs = np.where(np.isin(p2_sampleh["ox_code_sam"], sample_119S_list))[0]

p1_haploty_sub = np.take(a=p1_haploty_sub, indices = p1_sam_ixs, axis=1)
p2_haploty_sub = np.take(a=p2_haploty_sub, indices = p2_sam_ixs, axis=1)

# # get segregating from phase 2
# p2_hapalco_sub = p2_haploty_sub.count_alleles()
# p2_is_seg = p2_hapalco_sub.is_segregating()
# p2_haploty_sub_seg = p2_haploty_sub.compress(p2_is_seg, axis=0)

# # find problem haplotyes
# p2_haps_from_samples_with_mut = p2_sampleh[p2_sampleh["genotype"] != "wt"]["ox_code"].values
# p2_haps_from_samples_with_mut_boo = np.isin(p2_samlilh, p2_haps_from_samples_with_mut)

# find solution haplotypes
p1_focal_snp_ix = np.where(p1_genvars_sub["POS"] == ace_119S)[0][0]
p1_focal_snp_alt = p1_genvars_sub["ALT"][p1_focal_snp_ix].astype(str)
p1_focal_snp_ref = p1_genvars_sub["REF"][p1_focal_snp_ix].astype(str)

# lists of p1 haplotypes with REF and ALT alleles
p1_haps_with_alt_ix = np.where(p1_haploty_sub[p1_focal_snp_ix,] == 1)[0]
# p1_haps_with_alt = p1_samlilh[p1_haps_with_alt_ix]
# p1_haps_with_ref_ix = np.where(p1_haploty_sub[p1_focal_snp_ix,] == 0)[0]
# p1_haps_with_ref = p1_samlilh[p1_haps_with_ref_ix]

# p1 response array
p1_response = np.zeros(shape=p1_haploty_sub.shape[1]).astype(int)
p1_response[ p1_haps_with_alt_ix ] = 1
np.unique(p1_response, return_counts=True)

# homogeneise phase1 and phase2 (i.e. keep only phased variants present in both datasets)
p1_in_p2_boo = np.isin(p1_genvars_sub["POS"], p2_genvars_sub["POS"])
p2_in_p1_boo = np.isin(p2_genvars_sub["POS"], p1_genvars_sub["POS"])

# subsample p1 to fit p2
p1m_genvars_sub = p1_genvars_sub.compress(p1_in_p2_boo)
p1m_haploty_sub = p1_haploty_sub.compress(p1_in_p2_boo)

# subsample p2 to fit p1
p2m_genvars_sub = p2_genvars_sub.compress(p2_in_p1_boo)
p2m_haploty_sub = p2_haploty_sub.compress(p2_in_p1_boo)

# p1 train data
p1_dat = np.transpose(p1m_haploty_sub)

# p2 problem array
p2_dat = np.transpose(p2m_haploty_sub)

# supervised learning UMAP
import umap

# train UMAP embedding with test data
print("# UMAP")
# p1_mapper = umap.UMAP(n_neighbors=20, n_components=2, metric="hamming").fit(p1_dat, np.array(p1_response))
p1_mapper = umap.UMAP(n_neighbors=10, n_components=2, metric="hamming", target_weight=1.0).fit(p1_dat, np.array(p1_response))
p1_embedding = p1_mapper.embedding_

# map test data with the trained embedder
p2_embedding = p1_mapper.transform(p2_dat)

# save pdf
plt.figure(figsize=(5,5))
plt.scatter(p2_embedding[:,0], p2_embedding[:,1], s=2, color="gray")
plt.scatter(p1_embedding[:,0], p1_embedding[:,1], c=np.array(p1_response), s=4, cmap="bwr")
plt.title("UMAP (n var = %i)" % p1_dat.shape[1])
plt.savefig("%s/umap_classification.pdf" % (results_fo), bbox_inches='tight')
plt.close()

# train classifiers on UMAP coordinates
import sklearn
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier

# train decision tree
print("# Decision tree")
dtc = DecisionTreeClassifier(max_depth=None, min_samples_split=10, class_weight="balanced", criterion="entropy")
dtc.fit(p1_embedding, p1_response)
p2_dtc_prediction = dtc.predict(p2_embedding)

print("# accuracy score =", accuracy_score(p1_response, dtc.predict(p1_embedding)))
print("# confusion matrix")
print(confusion_matrix(p1_response, dtc.predict(p1_embedding)))
print("# classification report")
print(classification_report(p1_response, dtc.predict(p1_embedding)))
print("# prediction counts")
print(np.unique(p2_dtc_prediction, return_counts=True))

# try the same with a more complex random forest
print("# Random forest")
rfc = RandomForestClassifier(max_depth=None, min_samples_split=10, n_estimators=100, class_weight="balanced", bootstrap=True, criterion="entropy")
rfc.fit(p1_embedding, p1_response)
p2_rfc_prediction = rfc.predict(p2_embedding)

print("# accuracy score =", accuracy_score(p1_response, rfc.predict(p1_embedding)))
print("# confusion matrix")
print(confusion_matrix(p1_response, rfc.predict(p1_embedding)))
print(classification_report(p1_response, rfc.predict(p1_embedding)))
print("# prediction counts")
print(np.unique(p2_rfc_prediction, return_counts=True))

# try a random forest on the original data (prone to over-fitting)
print("# Random forest on original data")
rfcd = RandomForestClassifier(max_depth=None, min_samples_split=10, n_estimators=1000, class_weight="balanced", bootstrap=False, criterion="gini")
rfcd.fit(p1_dat, p1_response)
p2_rfcd_prediction = rfcd.predict(p2_dat)

print("# accuracy score =", accuracy_score(p1_response, rfcd.predict(p1_dat)))
print("# confusion matrix")
print(confusion_matrix(p1_response, rfcd.predict(p1_dat)))
print(classification_report(p1_response, rfcd.predict(p1_dat)))
print("# prediction counts")
print(np.unique(p2_rfcd_prediction, return_counts=True))

# store UMAP-based predictions
p2_prediction = pd.DataFrame()
p2_prediction["ox_code"] = p2_samlilh
p2_prediction["p2_dtc_prediction"] = 0
p2_prediction["p2_dtc_prediction"].loc[p2_sam_ixs] = p2_dtc_prediction
p2_prediction["p2_rfc_prediction"] = 0
p2_prediction["p2_rfc_prediction"].loc[p2_sam_ixs] = p2_rfc_prediction
p2_prediction["p2_rfcd_prediction"] = 0
p2_prediction["p2_rfcd_prediction"].loc[p2_sam_ixs] = p2_rfcd_prediction

# add p2 embeddings
for i in range(p2_embedding.shape[1]):
	p2_prediction["p2_umap_%i" % i] = np.inf
	p2_prediction["p2_umap_%i" % i].loc[p2_sam_ixs] = p2_embedding[:,i]

# print
p2_prediction.to_csv("%s/umap_classification.csv" % (results_fo), index=False, sep="\t")


# p2_prediction["p2_rfc_prediction"] = p2_rfc_prediction
# p2_prediction["p2_rfc_prediction_from_vars"] = p2_rfcd_prediction

# # add phase1 calls
# # MEANINGLESS BECAUSE PHASING IS INDEPENDENT AND THEREFORE ab LABELS CANNOT BE COMPARED
# p2_prediction["p1_calls"] = -1
# p2_prediction["p1_calls"] [ np.isin(p2_prediction["ox_code"] , p1_haps_with_ref) ] = 0
# p2_prediction["p1_calls"] [ np.isin(p2_prediction["ox_code"] , p1_haps_with_alt) ] = 1

# p2_prediction.to_csv("%s/umap_classification.csv" % (results_fo), index=False, sep="\t")


# # subset p2 to keep only haplotypes from 
# p2_response = np.zeros(shape=p2_samlilh.shape).astype(int)
# p2_response[:] = -1
# p2_response [ np.isin(element=p2_samlilh, test_elements=p2_haps_from_samples_with_mut) ] = 0
# p2_response [ np.isin(element=p2_samlilh, test_elements=p1_haps_with_mut) ] = 1
# np.unique(p2_response, return_counts=True)

# # training data
# dat_is_train = p2_response != -1
# dat_train = p2_haploty_sub_seg.compress(dat_is_train, axis=1)
# dat_train = np.transpose(dat_train)
# dat_train_response = p2_response[dat_is_train]

# # test data
# dat_is_test = p2_response == -1
# dat_test = p2_haploty_sub_seg.compress(dat_is_test, axis=1)
# dat_test = np.transpose(dat_test)
# dat_test_response = p2_response[dat_is_test]


# # subsample
# vars_random = np.random.choice(p2_haploty_sub_seg.shape[0], size=10000)
# vars_input  = np.transpose(p2_haploty_sub_seg[:][vars_random,:])

# # normal
# vars_input = np.transpose(p2_haploty_sub_seg)
# map_embedding = umap.UMAP().fit_transform(vars_input, y=p2_response)

# # plot embedding
# fig, ax = plt.subplots(1)
# plt.scatter(*map_embedding.T, c=p2_response, alpha=0.5, s=1, cmap='Spectral')
# cbar = plt.colorbar(boundaries=np.arange(4))
# cbar.set_ticklabels([-1,0,1])
# plt.show()





