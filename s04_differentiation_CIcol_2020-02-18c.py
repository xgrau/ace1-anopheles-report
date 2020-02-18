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
chromlist  = ["2R","2L","3R","3L","X"]
metasam_fn = "metadata/samples.meta_phenotypes.txt"
callset_fn = "/home/xavi/dades/Variation/phase2.AR1/variation/main/zarr2/ag1000g.phase2.ar1.pass/"
results_fo = "results_differentiation_CIcol"
anngene_fn = "metadata/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3"

# general settings
sns.set(style="ticks",
        font_scale=1.3,
		  rc={"lines.linewidth": 1},
        font="Arial")

logging.basicConfig(
	level=logging.DEBUG, 
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s"
	)

# define 0.001 threshld for standardised tests (two-sided)
zscore_thr = scipy.stats.norm.ppf((0.001)/2)


## LOAD METADATA
# population metadata
samples = pd.read_csv(metasam_fn, sep='\t')

# sample: haplotypes
sampleh = pd.DataFrame(data = {
	"haplotype" : np.array(list(itertools.chain(*[[s + ' a', s + ' b'] for s in samples["ox_code"] ]))),
	"population": np.array(list(itertools.chain(*[[s , s] for s in samples["population"] ]))),
	"phenotype": np.array(list(itertools.chain(*[[s , s] for s in samples["phenotype"] ])))
})

# load gene annotations: in pandas format, and in pyranges format
# pandas
anngene_pd     = allel.gff3_to_dataframe(path=anngene_fn, attributes=["ID","Parent"])
anngene_pd_gen = anngene_pd[anngene_pd["type"] == "gene"]

## POPULATION DICTIONARIES
# population dictionary
popdict = dict()
popdict["alive"] = samples[  np.logical_and( samples["phenotype"] == "alive", samples["population"] == "CIcol"  )  ].index.tolist()
popdict["dead"]  = samples[  np.logical_and( samples["phenotype"] == "dead",  samples["population"] == "CIcol"  )  ].index.tolist()
popdict["AOcol"] = samples[  samples["population"] == "AOcol"  ].index.tolist()
popdict["CIcol"] = popdict["alive"] + popdict["dead"]
popdict["all"]   = popdict["CIcol"] + popdict["AOcol"]

## LOAD DATA

def load_genotypes(chrom, zarr_fn, popdict, genotype_tag="GT", segregating_pop="all"):
	# load callset
	callset = zarr.open(zarr_fn)
	genvars = allel.VariantChunkedTable( callset["variants"], names=["POS","REF","ALT"],index="POS" )
	genotyp = allel.GenotypeChunkedArray(callset["calldata"][genotype_tag])
	# calculate allele counts
	genalco = genotyp.count_alleles_subpops(subpops=popdict)
	# filters: segregating, no singletons
	is_seg    = genalco[segregating_pop].is_segregating()[:]   # segregating in population of interest? (default: any pop aka "all")
	is_nosing = genalco[segregating_pop][:,:2].min(axis=1) > 1 # no singletons
	is_retain = np.logical_and(is_seg, is_nosing)
	# apply filters
	genotyp_seg = genotyp.subset(sel0=is_retain)
	genalco_seg = genalco.compress(is_retain)
	genvars_seg = genvars[:].compress(is_retain)
	# logging
	logging.info("chr %s seg variants: %i / %i ( %.1f )" % ( chrom, genotyp_seg.shape[0], genotyp.shape[0], 100 * genotyp_seg.shape[0] / genotyp.shape[0]) )
	logging.info("chr %s samples: %i" % (chrom,genotyp_seg.shape[1]))
	# return
	return genotyp_seg, genalco_seg, genvars_seg

# a function to detect whether the statistic exceeds a certain threshold (works with intervals and per-variant data, 
# and with  greather-than, less-than, equal, etc. comparisons)
def statistic_above_threshold(stat_array, index_array_start, thr, index_array_stop=None, compare="gt"):
	
	if index_array_stop is None:
		index_array_stop = index_array_start

	# select windows where the statistic exceeds the threshold
	if compare == "gt" :
		ix_thr = np.where(stat_array > thr)[0]
	if compare == "ge" :
		ix_thr = np.where(stat_array >= thr)[0]
	if compare == "lt" :
		ix_thr = np.where(stat_array < thr)[0]
	if compare == "le" :
		ix_thr = np.where(stat_array <= thr)[0]

	# find start and stop of select windows (in index array, i.e. positions)
	index_start = index_array_start[ix_thr]
	index_stop  = index_array_stop[ix_thr]
	
	# statistic at this position
	stat_thr = stat_array[ix_thr]

	# return: statistic at the selected window, start and end of window, indexes in the original array
	return stat_thr, index_start, index_stop, ix_thr


# empty dicts for phased variants
genvars_seg   = dict()
genotyp_seg   = dict()
genalco_seg   = dict()

# loop to load data per chromosome
for chrom in chromlist:
	
	# load phased variants
	genotyp_seg[chrom], genalco_seg[chrom], genvars_seg[chrom] = load_genotypes(
		chrom=chrom, 
		zarr_fn="%s/%s" % (callset_fn, chrom), 
		popdict=popdict, 
		genotype_tag="genotype",
		segregating_pop="CIcol")


#### Fst ####
var_block = 1000
var_step  = 1000

def fraction_above(array, thr=0.05):
	l = array.shape[0]
	a = np.sum(array > thr)
	f = a/l
	return f


## Fst differentiation dead~alive
# loop per chromosome
plt.figure(figsize=(24,2))
plt.subplots_adjust(wspace=0.3,hspace=0.6)
fst_b_pos = dict()
fst_b_est = dict()
for i,chrom in enumerate(chromlist):

	# Hudson Fst, block-wise
	fst_b_pos[chrom] = allel.moving_statistic(genvars_seg[chrom]["POS"], np.min, size=var_block, step=var_step)
	fst_b_est[chrom] = allel.moving_hudson_fst(genalco_seg[chrom]["alive"], genalco_seg[chrom]["dead"], size=var_block, step=var_step)

	# average hudson Fst
	fst_per_chrom = allel.average_hudson_fst(ac1=genalco_seg[chrom]["alive"], ac2=genalco_seg[chrom]["dead"], blen=var_block)

	# fraction of windows above Fst = 5%
	fst_b_est_fraction_above = fraction_above(fst_b_est[chrom], thr=0.05)

	# plot
	ax=plt.subplot(1,5,i+1)
	sns.despine(ax=ax,offset=5)
	plt.step(fst_b_pos[chrom]/1e6, fst_b_est[chrom], color="blue", where="pre")
	plt.title("%s\nFst=%.2E +/- %.2E\nFr win Fst>0.05=%.2E" % (chrom, fst_per_chrom[0], fst_per_chrom[1], fst_b_est_fraction_above))
	plt.axhline(y=0, color='black', linestyle='--')
	plt.xlabel("Mb")
	plt.ylabel("Fst")
	plt.xlim(0,65)
	plt.ylim(-0.05,1)

	### TODO: ADD TRACK WITH EXTREME VALUE DENSITY
	### TODO: OUTPUT TOP VALUES

# save
plt.savefig("%s/differentiation_Fst.pdf" % (results_fo), bbox_inches='tight')
plt.close()


## PBS differentiation dead~alive~AOcol
# loop per chromosome
plt.figure(figsize=(24,2))
plt.subplots_adjust(wspace=0.3,hspace=0.6)
pbs_b_pos = dict()
pbs_b_est = dict()
pbs_b_est_std = dict()
for i,chrom in enumerate(chromlist):

	# PBS, block-wise
	pbs_b_pos[chrom] = allel.moving_statistic(genvars_seg[chrom]["POS"], np.min, size=var_block, step=var_step)
	pbs_b_est[chrom] = allel.pbs(
		ac1=genalco_seg[chrom]["alive"], 
		ac2=genalco_seg[chrom]["dead"], 
		ac3=genalco_seg[chrom]["AOcol"], 
		window_size=var_block, 
		window_step=var_step,
		normed=True)

	# standardise to unit variance
	pbs_b_est_std[chrom] = allel.standardize(pbs_b_est[chrom])

	# find top windows
	pbs_b_est_top = statistic_above_threshold(stat_array=pbs_b_est[chrom], index_array_start=pbs_b_pos[chrom], thr=0.05, compare="gt")

	# plot
	ax=plt.subplot(1,5,i+1)
	sns.despine(ax=ax,offset=5)
	plt.step(pbs_b_pos[chrom]/1e6, pbs_b_est[chrom], color="blue", where="pre")
	plt.plot(pbs_b_pos[chrom][pbs_b_est_top[3]]/1e6, np.repeat(0.19, pbs_b_est_top[3].shape), color="red", marker="v", linewidth=0, markersize=5)
	plt.title("%s" % chrom)
	plt.axhline(y=0, color='black', linestyle='--')
	plt.xlabel("Mb")
	plt.ylabel("PBS")
	plt.xlim(0,65)
	plt.ylim(-0.05,0.2)

# save
plt.savefig("%s/differentiation_PBS_alive-dead-AOcol.pdf" % (results_fo), bbox_inches='tight')
plt.close()



# store results table
df_out = pd.DataFrame()
for i,chrom in enumerate(chromlist):

	df_start = allel.moving_statistic(genvars_seg[chrom]["POS"], np.min, size=var_block, step=var_step)
	df_stop = allel.moving_statistic(genvars_seg[chrom]["POS"], np.max, size=var_block, step=var_step)
	df_oui = pd.DataFrame({
		"chrom" : chrom,
		"start": df_start,
		"end": df_stop,
		"Fst": fst_b_est[chrom],
		"PBS": pbs_b_est[chrom],
		"PBS_s": pbs_b_est_std[chrom],
		"PBS_p": scipy.stats.norm.sf(abs(pbs_b_est_std[chrom]))*2
	})
	df_out = pd.concat([df_out, df_oui])
	

# write output
df_out.to_csv("%s/differentiation_output.csv" % results_fo, sep="\t", index=False)



## PCA
# remove SNPs in linkage disequilibrium from a genotype array
def ld_prune(
	gn, size, step, 
	threshold=.1, n_iter=1):
	for i in range(n_iter):
		loc_unlinked = allel.locate_unlinked(gn,  size=size,  step=step,  threshold=threshold)
		n = np.count_nonzero(loc_unlinked)
		n_remove = gn.shape[0] - n
		print('# iteration', i+1, 'retaining', n,  'removing', n_remove, 'variants')
		gn = gn.compress(loc_unlinked, axis=0)
	return gn


genotyp_seg_nalt_ldp_all = np.empty((0,len(popdict["CIcol"])),dtype="int8")
for chrom in ["3R","3L"]:
	# concatenate all genotypes in all chromosomes
	genotyp_seg_nalt = genotyp_seg[chrom].subset(sel1=popdict["CIcol"]).to_n_alt()[:]
	genotyp_seg_nalt_ldp = ld_prune( genotyp_seg_nalt, size=500, step=200,threshold=0.1, n_iter=10 )
	genotyp_seg_nalt_ldp_all = np.concatenate( ( genotyp_seg_nalt_ldp_all, genotyp_seg_nalt_ldp ) )
	print(chrom, "variants:", genotyp_seg_nalt_ldp.shape)

# PCA of LD-pruned genotypes in all chromosomes
# PCA plot function
def plot_pca_coords(
	coords, model, pc1, pc2,
	pop1,pop2,
	colorpop1,colorpop2,whichpop1, whichpop2, 
	title="PCA",figsizex=6, figsizey=6):

	fig, ax = plt.subplots(figsize=(figsizex, figsizey))
	sns.despine(ax=ax, offset=10)
	x = coords[:, pc1]
	y = coords[:, pc2]
	ax.plot(x[whichpop1], y[whichpop1], marker='o', linestyle=' ', label=pop1, markersize=6, color=colorpop1, mfc='none')
	ax.plot(x[whichpop2], y[whichpop2], marker='o', linestyle=' ', label=pop2, markersize=6, color=colorpop2, mfc='none')
	ax.set_xlabel('PC%s (%.2f%%)' % (pc1+1,model.explained_variance_ratio_[pc1]*100))
	ax.set_ylabel('PC%s (%.2f%%)' % (pc2+1,model.explained_variance_ratio_[pc2]*100))
	ax.set_title(title)
	ax.legend(loc='best')
	return fig,ax


# PCA variance explained histogram
def plot_pcavarexp(
	pcamodel,
	color="blue",
	title="Variance explained",figsizex=5, figsizey=4):
	fig, ax = plt.subplots(figsize=(figsizex, figsizey))
	sns.despine(ax=ax, offset=10)
	y = 100 * pcamodel.explained_variance_ratio_
	x = np.arange(len(y))
	ax.set_xticks(x + .4)
	ax.set_xticklabels(x + 1)
	ax.bar(x, y,color=color)
	ax.set_xlabel('PC')
	ax.set_ylabel('Variance (%)')
	ax.set_title(title)
	return fig

# define sample lables for nice PCA
samples_subset = samples[samples["population"] == "CIcol"]
samples_subset = samples_subset.reset_index(drop=True)
ix_pop1_subset = np.where(samples_subset["phenotype"] == "alive")[0]
ix_pop2_subset = np.where(samples_subset["phenotype"] == "dead")[0]

# PCA1-2
pdf_pages = PdfPages("%s/pca.pdf" % results_fo)
pca_coo, pca_mod = allel.pca(genotyp_seg_nalt_ldp_all)
fig,ax = plot_pca_coords( coords=pca_coo,model=pca_mod,pc1=0,pc2=1,pop1="alive", pop2="dead", colorpop1="green", colorpop2="magenta", whichpop1=ix_pop1_subset, whichpop2=ix_pop2_subset, title="PCA %i SNPs, LD-pruned" % (len(genotyp_seg_nalt_ldp_all)) )
fig,ax = plot_pca_coords( coords=pca_coo,model=pca_mod,pc1=0,pc2=2,pop1="alive", pop2="dead", colorpop1="green", colorpop2="magenta", whichpop1=ix_pop1_subset, whichpop2=ix_pop2_subset, title="PCA %i SNPs, LD-pruned" % (len(genotyp_seg_nalt_ldp_all)) )
pdf_pages.savefig(fig,bbox_inches='tight')

# variance explained
fig = plot_pcavarexp(pcamodel=pca_mod,color="slategray",title="Variance explained")
pdf_pages.savefig(fig,bbox_inches='tight')

pdf_pages.close()


