#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
import allel
import h5py
import zarr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import itertools
from scipy.spatial.distance import squareform
# load plot settings
sns.set(context="notebook",style="ticks",
        font_scale=1,font="Arial",palette="bright")

# # Haplotype analysis for genes in high-$PBS$ regions
# ## Input
# Input files and parameters:
# output
outdir   = "results_peaks" # where to store output
popc     = "population"

# input data phase1
oc_metasam_fn = "../metadata/samples.meta_phenotypes_acegenotype.simple.txt"
oc_hapcall_fn = "/home/xavi/Documents/VariationAg1k/data/phase2.AR1/haplotypes/zarr2/ag1000g.phase2.ar1.samples/"
oc_accessi_fn = "/home/xavi/Documents/VariationAg1k/data/phase2.AR1/accessibility/accessibility.h5"
oc_popc       = popc
oc_popl       = ["CIcol"]
# gff
gffann_fn  = "../metadata/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3"
genes_fn = "differentiation_output_post.PBStop_genes.csv"
genes = pd.read_csv(genes_fn, sep="\t")
genes = genes[["chrom","window_s","window_e"]]
genes = genes.drop_duplicates(ignore_index=True)

# manually add the odd kmer region
add_kmer_region =  pd.DataFrame([["orphan_kmer","2L",8662023,8662023+1e4]], columns=["gene","chrom","window_s","window_e"])
genes = genes.append(add_kmer_region, ignore_index=True)

# gene list
genelist = genes["window_s"].values

prev_chrom = "any"
for n,gene in enumerate(genelist):

	# gene of interest
	chrom     = genes.iloc[n]["chrom"]
	l_nom     = gene # nom loci
	loc_start = genes.iloc[n]["window_s"]     # start gene
	loc_end   = genes.iloc[n]["window_e"]       # end gene

	if prev_chrom != chrom:

		# ## Load data
		# ### Phase1 variants
		# Population and sample data:
		# load samples list with sample code, groupings, locations etc.
		oc_samples_df   = pd.read_csv(oc_metasam_fn, sep='\t')
		oc_samples_bool = (oc_samples_df[oc_popc].isin(oc_popl).values)
		oc_samples      = oc_samples_df[oc_samples_bool]
		oc_samples.reset_index(drop=True, inplace=True)
		# indexed dictionary of populations
		oc_popdict = dict()
		for popi in oc_popl: 
			oc_popdict[popi]  = oc_samples[oc_samples[oc_popc] == popi].index.tolist()
		# add an extra population composed of all other locations
		oc_popdict["all"] = []
		for popi in oc_popl:
			oc_popdict["all"] = oc_popdict["all"] + oc_popdict[popi]
		# report
		print("Data:")
		print("* Samples     = ", oc_samples.shape[0])
		print("* Populations = ", set(oc_samples[oc_popc]))
		print(oc_samples.groupby(("population")).size())
		# Phased variants and genotypes:
		# declare objects with variant data
		oc_hapcall   = zarr.open(oc_hapcall_fn)
		oc_hapcall_goodsamples = oc_hapcall["X"]["samples"][:].astype(str)
		oc_hapcall_goodsamples_ix = np.where(np.isin(element=oc_hapcall_goodsamples , test_elements=oc_samples["ox_code"].values))[0]
		# variants of genotypes
		print("Variants phased...")
		oc_hapcall_var = oc_hapcall[chrom]["variants"]
		oc_hapvars = allel.VariantChunkedTable(oc_hapcall_var,names=["POS","REF","ALT"],index="POS") 
		print(oc_hapvars.shape)
		# genotype data
		print("Genotypes phased...")
		oc_hapcall_hap = oc_hapcall[chrom]["calldata"]["genotype"]
		oc_haploty     = allel.GenotypeChunkedArray(oc_hapcall_hap) 
		oc_haploty     = oc_haploty.subset(sel1=oc_hapcall_goodsamples_ix)
		print(oc_haploty.shape)
		# Get haplotypes from phased variants:

		# recast haplotypes: drop ploidy
		print("Expand phase haplotypes...")
		oc_haploty_hap = oc_haploty.to_haplotypes()
		print(oc_haploty_hap.shape)
		# ### Sample data
		# Get dataframe from metadata file, with sample codes, species and populations:

		oc_samples = pd.DataFrame(data={
			"ox_code"    :  oc_samples["ox_code"].values.tolist() ,
			"species"    :  oc_samples["m_s"].values.astype(str).tolist(),
			"population" :  oc_samples[oc_popc].values.tolist() ,
			"sex" :  oc_samples["sex"].values.tolist() ,
			"phenotype"  :  oc_samples["phenotype"].values.tolist()
		})
		print(oc_samples.shape)
		# rename species...
		oc_samples["species"].values[oc_samples["species"].values == "M"]   = "col"
		oc_samples["species"].values[oc_samples["species"].values == "S"]   = "gam"
		oc_samples["species"].values[oc_samples["species"].values == "M/S"] = "gamcol"
		oc_samples["species"].values[oc_samples["species"].values == "M-S"] = "gamcol"
		oc_samples["species"].values[oc_samples["species"].values == "nan"] = "gamcol"
		# obtain population & species list
		oc_popl = np.unique(oc_samples["population"].values)
		oc_spsl = np.unique(oc_samples["species"].values)
		# Duplicate rows in metadata dataframe, to get population metadata from each haplotype:

		oc_sampleh = pd.DataFrame(data={
			"ox_code"    :  list(itertools.chain(*[[ s + 'a', s + 'b'] for s in oc_samples["ox_code"].values.tolist()])),    # takes col from oc_samples and duplicates it, a/b
			"species"    :  list(itertools.chain(*[[ s      , s      ] for s in oc_samples["species"].values.tolist()])),
			"population" :  list(itertools.chain(*[[ s      , s      ] for s in oc_samples["population"].values.tolist()])),
			"sex" :  list(itertools.chain(*[[ s      , s      ] for s in oc_samples["sex"].values.tolist()])),
			"phenotype"  :  list(itertools.chain(*[[ s      , s      ] for s in oc_samples["phenotype"].values.tolist()]))
		})
		print(oc_sampleh.shape)

		print("Population dict...")
		oc_popdict = dict()
		oc_popdict["CIcol"] = oc_samples[oc_samples["population"] == popi].index.tolist()
		print("Population dict phased...")
		oc_popdich = dict()
		oc_popdich["CIcol"] = oc_sampleh[oc_sampleh["population"] == popi].index.tolist()
		# ### Allele counts
		# Using both dictionaries:

		print("Genotypes phased to allele counts (population)...")
		oc_hapalco_pop = oc_haploty.count_alleles_subpops(subpops=oc_popdict)
		print(oc_hapalco_pop.shape)
		print("Haplotypes phased to allele counts (population)...")
		oc_hapalco_hap_pop = oc_haploty_hap.count_alleles_subpops(subpops=oc_popdich)
		print(oc_hapalco_hap_pop.shape)
		# ### Filters
		# #### Retain segregating and non-singletons
		# Define which phased variants to retain from phase1:

		# subset data: segregating alleles & no singletons
		print("Filters phased...")
		oc_is_seg_h    = oc_hapalco_hap_pop["CIcol"].is_segregating()[:] # segregating
		oc_is_nosing_h = oc_hapalco_hap_pop["CIcol"][:,:2].min(axis=1)>1 # no singletons
		# subset phase2 to segregating & no singletons
		oc_hapvars_seg         = oc_hapvars.compress((oc_is_seg_h & oc_is_nosing_h))
		oc_haploty_seg         = oc_haploty.compress((oc_is_seg_h & oc_is_nosing_h))
		oc_hapalco_pop_seg     = oc_hapalco_pop.compress((oc_is_seg_h & oc_is_nosing_h))
		oc_haploty_hap_seg     = oc_haploty_hap.compress((oc_is_seg_h & oc_is_nosing_h))
		oc_hapalco_hap_pop_seg = oc_hapalco_hap_pop.compress((oc_is_seg_h & oc_is_nosing_h))
		# report
		print(oc_haploty_seg.shape,"/", oc_haploty.shape)
		# ### Other data
		# Accessibility:

		# Accessibility
		print("Load accessibility array...")
		accessi_df  = h5py.File(oc_accessi_fn,mode="r")
		accessi_arr = accessi_df[chrom]["is_accessible"][:]
		# ## Selection signals in clusters
		# We want to see if the resistant individuals have positive selection. 
		# Create dictionary:

		popdich_clu = dict()
		# clusters, and non clustered
		popdich_clu["alive"] = np.where(oc_sampleh["phenotype"] == "alive")[0]
		popdich_clu["dead"]  = np.where(oc_sampleh["phenotype"] == "dead")[0]
		popdich_clu["all"]  =  np.where(oc_sampleh["population"] == "CIcol")[0]
		# allele counts in these clusters:
		oc_hapalco_hap_clu_seg = oc_haploty_hap_seg.count_alleles_subpops(subpops=popdich_clu)
		oc_hapalco_hap_clu_seg.shape

	# Colors for plot:
	# list colors for each haplotype in popdich_clu.keys()
	colors = ["forestgreen","deeppink","lightsteelblue"]


	# ### EHH decay
	# Now calculate **EHH decay** on the region of interest, using phased variants around various variants.
	# Common parameters for EHH plots:
	ehh_above_thr = 0.50
	ehh_below_thr = 0.05
	flank_bp_EHH  = 2e5

	# variants to retain
	clu_varbool_up = np.logical_and(oc_hapvars_seg["POS"] >= loc_start-flank_bp_EHH, oc_hapvars_seg["POS"] < loc_start)
	clu_varbool_do = np.logical_and(oc_hapvars_seg["POS"] > loc_end, oc_hapvars_seg["POS"] <= loc_end+flank_bp_EHH)
	clu_varbool    = np.logical_or(clu_varbool_up,clu_varbool_do)
	# samples to remove from analysis (EHH function can't handle missing -1 data)
	rmv_miss_ix   = np.unique(np.where(oc_haploty_hap_seg.subset(sel0=clu_varbool) == -1)[1]).tolist()
	rmv_miss_bool = np.invert(np.isin(range(0,oc_haploty_hap_seg.n_haplotypes),test_elements=rmv_miss_ix))
	# positions
	clu_ehh_pos = oc_hapvars_seg["POS"].subset(sel0=clu_varbool)
	# plot
	pdf = PdfPages("%s/sel_%s_%s_EHHdecay.pdf" % (outdir,chrom,l_nom))
	fig = plt.figure(figsize=(5,3))
	ax3 = plt.subplot(1, 1, 1)
	for i,clu_key in enumerate(popdich_clu.keys()):
		print("EHH %s" % clu_key)
		# which variants include in the cluster-wise analysis of selection?
		clu_sambool = np.isin(range(0,oc_haploty_hap_seg.n_haplotypes),test_elements=popdich_clu[clu_key])
		clu_sambool = np.logical_and(clu_sambool,rmv_miss_bool)
		
		# calculate actual EHH
		clu_ehh_up_i = allel.ehh_decay(h=oc_haploty_hap_seg.subset(sel0=clu_varbool_up,sel1=clu_sambool))
		clu_ehh_do_i = allel.ehh_decay(h=oc_haploty_hap_seg.subset(sel0=clu_varbool_do,sel1=clu_sambool))
		clu_ehh_i    = np.concatenate((clu_ehh_up_i[::-1],clu_ehh_do_i))
		clu_ehh_i_ar = np.trapz(clu_ehh_i)
		ehh_above_start = clu_ehh_pos.compress(clu_ehh_i > ehh_above_thr)[0]
		ehh_above_end   = clu_ehh_pos.compress(clu_ehh_i > ehh_above_thr)[-1]
		ehh_below_start = clu_ehh_pos.compress(clu_ehh_i < ehh_below_thr)[0]
		ehh_below_end   = clu_ehh_pos.compress(clu_ehh_i < ehh_below_thr)[-1]
		# lab is data
		clu_lab    = "%s, n=%i, a=%.3f\nEHH>%.2f: %i bp %i-%i\nEHH<%.2f: %i bp %i-%i" % (
			clu_key, len(popdich_clu[clu_key]),clu_ehh_i_ar, 
			ehh_above_thr, ehh_above_end-ehh_above_start, ehh_above_start, ehh_above_end,
			ehh_below_thr, ehh_below_end-ehh_below_start, ehh_below_start, ehh_below_end
		)
		
		# plot EHH background & foreground
		ax3.plot(clu_ehh_pos/1e6,clu_ehh_i,color=colors[i],label=clu_lab,mfc='none')
	sns.despine(ax=ax3,offset=10)
	ax3.set_title("EHH decay %s, %s:%i-%i +/- %i, n=%s vars" % (l_nom,chrom,loc_start,loc_end,flank_bp_EHH,clu_ehh_pos.shape[0]))
	ax3.set_xlabel("Mb")
	ax3.set_ylabel("EHH")
	ax3.set_ylim(0,1)
	plt.axhline(ehh_above_thr, color='lightgray',linestyle=":",label=str(ehh_above_thr))
	plt.axhline(ehh_below_thr, color='lightgray',linestyle=":",label=str(ehh_below_thr))
	plt.axvline(loc_start/1e6, color='magenta',linestyle=":",label="gene")
	plt.axvline(loc_end/1e6, color='magenta',linestyle=":",label="")
	# ax3.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
	# save
	pdf.savefig(fig,bbox_inches='tight')
	pdf.close()


	# ### Garud H and haplotype diversity
	# Compute **Garud H statistics and haplotype diversity** for each cluster and estimates in the region of interest. Plots represent a wide region around the gene, and statistics are estimated from variants within the cluster only.
	# variants to examine
	block_size = 100
	# region to plot
	flanking_bp = 5e5
	clu_varbool = np.logical_and(
		oc_hapvars_seg["POS"] >= loc_start-flanking_bp,
		oc_hapvars_seg["POS"] <= loc_end+flanking_bp)
	# region to focus: statistics will be calcualted in this region
	clu_varbool_focus = np.logical_and(oc_hapvars_seg["POS"] > loc_start-1e4, oc_hapvars_seg["POS"] <= loc_end+1e4)


	# First, for H12 plot:
	# open PDF
	pdf = PdfPages("%s/sel_%s_%s_GarudH12.pdf" % (outdir,chrom,l_nom))
	fig = plt.figure(figsize=(4,3))
	ax9 = plt.subplot(1, 1, 1)
	for i,clu_key in enumerate(popdich_clu.keys()):
		# which variants include in the cluster-wise analysis of selection?
		clu_sambool = np.isin(range(0,oc_haploty_hap_seg.n_haplotypes),test_elements=popdich_clu[clu_key])
		clu_sambool = np.logical_and(clu_sambool,rmv_miss_bool)
		# Garud H along chromosome
		clu_pos_wib = allel.moving_statistic(
			oc_hapvars_seg["POS"].subset(sel0=clu_varbool), statistic=lambda v: v[0], size=block_size)
		clu_gah_wib = allel.moving_garud_h(
			oc_haploty_hap_seg.subset(sel0=clu_varbool,sel1=clu_sambool), size=block_size)
		
		# garud in focus region
		gah_focus_est = allel.moving_garud_h(oc_haploty_hap_seg.subset(sel0=clu_varbool_focus, sel1=clu_sambool), size=block_size)
		gah_focus_est_jack = allel.stats.misc.jackknife(gah_focus_est[1], statistic=np.nanmean)
		clu_label = "%s\nH12 = %.6f +/- %.6f SE, n = %i" % (clu_key, gah_focus_est_jack[0], gah_focus_est_jack[1],np.sum(clu_sambool))
		print(clu_label)
		# plot
		plt.step(clu_pos_wib/1e6, clu_gah_wib[1], color=colors[i], label=clu_label)
		
	sns.despine(ax=ax9,offset=10)
	ax9.set_title("Garud H12 %s:%i-%i %s" % (chrom,loc_start,loc_end,l_nom))
	ax9.set_ylim(0,1)
	ax9.set_xlabel("Mb")
	ax9.set_ylabel("h")
	plt.axvline(loc_start/1e6, color='magenta',linestyle=":",label="gene")
	plt.axvline(loc_end/1e6, color='magenta',linestyle=":",label="")
	# ax9.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
	# save
	pdf.savefig(fig,bbox_inches='tight')
	pdf.close()


	# Now, same with H2/H1 plot:
	# open PDF
	pdf = PdfPages("%s/sel_%s_%s_GarudH2H1.pdf" % (outdir,chrom,l_nom))
	fig = plt.figure(figsize=(4,3))
	ax9 = plt.subplot(1, 1, 1)
	for i,clu_key in enumerate(popdich_clu.keys()):
		# which variants include in the cluster-wise analysis of selection?
		clu_sambool = np.isin(range(0,oc_haploty_hap_seg.n_haplotypes),test_elements=popdich_clu[clu_key])
		clu_sambool = np.logical_and(clu_sambool,rmv_miss_bool)
		# Garud H along chromosome
		clu_pos_wib = allel.moving_statistic(
			oc_hapvars_seg["POS"].subset(sel0=clu_varbool), statistic=lambda v: v[0], size=block_size)
		clu_gah_wib = allel.moving_garud_h(
			oc_haploty_hap_seg.subset(sel0=clu_varbool,sel1=clu_sambool), size=block_size)
		
		# garud in focus region
		gah_focus_est = allel.moving_garud_h(oc_haploty_hap_seg.subset(sel0=clu_varbool_focus, sel1=clu_sambool), size=block_size)
		gah_focus_est_jack = allel.stats.misc.jackknife(gah_focus_est[3], statistic=np.nanmean)
		clu_label = "%s\nH2H1 = %.6f +/- %.6f SE, n = %i" % (clu_key, gah_focus_est_jack[0], gah_focus_est_jack[1],np.sum(clu_sambool))
		print(clu_label)
		# plot
		plt.subplot(1, 1, 1)
		plt.step(clu_pos_wib/1e6, clu_gah_wib[3], color=colors[i], label=clu_label)
		
	sns.despine(ax=ax9,offset=10)
	ax9.set_title("Garud H2H1 %s:%i-%i %s" % (chrom,loc_start,loc_end,l_nom))
	ax9.set_ylim(0,1)
	ax9.set_xlabel("Mb")
	ax9.set_ylabel("H2H1")
	plt.axvline(loc_start/1e6, color='magenta',linestyle=":",label="gene")
	plt.axvline(loc_end/1e6, color='magenta',linestyle=":",label="")
	# ax9.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
	# save
	pdf.savefig(fig,bbox_inches='tight')
	pdf.close()


	# Finally, haplotype diversity:
	# open PDF
	pdf = PdfPages("%s/sel_%s_%s_hapdiv.pdf" % (outdir,chrom,l_nom))
	fig = plt.figure(figsize=(4,3))
	ax9 = plt.subplot(1, 1, 1)
	for i,clu_key in enumerate(popdich_clu.keys()):
		# which variants include in the cluster-wise analysis of selection?
		clu_sambool = np.isin(range(0,oc_haploty_hap_seg.n_haplotypes),test_elements=popdich_clu[clu_key])
		clu_sambool = np.logical_and(clu_sambool,rmv_miss_bool)
		# hap div along chromosome
		clu_pos_wib = allel.moving_statistic(
			oc_hapvars_seg["POS"].subset(sel0=clu_varbool), statistic=lambda v: v[0], size=block_size)
		clu_hdi_wib = allel.moving_haplotype_diversity(
			oc_haploty_hap_seg.subset(sel0=clu_varbool,sel1=clu_sambool), size=block_size)
		
		# garud in focus region
		gah_focus_est = allel.moving_haplotype_diversity(oc_haploty_hap_seg.subset(sel0=clu_varbool_focus, sel1=clu_sambool), size=block_size)
		gah_focus_est_jack = allel.stats.misc.jackknife(gah_focus_est, statistic=np.nanmean)
		clu_label = "%s\nh = %.6f +/- %.6f SE, n = %i" % (clu_key, gah_focus_est_jack[0], gah_focus_est_jack[1],np.sum(clu_sambool))
		print(clu_label)
		# plot
		plt.subplot(1, 1, 1)
		plt.step(clu_pos_wib/1e6, clu_hdi_wib, color=colors[i], label=clu_label)
		
	sns.despine(ax=ax9,offset=10)
	ax9.set_title("Hap diversity %s:%i-%i %s" % (chrom,loc_start,loc_end,l_nom))
	ax9.set_ylim(0,1)
	ax9.set_xlabel("Mb")
	ax9.set_ylabel("h")
	plt.axvline(loc_start/1e6, color='magenta',linestyle=":",label="gene")
	plt.axvline(loc_end/1e6, color='magenta',linestyle=":",label="")
	# ax9.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
	# save
	pdf.savefig(fig,bbox_inches='tight')
	pdf.close()



	# export fasta

	# refs and alts
	export_start = loc_start
	export_end = loc_end
	hap_bool    = np.logical_and(oc_hapvars_seg["POS"][:] >= export_start, oc_hapvars_seg["POS"][:] <= export_end)
	hapvars_sub = oc_hapvars_seg.compress(hap_bool)
	haploty_sub = oc_haploty_hap_seg.compress(hap_bool)
	hapvars_sub_REF = hapvars_sub["REF"][:].astype(str)
	hapvars_sub_ALT = hapvars_sub["ALT"][:].astype(str)

	# output
	print("FASTA...")
	happhy = pd.DataFrame({
		"hap": ">"+ oc_sampleh["ox_code"] + "_" + oc_sampleh["population"] + "_" + oc_sampleh["phenotype"],
		"seq": np.nan},    
		columns=["hap", "seq"])

	

	for pn,popi in enumerate(oc_sampleh["ox_code"]):
		
		popi_gen = np.ndarray.tolist(haploty_sub[:,pn])
		popi_seq = [hapvars_sub_REF[gn] if gei == 0 else hapvars_sub_ALT[gn] for gn,gei in enumerate(popi_gen)]
		happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

	happhy.to_csv("%s/sel_%s_%s_alignment.fasta" % (outdir,chrom,l_nom),sep="\n",index=False, header=False)


	prev_chrom = chrom

