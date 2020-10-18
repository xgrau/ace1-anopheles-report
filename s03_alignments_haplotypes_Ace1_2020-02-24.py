#!/usr/bin/env python

# # Haplotype alignments
# Alignment of haplotypes around the *Ace1* duplication for further phylogenetic analysis.

# ## Input
# input data
outdir      = "results_admixture_phylo/"
metasam_fn  = "metadata/samples.meta_phenotypes_acegenotype.simple.txt"
callset_fn  = "/media/xavi/Saigon/VariationAg1k/data/phase2.AR1/variation/main/zarr2/ag1000g.phase2.ar1.pass/"
accessi_fn  = "/media/xavi/Saigon/VariationAg1k/data/phase2.AR1/accessibility/accessibility.h5"
haploty_fn  = "/media/xavi/Saigon/VariationAg1k/data/phase2.AR1/haplotypes/zarr2/ag1000g.phase2.ar1.samples/"
snpeff_fn   = "/media/xavi/Saigon/VariationAg1k/data/phase2.AR1/snpeff/zarr2/"
gffann_fn   = "metadata/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3"

# define populations
popl    = ["BFcol","BFgam","CIcol","GHcol","GHgam","GNgam"]
popc    = "population"
sub1l   = popl
sub1c   = popc
chrom   = "2R"

# Libraries:
import os
import numpy as np
import zarr
import pandas as pd
import allel
import itertools
import matplotlib.pyplot as plt


# ## Load data
# ### Genotypes, haplotypes, variants & samples
# Load data for all variants & genotypes. Population and sample structure:

# load samples list with sample code, groupings, locations etc.
samples_df   = pd.read_csv(metasam_fn, sep='\t')
samples_bool = (
    samples_df[popc].isin(popl).values & 
    samples_df[sub1c].isin(sub1l).values
)
samples_sub  = samples_df[samples_bool]
samples_sub.reset_index(drop=True, inplace=True)

# indexed dictionary of populations
popdict = dict()
for popi in popl: 
    popdict[popi]  = samples_sub[samples_sub[popc] == popi].index.tolist()

# add an extra population composed of all other locations
popdict["all"] = []
for popi in popl:
    popdict["all"] = popdict["all"] + popdict[popi]
    
# report
print("Data:")
print("* Samples     = ", samples_sub.shape[0])
print("* Populations = ", set(samples_sub[popc]))
print(samples_sub.groupby(("population")).size())


# Load variant, genotypes, haplotypes, accessibility, etc. data:
# Accessibility
import h5py
print("Load accessibility array...")
accessi_df  = h5py.File(accessi_fn,mode="r")
accessi_arr = accessi_df[chrom]["is_accessible"][:]


# ## Divergence
# Calculate divergence of genotype frequencies between species and duplicated sequences:

has_dup = np.array([any(b in s for b in ["TRUE"]) for s in samples_sub["pop_dup"].values])

# dictionary with species (nondup) and dups
popdict_div = dict()
popdict_div["dup"] = np.where(has_dup)[0].tolist()
popdict_div["col"] = samples_sub[  np.logical_and( samples_sub["m_s"] == "M" , np.logical_not(has_dup) )  ].index.tolist()
popdict_div["gam"] = samples_sub[  np.logical_and( samples_sub["m_s"] == "S" , np.logical_not(has_dup) )  ].index.tolist()
popdict_div["cod"] = samples_sub[  np.logical_and( samples_sub["m_s"] == "M" , has_dup )  ].index.tolist()
popdict_div["gad"] = samples_sub[  np.logical_and( samples_sub["m_s"] == "S" , has_dup )  ].index.tolist()
popdict_div["all"] = popdict_div["dup"] + popdict_div["gam"] + popdict_div["col"]

# Region of interest:
ace_start = 3489213 # start gene
ace_end   = 3493788 # end gene
ace_119S  = 3492074 # resistant variant
ace_dups  = 3436800 # start duplication
ace_dupe  = 3639600 # end duplication


# Load data from region of interest:
export_start = ace_dups
export_end   = ace_dupe
export_name  = "duplication"

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# calculate population allele counts
hapalco_sub = haploty_sub.count_alleles_subpops(subpops=popdict_div)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)


# Pairwise differences:
# pairwise difference
mpd_dg = allel.mean_pairwise_difference_between(ac1=hapalco_seg["dup"], ac2=hapalco_seg["gam"])
mpd_dc = allel.mean_pairwise_difference_between(ac1=hapalco_seg["dup"], ac2=hapalco_seg["col"])
# jack-knife per-SNP estimates
mpd_dg_est = allel.stats.misc.jackknife(mpd_dg , np.mean)
mpd_dc_est = allel.stats.misc.jackknife(mpd_dc , np.mean)
# report
print("dif gam to dup = %.7f +/- %.7f " % (mpd_dg_est[0], mpd_dg_est[1]))
print("dif col to dup = %.7f +/- %.7f " % (mpd_dc_est[0], mpd_dc_est[1]))


# Dxy divergence:
# dxy per window
dxy_dg = allel.windowed_divergence(ac1=hapalco_seg["dup"], ac2=hapalco_seg["gam"], pos = hapvars_seg["POS"], size=100, is_accessible = accessi_arr)
dxy_dc = allel.windowed_divergence(ac1=hapalco_seg["dup"], ac2=hapalco_seg["col"], pos = hapvars_seg["POS"], size=100, is_accessible = accessi_arr)

# jack-knife per-SNP estimate
dxy_dg_est = allel.stats.misc.jackknife(dxy_dg[0] , np.nanmean)
dxy_dc_est = allel.stats.misc.jackknife(dxy_dc[0] , np.nanmean)

# report
print("dxy gam to dup = %.7f +/- %.7f " % (dxy_dg_est[0], dxy_dg_est[1]))
print("dxy col to dup = %.7f +/- %.7f " % (dxy_dc_est[0], dxy_dc_est[1]))


# PBS relative to duplicated sequences:
pbs_col = allel.pbs(ac1=hapalco_seg["col"], ac2=hapalco_seg["gam"], ac3=hapalco_seg["dup"], window_size=100)
pbs_gam = allel.pbs(ac1=hapalco_seg["gam"], ac2=hapalco_seg["col"], ac3=hapalco_seg["dup"], window_size=100)

# jack-knifing
pbs_col_est = allel.stats.misc.jackknife(pbs_col , np.nanmean)
pbs_gam_est = allel.stats.misc.jackknife(pbs_gam , np.nanmean)

# report
print("PBS col = %.7f +/- %.7f" % (pbs_col_est[0], pbs_col_est[1]))
print("PBS gam = %.7f +/- %.7f" % (pbs_gam_est[0], pbs_gam_est[1]))


# Same, using only duplicated sequences from *col*:
pbs_col = allel.pbs(ac1=hapalco_seg["col"], ac2=hapalco_seg["gam"], ac3=hapalco_seg["cod"], window_size=100)
pbs_gam = allel.pbs(ac1=hapalco_seg["gam"], ac2=hapalco_seg["col"], ac3=hapalco_seg["cod"], window_size=100)

# jack-knifing
pbs_col_est = allel.stats.misc.jackknife(pbs_col , np.nanmean)
pbs_gam_est = allel.stats.misc.jackknife(pbs_gam , np.nanmean)

# report
print("PBS col = %.7f +/- %.7f" % (pbs_col_est[0], pbs_col_est[1]))
print("PBS gam = %.7f +/- %.7f" % (pbs_gam_est[0], pbs_gam_est[1]))


# Same, using only duplicated sequences from *gam*:
pbs_col = allel.pbs(ac1=hapalco_seg["col"], ac2=hapalco_seg["gam"], ac3=hapalco_seg["gad"], window_size=100)
pbs_gam = allel.pbs(ac1=hapalco_seg["gam"], ac2=hapalco_seg["col"], ac3=hapalco_seg["gad"], window_size=100)

# jack-knifing
pbs_col_est = allel.stats.misc.jackknife(pbs_col , np.nanmean)
pbs_gam_est = allel.stats.misc.jackknife(pbs_gam , np.nanmean)

# report
print("PBS col = %.7f +/- %.7f" % (pbs_col_est[0], pbs_col_est[1]))
print("PBS gam = %.7f +/- %.7f" % (pbs_gam_est[0], pbs_gam_est[1]))


# Plot, with a bit of a rolling window and including flanking regions
# reload data to include flanks
export_start = ace_dups-1e6
export_end   = ace_dupe+1e6
export_name  = "duplication"

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# calculate population allele counts
hapalco_sub = haploty_sub.count_alleles_subpops(subpops=popdict_div)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)



# window bp
winsize=5000
stepsize=500
windows_pos = allel.moving_statistic(
    hapvars_seg["POS"][:],
    statistic=lambda v: v[0],
    size=winsize,step=stepsize)

# based on all duplications
pbs_col = allel.pbs(ac1=hapalco_seg["col"], ac2=hapalco_seg["gam"], ac3=hapalco_seg["dup"], window_size=winsize, window_step=stepsize)
pbs_gam = allel.pbs(ac1=hapalco_seg["gam"], ac2=hapalco_seg["col"], ac3=hapalco_seg["dup"], window_size=winsize, window_step=stepsize)

plt.figure(figsize=(5,1.5))
plt.step(windows_pos/1e6, pbs_gam, color="orangered", label="gam to dup", where="post")
plt.step(windows_pos/1e6, pbs_col, color="blue", label="col to dup", where="post")
plt.title("Distance from duplication genotypes")
plt.ylim([-0.01,0.15])
plt.ylabel("Distance")
plt.xlabel("Mb")
plt.axhline(0, color='k',linestyle="--",label="")
plt.axvline(ace_dups/1e6, color='r',linestyle=":",label="")
plt.axvline(ace_dupe/1e6, color='r',linestyle=":",label="")
plt.axvline(ace_start/1e6, color='orange',linestyle=":",label="")
plt.axvline(ace_end/1e6, color='orange',linestyle=":",label="")
plt.legend()
plt.savefig("%s/dist_3pops_alldup.pdf" % (outdir), bbox_inches='tight')
plt.close()



# based on col duplications
pbs_col = allel.pbs(ac1=hapalco_seg["col"], ac2=hapalco_seg["gam"], ac3=hapalco_seg["cod"], window_size=winsize, window_step=stepsize)
pbs_gam = allel.pbs(ac1=hapalco_seg["gam"], ac2=hapalco_seg["col"], ac3=hapalco_seg["cod"], window_size=winsize, window_step=stepsize)

plt.figure(figsize=(5,1.5))
plt.step(windows_pos/1e6, pbs_gam, color="orangered", label="gam to dup", where="post")
plt.step(windows_pos/1e6, pbs_col, color="blue", label="col to dup", where="post")
plt.title("Distance from coluzzii duplication genotypes")
plt.ylim([-0.01,0.15])
plt.ylabel("Distance")
plt.xlabel("Mb")
plt.axhline(0, color='k',linestyle="--",label="")
plt.axvline(ace_dups/1e6, color='r',linestyle=":",label="")
plt.axvline(ace_dupe/1e6, color='r',linestyle=":",label="")
plt.axvline(ace_start/1e6, color='orange',linestyle=":",label="")
plt.axvline(ace_end/1e6, color='orange',linestyle=":",label="")
plt.legend()
plt.savefig("%s/dist_3pops_coldup.pdf" % (outdir), bbox_inches='tight')
plt.close()





# based on gam duplications
pbs_col = allel.pbs(ac1=hapalco_seg["col"], ac2=hapalco_seg["gam"], ac3=hapalco_seg["gad"], window_size=winsize, window_step=stepsize)
pbs_gam = allel.pbs(ac1=hapalco_seg["gam"], ac2=hapalco_seg["col"], ac3=hapalco_seg["gad"], window_size=winsize, window_step=stepsize)

plt.figure(figsize=(5,1.5))
plt.step(windows_pos/1e6, pbs_gam, color="orangered", label="dup to gam", where="post")
plt.step(windows_pos/1e6, pbs_col, color="blue", label="dup to col", where="post")
plt.title("distance from gambiae duplication genotypes")
plt.ylim([-0.01,0.15])
plt.ylabel("Distance")
plt.xlabel("Mb")
plt.axhline(0, color='k',linestyle="--",label="")
plt.axvline(ace_dups/1e6, color='r',linestyle=":",label="")
plt.axvline(ace_dupe/1e6, color='r',linestyle=":",label="")
plt.axvline(ace_start/1e6, color='orange',linestyle=":",label="")
plt.axvline(ace_end/1e6, color='orange',linestyle=":",label="")
plt.legend()
plt.savefig("%s/dist_3pops_gamdup.pdf" % (outdir), bbox_inches='tight')
plt.close()


# ## Export alignment
# Export `Phylip` alignment of haplotypes.
# First, the entire duplication:
export_start = ace_dups
export_end   = ace_dupe
export_name  = "duplication"
poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)


# now, do the same thing but subsetting by haplotype score
hsd_fn  = "results_tables/dupcoverage_stats.HaplotypeScore.csv"
hsd = pd.read_csv(hsd_fn, sep="\t")
hsd_fil = hsd[hsd["HaplotypeScore"]<2.5]["pos"]
export_name  = "duplication_HSdiploid"
poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hap_bool    = np.logical_and(hap_bool , np.isin(hapvars["POS"], hsd_fil) )
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)






# Now, a region **just downstream** of the duplication breakpoint:
export_start = ace_dupe
export_end   = ace_dupe+5e4
export_name  = "breakdodu"
poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)





export_start = ace_dupe
export_end   = ace_dupe+5e3
export_name  = "breakdodu5k"
poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)





export_start = ace_dupe
export_end   = ace_dupe+1e3
export_name  = "breakdodu1k"





poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)


# Now, a region **just upstream** of the duplication breakpoint:




export_start = ace_dups-5e4
export_end   = ace_dups
export_name  = "breakupdu"





poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)





export_start = ace_dups-5e3
export_end   = ace_dups
export_name  = "breakupdu5k"





poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)





export_start = ace_dups-1e3
export_end   = ace_dups
export_name  = "breakupdu1k"





poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)


# Now, an unadmixed region upstream of the duplication:




export_start = ace_dups - 1e6
export_end   = export_start + 5e4
export_name  = "upstream"





poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)


# Now, a region downstream of the duplication:




export_start = ace_dupe + 1e6
export_end   = export_start + 5e4
export_name  = "downstream"





poplhap = popl

# haplotypes: variants
hapcall     = zarr.open(haploty_fn)
print("Load haplotype variants...")
hapcall_var = hapcall[chrom]["variants"]
hapvars     = allel.VariantChunkedTable(hapcall_var,names=["POS","REF","ALT"],index="POS")
hap_bool    = np.logical_and(hapvars["POS"][:] >= export_start, hapvars["POS"][:] <= export_end)
hapvars_sub = hapvars.compress(hap_bool)

# haplotypes: phased genotypes
print("Load haplotype haplotypes...")
hapcall_gen = hapcall[chrom]["calldata/genotype"]
haploty_gen = allel.GenotypeChunkedArray(hapcall_gen)
# find samples in haplotype dataset that coincide with genotypes
haploty_sam = hapcall[chrom]["samples"][:].astype(str)
hapsam_bool = np.isin(haploty_sam, np.array(samples_sub["ox_code"]))
haploty_sub = haploty_gen.subset(sel0=hap_bool,sel1=hapsam_bool)

# recast haplotypes: drop ploidy
print("Drop ploidy haplotypes...")
haploty_sub_hap = haploty_sub.to_haplotypes()

# haplotype dicts
# arrays of hap ids and populations of each hap (double the size of genotype arryays: 2 haps per individual except in X chromosome)
print("Samples dictionary for haps...")
is_samp_in_hap = np.isin(np.array(samples_sub["ox_code"]),haploty_sam)
hap_ids        = np.array(list(itertools.chain(*[[s + 'a', s + 'b'] for s in haploty_sam[hapsam_bool]])))
hap_pops       = np.array(list(itertools.chain(*[[s, s] for s in np.array(samples_sub[popc][is_samp_in_hap])])))
hap_pops_df    = pd.DataFrame(data={ popc : hap_pops , "ids" : hap_ids})

# pop dicts for haplotype data
popdicthap = dict()
for popi in poplhap: 
    popdicthap[popi]  = hap_pops_df[hap_pops_df[popc] == popi].index.tolist()

popdicthap["all"] = []
for popi in poplhap:
    popdicthap["all"] = popdicthap["all"] + popdicthap[popi]

# haplotypes: allele counts
print("Allele counts haplotypes...")
hapalco_sub = haploty_sub_hap.count_alleles_subpops(subpops=popdicthap)

# filter haplotypes: segregating alleles, no singletons
is_hapseg   = hapalco_sub["all"].is_segregating()[:] # segregating
is_hapnosing= hapalco_sub["all"][:,:2].min(axis=1)>2 # no singletons
filhap_bool = (is_hapseg[:] & is_hapnosing[:])

# subset
print("Subset haps...")
haploty_seg = haploty_sub_hap.compress(filhap_bool)
hapvars_seg = hapvars_sub.compress(filhap_bool)
hapalco_seg = hapalco_sub.compress(filhap_bool)

# refs and alts
hapvars_seg_REF = hapvars_seg["REF"][:].astype(str)
hapvars_seg_ALT = hapvars_seg["ALT"][:].astype(str)


# output
print("FASTA...")
happhy = pd.DataFrame({
    "hap": ">"+hap_pops_df["ids"]+"_"+hap_pops_df["population"],
    "seq": np.nan},    
    columns=["hap", "seq"])

for pn,popi in enumerate(hap_pops_df["ids"]):
    
    popi_gen = np.ndarray.tolist(haploty_seg[:,pn])
    popi_seq = [hapvars_seg_REF[gn] if gei == 0 else hapvars_seg_ALT[gn] for gn,gei in enumerate(popi_gen)]
    happhy["seq"][pn] = ''.join(str(e) for e in popi_seq)

happhy.to_csv("%s/hapalignment_%s.fasta" % (outdir,export_name),sep="\n",index=False, header=False)

