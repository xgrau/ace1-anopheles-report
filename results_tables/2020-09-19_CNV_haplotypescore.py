# import libraries
import allel
import zarr
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats

# input files
metasam_fn = "../metadata/samples.meta_phenotypes.txt"
callset_fn = "/media/xavi/Saigon/VariationAg1k/data/phase2.AR1/variation/main/zarr2/ag1000g.phase2.ar1.pass/"

## LOAD DATA

chrom="2R"	
zarr_fn = "%s/%s" % (callset_fn, chrom)
callset = zarr.open(zarr_fn)
genvars = allel.VariantChunkedTable( callset["variants"],index="POS" )

# region of interest
ace_start = 3489213 # start gene
ace_end   = 3493788 # end gene
ace_119S  = 3492074 # resistant variant
ace_dups  = 3436800 # start duplication
ace_dupe  = 3639600 # end duplication

# tagging variants
tags = [3465693,3469441,3481632,3504796]


pos_start = min(tags)-5e3
pos_end = max(tags)+5e3
pos_bool = np.logical_and(genvars["POS"] >= pos_start , genvars["POS"] < pos_end)



winsize = 100
winstep = 5


pdf_pages = PdfPages("dupcoverage_stats.pdf")

# read depth
var = "DP"
ar_pos = allel.moving_statistic(genvars["POS"].compress(pos_bool)/1e6, np.median, size=winsize, step=winstep)
ar_loc = allel.moving_statistic(genvars[var].compress(pos_bool), np.median, size=winsize, step=winstep)
plt.figure(figsize=(12,2))
plt.step(ar_pos, ar_loc, color="blue", where="mid")
for tag in tags:
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
plt.axvline(x=ace_119S/1e6, color='orange', linestyle='--')
plt.axvline(x=ace_start/1e6, color='gold', linestyle='--')
plt.axvline(x=ace_end/1e6, color='gold', linestyle='--')
plt.title(var)
plt.xlabel("Mb")
plt.ylabel(var)
# save
pdf_pages.savefig(bbox_inches='tight')
plt.close()

# hap score (2=diploid)
var="HaplotypeScore"
ar_pos = allel.moving_statistic(genvars["POS"].compress(pos_bool)/1e6, np.mean, size=winsize, step=winstep)
ar_loc = allel.moving_statistic(genvars[var].compress(pos_bool), np.median, size=winsize, step=winstep)
plt.figure(figsize=(12,2))
plt.step(ar_pos, ar_loc, color="blue", where="mid")
for tag in tags:
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
plt.axvline(x=ace_119S/1e6, color='orange', linestyle='--')
plt.axvline(x=ace_start/1e6, color='gold', linestyle='--')
plt.axvline(x=ace_end/1e6, color='gold', linestyle='--')
plt.title(var)
plt.xlabel("Mb")
plt.ylabel(var)
plt.ylim((0,12))
plt.axhline(2, linestyle="--", color="black")
# save
pdf_pages.savefig(bbox_inches='tight')
plt.close()

ar_dat = pd.DataFrame(data={
	"pos" : genvars["POS"].compress(pos_bool),
	"HaplotypeScore": genvars[var].compress(pos_bool)
})
ar_dat.to_csv("dupcoverage_stats.HaplotypeScore.csv", sep="\t", index=False)

# fisher test (strand bias) (low = bias)
var="FS"
ar_pos = allel.moving_statistic(genvars["POS"].compress(pos_bool)/1e6, np.mean, size=winsize, step=winstep)
ar_loc = allel.moving_statistic(genvars[var].compress(pos_bool), np.median, size=winsize, step=winstep)
plt.figure(figsize=(12,2))
plt.step(ar_pos, ar_loc, color="blue", where="mid")
for tag in tags:
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
plt.axvline(x=ace_119S/1e6, color='orange', linestyle='--')
plt.axvline(x=ace_start/1e6, color='gold', linestyle='--')
plt.axvline(x=ace_end/1e6, color='gold', linestyle='--')
plt.title(var)
plt.xlabel("Mb")
plt.ylabel(var)
plt.ylabel(var)
# save
pdf_pages.savefig(bbox_inches='tight')
plt.close()

# ReadPosRankSum
var="ReadPosRankSum"
ar_pos = allel.moving_statistic(genvars["POS"].compress(pos_bool)/1e6, np.mean, size=winsize, step=winstep)
ar_loc = allel.moving_statistic(genvars[var].compress(pos_bool), np.median, size=winsize, step=winstep)
plt.figure(figsize=(12,2))
plt.step(ar_pos, ar_loc, color="blue", where="mid")
for tag in tags:
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
plt.axvline(x=ace_119S/1e6, color='orange', linestyle='--')
plt.axvline(x=ace_start/1e6, color='gold', linestyle='--')
plt.axvline(x=ace_end/1e6, color='gold', linestyle='--')
plt.title(var)
plt.xlabel("Mb")
plt.ylabel(var)
plt.ylim((-1,1))
plt.axhline(0, linestyle="--", color="black")
# save
pdf_pages.savefig(bbox_inches='tight')
plt.close()

# qual score
var="QUAL"
ar_pos = allel.moving_statistic(genvars["POS"].compress(pos_bool)/1e6, np.mean, size=winsize, step=winstep)
ar_loc = allel.moving_statistic(genvars[var].compress(pos_bool), np.median, size=winsize, step=winstep)
plt.figure(figsize=(12,2))
plt.step(ar_pos, ar_loc, color="blue", where="mid")
for tag in tags:
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
	plt.axvline(x=(tag)/1e6, color='magenta', linestyle='--')
plt.axvline(x=ace_119S/1e6, color='orange', linestyle='--')
plt.axvline(x=ace_start/1e6, color='gold', linestyle='--')
plt.axvline(x=ace_end/1e6, color='gold', linestyle='--')
plt.title(var)
plt.xlabel("Mb")
plt.ylabel(var)
# save
pdf_pages.savefig(bbox_inches='tight')
plt.close()

# pdf close
pdf_pages.close()


# REPEATS
pos_start = ace_dups-5e4
pos_end = ace_dupe+5e4
pos_bool = np.logical_and(genvars["POS"] >= pos_start , genvars["POS"] < pos_end)

pdf_pages = PdfPages("dup_repeats.pdf")
plt.figure(figsize=(8,2))

var="RepeatMasker"
ar_pos = genvars["POS"].compress(pos_bool)/1e6
ar_loc = genvars[var].compress(pos_bool)
plt.step(ar_pos, ar_loc, color="blue", where="pre", label = "RepeatMasker (Repbase & Dfam)")

# var="RepeatTRF"
# ar_loc = genvars[var].compress(pos_bool)
# plt.step(ar_pos, ar_loc, color="orange", where="pre", label="TRF")

var="RepeatDUST"
ar_loc = genvars[var].compress(pos_bool)
plt.step(ar_pos, ar_loc, color="grey", where="pre", label="Dust")

#other
plt.axvline(x=ace_119S/1e6, color='orange', linestyle='--')
plt.axvline(x=ace_start/1e6, color='gold', linestyle='--', label="Ace1")
plt.axvline(x=ace_end/1e6, color='gold', linestyle='--')
plt.axvline(x=ace_dups/1e6, color='red', linestyle='--', label="duplication")
plt.axvline(x=ace_dupe/1e6, color='red', linestyle='--')
plt.legend(loc='upper left',fontsize='x-small', bbox_to_anchor=(1, 1))
plt.title("Repeats")
plt.xlabel("Mb")
plt.ylabel("Presence")
# save
pdf_pages.savefig(bbox_inches='tight')
plt.close()

# pdf close
pdf_pages.close()
