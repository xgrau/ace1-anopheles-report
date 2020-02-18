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
kmerdat_fn = "windowed_kmer_counts_table.csv"
results_fo = "."

# general settings
sns.set(style="ticks",
        font_scale=1.3,
		  rc={"lines.linewidth": 1},
        font="Arial")

logging.basicConfig(
	level=logging.DEBUG, 
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s"
	)

## LOAD METADATA
# population metadata
kmerdat = pd.read_csv(kmerdat_fn, sep='\t')

plt.figure(figsize=(24,2))
plt.subplots_adjust(wspace=0.3,hspace=0.6)
for i,chrom in enumerate(chromlist):

	kmerdai = kmerdat[kmerdat["Chrom"] == chrom]
	kmerdai_ix = np.where(kmerdai["Kmer.num"].values>0)[0]
	# plot
	ax=plt.subplot(1,5,i+1)
	sns.despine(ax=ax,offset=5)
	plt.plot(kmerdai["Window.midpoint"].values/1e6, kmerdai["Kmer.num"].values, color="blue")
	plt.plot(kmerdai["Window.midpoint"].values[kmerdai_ix]/1e6, kmerdai["Kmer.num"].values[kmerdai_ix] + 30, color="red", marker="v", linewidth=0, markersize=5)
	plt.title("%s" % chrom)
	plt.xlabel("Mb")
	plt.ylabel("Frequency")
	plt.xlim(0,65)
	plt.ylim(0,500)
	
# save
plt.savefig("%s/windowed_kmer_counts_table.pdf" % (results_fo), bbox_inches='tight')
plt.close()




