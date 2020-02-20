# import libraries
import allel
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# input files
chromlist  = ["2R","2L","3R","3L","X"]
kmerbat_fn = "full_sig_kmers.csv"
kmerdat_fn = "windowed_kmer_counts_table.csv"
chrsize_fn = "../metadata/Anogam_gDNA.fasta.fai"
results_fo = "."

# general settings
sns.set(style="ticks",
        font_scale=1.3,
		  rc={"lines.linewidth": 1},
        font="Arial")

## LOAD METADATA
# kmer mapping coordinates
kmerbat = pd.read_csv(kmerbat_fn, sep='\t')
kmerdat = pd.read_csv(kmerdat_fn, sep='\t')


# chr sizes (fai index)
chrsize = pd.read_csv(chrsize_fn, sep="\t", header=None)
chrsize.columns = ["chrom","length","offset","linebases","linewidth"]

size_block = int(5e5)
step_block = int(1e5)

plt.figure(figsize=(24,2))
plt.subplots_adjust(wspace=0.3,hspace=0.6)
for i,chrom in enumerate(chromlist):

	chrlen = chrsize[chrsize["chrom"] == chrom]["length"].values[0]
	chrpos = np.arange(chrlen).astype(int)
	kmerbai = kmerbat[kmerbat["Chrom"] == chrom]

	chrstarts = allel.moving_statistic(chrpos, np.min, size=size_block, step=step_block)
	chrstops  = allel.moving_statistic(chrpos, np.max, size=size_block, step=step_block)
	mapfreqs  = np.zeros(chrstarts.shape)

	for p,start in enumerate(chrstarts):
		count = len(
			kmerbai[ np.logical_and( 
				kmerbai["Pos"] >= chrstarts[p] , 
				kmerbai["Pos"] <= chrstops[p] ) ] )
		mapfreqs[p] = count

	# windows
	kmerdai = kmerdat[kmerdat["Chrom"] == chrom]
	kmerdai_ix = np.where(kmerdai["Kmer.num"].values>0)[0]

	# plot
	ax=plt.subplot(1,5,i+1)
	sns.despine(ax=ax,offset=5)
	plt.step(chrstarts/1e6, mapfreqs, color="blue", where="pre")
	# plt.plot(kmerdai["Window.midpoint"].values/1e6, kmerdai["Kmer.num"].values, color="blue")
	plt.plot(kmerdai["Window.midpoint"].values[kmerdai_ix]/1e6, kmerdai["Kmer.num"].values[kmerdai_ix] + 30, color="red", marker="v", linewidth=0, markersize=5)
	plt.title("%s" % chrom)
	plt.xlabel("Mb")
	plt.ylabel("Frequency")
	plt.xlim(0,65)
	plt.ylim(0,500)
	
# save
plt.savefig("%s/windowed_kmer_counts_table.pdf" % (results_fo), bbox_inches='tight')
plt.close()




