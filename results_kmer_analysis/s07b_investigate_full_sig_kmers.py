import numpy as np
import pandas as pd
import pysam

# Write a function to extract the information of interest from a given read alignment
def chrom_and_pos(al):
	if al.is_unmapped:
		return('NA', al.pos, 'NA',  al.seq, al.query_name, 0, False)
	else:
		return(al.reference_name, al.pos, al.cigarstring, al.seq, al.query_name, al.reference_length, al.is_supplementary)

# Load the sam file
alignment = pysam.AlignmentFile('full_sig_merged_kmers.sam')

# Create the output table
summary_table = pd.DataFrame(np.array([chrom_and_pos(x) for x in alignment]), columns = ['Chrom', 'Pos', 'Cigar', 'Seq', 'Kmer', 'Ref_length', 'Supplementary'])
summary_table['Pos'] = summary_table['Pos'].astype('int')

chrom_counts = np.unique(summary_table['Chrom'], return_counts = True)

summary_table['In_Dup1'] = (summary_table['Chrom'] == '2R') & (summary_table['Pos'] > 3436500) & (summary_table['Pos'] < 3640000)
summary_table['In_Ace1'] = (summary_table['Chrom'] == '2R') & (summary_table['Pos'] >= 3484107) & (summary_table['Pos'] <= 3495790)

# Save the table 
summary_table.to_csv('full_sig_kmers.csv', sep = '\t', index = False)


