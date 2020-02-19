from Bio import SeqIO, SeqRecord
from collections import OrderedDict
import re

# Load the fasta file of significant kmers
all_kmers = list(SeqIO.parse('full_sig_kmers.fa', 'fasta'))

# Write a function that checks two sequences to see if they overlap. If so, return the merged sequence.
# If not, return None
def merge_sequences(seq1, seq2, min_overlap = 10, verbose = False):
	# If either string is fully subsumed in the other, output the longer one
	if str(seq1) in str(seq2):
		if verbose:
			print('Found seq1 fully contained in seq2.')
		return seq2
	elif str(seq2) in str(seq1):
		if verbose:
			print('Found seq2 fully contained in seq1.')
		return seq1
	elif str(seq1) in str(seq2.reverse_complement()):
		if verbose:
			print('Found seq1 fully contained in reverse complement of seq2.')
		return seq2.reverse_complement()
	elif str(seq2.reverse_complement()) in str(seq1):
		if verbose:
			print('Found reverse complement of seq2 fully contained in seq2.')
		return seq1
	# If one sequence is not subsumed by the other, then the maximum possible overlap is the shorter of
	# the two sequence lengths - 1
	max_overlap = min(len(seq1), len(seq2)) - 1
	for i in range(max_overlap, min_overlap - 1, -1):
		if seq1[-i:] == seq2[:i]:
			if verbose:
				print('Found ' + str(i) + 'bp overlap at end of seq1 and start of seq2.')
			return seq1 + seq2[i:]
		elif seq2[-i:] == seq1[:i]:
			if verbose:
				print('Found ' + str(i) + 'bp overlap at end of seq2 and start of seq1.')
			return seq2 + seq1[i:]
		elif seq1[-i:] == seq2.reverse_complement()[:i]:
			if verbose:
				print('Found ' + str(i) + 'bp overlap at end of seq1 and start of reverse complement of seq2.')
			return seq1 + seq2.reverse_complement()[i:]
		elif seq2.reverse_complement()[-i:] == seq1[:i]:
			if verbose:
				print('Found ' + str(i) + 'bp overlap at end of reverse complement of seq2 and start of seq1.')
			return seq2.reverse_complement() + seq1[i:]
	if verbose:
		print('No overlap found')
	return None

# Write a function that runs through all of the kmers in a list and joins any that it can. 
def merge_all_kmers(kmer_list, verbose = True):
	merged_kmers = []
	merged_kmer_names = []
	while(1):
		if verbose:
			print('\tkmers_remaining_to_merge: ' + str(len(kmer_list)))
		still_to_merge = []
		this_kmer = kmer_list[0].seq
		this_namelist = [kmer_list[0].id]
		for i in range(1, len(kmer_list)):
			merged_kmer = merge_sequences(this_kmer, kmer_list[i].seq)
			if merged_kmer is None:
				still_to_merge += [kmer_list[i]]
			else:
				this_kmer = merged_kmer
				this_namelist += [kmer_list[i].id]
		merged_kmers += [this_kmer]
		merged_kmer_names += [this_namelist]
		if len(still_to_merge) == 0:
			break
		else:
			kmer_list = still_to_merge
	return [SeqRecord.SeqRecord(merged_kmers[i], id = '.'.join(merged_kmer_names[i]), description = '') for i in range(len(merged_kmers))]

kmers = all_kmers.copy()
merge_pass = 1
# Keep looping through all the k-mers until no more can be merged.
while(1):
	print('Merge loop number ' + str(merge_pass))
	merged_kmers = merge_all_kmers(kmers)
	if len(merged_kmers) == len(kmers):
		print('Done!')
		break
	else:
		kmers = merged_kmers.copy()
		merge_pass += 1

# The long names for the sequences give problems for pysam later on (which doesn't seem to be able to cope with read
# ids that are too long), so we shorten them and also save a dictionary that will relate the short names to the full
# names (we still want the full names because this contains the information about all the kmers that make it up). 
name_dict = OrderedDict()
for k in merged_kmers:
	new_id = re.sub('\..*', '', k.id)
	name_dict[new_id] = k.id
	k.id = new_id

# Now output a fasta file of the merged kmers
with open('full_sig_merged_kmers.fa', 'w') as f:
	SeqIO.write(merged_kmers, f, 'fasta')

# Now write the dictionary to file
with open('full_sig_kmer_dictionary.txt', 'w') as f:
	for short_name, long_name in name_dict.items():
		f.write(short_name + '\t' + long_name + '\n')


