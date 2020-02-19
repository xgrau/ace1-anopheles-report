import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
import pandas as pd
import numpy as np
import multiprocessing as mp
import os
import re
from subprocess import call

if (len(argv) == 3):
	file_list = argv[1]
	path_to_jellyfish = argv[2].rstrip('/')
	threads = '15'
elif (len(argv) == 4):
	file_list = argv[1]
	path_to_jellyfish = argv[2].rstrip('/')
	threads = argv[3]
else:
	raise Exception("Fail. There should be two or three command line arguments (file_list, path_to_jellyfish (, threads)).")


print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')

print('Script was run using the following arguments:\n')
print('\tfile_list = ' + file_list + '\n')
print('\tpath_to_jellyfish = ' + path_to_jellyfish + '\n')
print('\tthreads = ' + threads + '\n\n')
stdout.flush()

# Write a function that will execute a string in bash, but return an error in Python if the bash process does
# not return 0
def bash_process(command_string, process_name, shell = False):
	print('\tUsing command ' + command_string)
	stdout.flush()
	if shell:
		result = call(command_string, shell = shell, executable = '/bin/bash')
	else:
		result = os.system(command_string)
	if result != 0:
		raise Exception('bash process ' + process_name + ' failed with error code ' + str(result))


# Make sure the following awk scripts exist in the directory. If not, create them
awk_code = """#!/usr/bin/awk -f

BEGIN{f = ARGV[1];
      gsub("\\\\..*", "", f)}
{gsub("A", "1", $0);
gsub("C", "2", $0);
gsub("G", "3", $0);
gsub("T", "4", $0);
sub(".{15}", "& ", $0);
print > f _ substr($0, 1, 2) ".txt"}

"""

if not os.path.exists('recode_and_split.sh'):
	with open('recode_and_split.sh', 'w') as f:
		f.write(awk_code)

awk_code2 = """#!/usr/bin/awk -f

BEGIN{f = ARGV[1];
      gsub("[0-9]+\\\\..*", "", f)}
{print > f _ substr($0, 1, 3) ".txt"}

"""

if not os.path.exists('split_further.sh'):
	with open('split_further.sh', 'w') as f:
		f.write(awk_code2)


# Read the table of sample names. This has been written so that the table of fastq files will work just as
# well as a sample manifest. It should work with any spce-delimited file where the first column contains 
# the sample names. 
fastq_files = pd.read_csv(file_list, header = None, sep = ' ')

sample_names = np.unique(fastq_files[0])

def process_jf_sample(s):
	ss = s + '/' + s
	# Dump the jellyfish counts to a table
	jf_dump_command = path_to_jellyfish + '/jellyfish dump -c ' + ss + '.jf -o ' + ss + '-table.csv > ' + ss + '-dump.log'
	print('\tDumping kmer counts for sample ' + s)
	stdout.flush()
	bash_process(jf_dump_command, 'kmer dump command for sample ' + s, shell = True)
	# Recode these files as numerics and split them based on the leading bp of each k-mer. 
	split_recode_command = 'awk -f recode_and_split.sh ' + ss + '-table.csv' 
	print('\tRecoding and splitting kmer count files for sample ' + s)
	stdout.flush()
	bash_process(split_recode_command, 'split_recode command for sample ' + s)
	# Then further split the files starting with A or C (1 or 2)
	files_to_split_further = [ss + '-table' + x + y + '.txt' for x in ['1', '2'] for y in ['1', '2', '3', '4']]
	for f in [ss + '-table' + x + y + '.txt' for x in ['1', '2'] for y in ['1', '2', '3', '4']]:
		further_split_command = 'awk -f split_further.sh ' + f
		bash_process(further_split_command, 'further splitting for sample ' + s)
		# Delete the file so that we don't end up sorting it later
		bash_process('rm ' + f, 'further split file deletion for sample ' + s)
	# Then sort all of the files
	files_to_sort = [s + '/' + x for x in os.listdir(s) if re.search('table\d+.txt', x)]
	print('\tSorting split kmer count files for sample ' + s)
	stdout.flush()
	for f in files_to_sort:
		f_root = re.sub('\.txt', '', f)
		sort_command = 'sort ' + f + ' > ' + f_root + '-sorted.txt'
		bash_process(sort_command, 'file sorting for sample ' + s)
		bash_process('rm ' + f, 'unsorted file deletion for sample ' + s)
	# Once all of this is done, delete the original file 
	print('\tDeleting original count file for sample ' + s)
	stdout.flush()
	delete_command = 'rm ' + ss + '-table.csv'
	bash_process(delete_command, 'original count file deletion for sample ' + s)

# We now want to run that function in parallel
print('Beginning parallel processing of kmer counts.\n')
pool = mp.Pool(min([int(threads), len(sample_names)]))
process_result = pool.map_async(process_jf_sample, sample_names)
process_result.get()
print('\nAll processes completed.')
stdout.flush()


print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')

