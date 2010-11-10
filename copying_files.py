#! usr/bin/python


import os
import shutil

destination = '/home/jaggu/documents/orthomcl/complaint_fasta_files'

try:
    os.mkdir(destination)
except OSError:
    pass
    

for root, dirs, files in os.walk(os.getcwd() + '/parsed_gbk_files_genomes'):
    for file_name in files:
	if file_name.endswith('.fasta'):
	    source = root + '/' + file_name
	    print source, destination
	    shutil.copy(source,destination)

print "Programme : copying_files.py completed!"	
