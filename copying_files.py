#! usr/bin/python


import os
import shutil


destination = os.getcwd() + '/orthomcl_complaint_fasta'

try:
    os.mkdir(destination)
except OSError:
    pass
    

for root, dirs, files in os.walk(os.getcwd() + '/parsed_gbk_files'):
    for file_name in files:
	if file_name.endswith('.fasta'):
	    source = root + '/' + file_name
	    print source, destination
	    shutil.copy(source,destination)
	
