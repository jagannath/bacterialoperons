#! usr/bin/python
# This script reads all the blast results generated in the folder - Blast results, checks if eval is < E-5.  All these are good proteins and the blast result is stored in a file goodproteins.bst. The other proteins are stored in poorproteins.bst. I didnt understand the orthomcl concept of percentage stop codons. The orthomcl in the next steps take care of the evalue which they set for E-5. Thus this file just appends the results of all the blast results file into a file - good proteins

import os
import sys
import sqlite3
from timeit import Timer
from time import localtime, strftime


def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	sys.exit(0)
    
    return file



def append_blast_file(file_name):
    """ This function opens the result file, reads all the lines and returns back the lines. This will be appended to the longer list. I may may a function parse_blast_file if there is a need to parse the blast results based on the Evalue (but for now i am not doing it) 
    """
    
    #Initializing
    path = os.getcwd() + '/blast_200/blast_results'
    file_name_path = path + '/' + file_name
    
    ifile = open_file(file_name_path)
    ofile = open_file('goodproteins.blast','a')
    
    lines = ifile.readlines()
    for line in lines:
	ofile.write(line)
    
    ifile.close()
    ofile.close()
    
    return True
    
 

def open_blast_files(blast_dir = os.getcwd() + '/blast_200/blast_results'):
    """ This function opens the blast result directory that contains 1000s of files, parses each one of them for the criteria of sequence length and evalue (which personally i think orthomcl will handle). 
    """
    
    #Initializing
    all_good_proteins = []
    
    
    
    try:
	file_names = [files for files in os.listdir(blast_dir) if files.endswith('.blast')]
    except:
	raise AssertionError ('Directory %s doesnt exist : %s'%(blast_dir, err))
    
    for file_name in file_names:
	directory_name, blast_file = os.path.split(file_name)
	if append_blast_file(file_name): print "Processing %s .. "%(blast_file)
	

   
    return
    

if __name__ == '__main__':
    

    if os.path.isfile('goodproteins.blast'):
	os.remove('goodproteins.blast')

    open_blast_files()

    t = Timer(stmt = "open_blast_files()",setup = "from __main__ import open_blast_files",)
    print " File Execution Completed : parsing_blast_results.py "
    print t.timeit(0)
    print strftime("%d %b %Y %H:%M:%S",localtime())
    















