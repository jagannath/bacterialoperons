#! usr/bin/python

# The script parses the .gbk files of all the genomes including plasmids. It looks for the entire DNA sequence, calculates the GC content. It writes to a database for easy retrieval. It will be another database - genome_sequences.db. It also writes a file accession_number.dnaseq into each of its directory

import os
import sys
import sqlite3
from warnings import filterwarnings
import timeit
import re

class




















def open_list(listfile):
    """ This function opens the listfile which is the file containing the path of all the gbk files. It returns a list (i.e. all the file lines)
    """
    # Initializing
    paths = []
    path = ''
    all_paths = []
    
    ifile = open(listfile,'r')
    paths = ifile.readlines()
    ifile.close()
    
    for path in paths:
	all_paths.append(path[:-1])

    
    return all_paths


def parse_gbks():
    """
    @function:(1)Obtains a list of all .gbk files from the file 'all_filenames.lst'. (2) Creates a class object for each of the path (2) Calls class function for entire DNA sequence (3) Calls class function for GC content (4) Calls function to write the DNA sequence into the file 'accession_number.dnaseq' in the directory parsed_gbk_files_genome. Note that it writes only to those accession numbers with the complete genome (not plasmid)
    @files used: all_filenames.lst
    @return: status = complete
    """
    
    #Initializing
    status = 'incomplete'
    
    listfile = 'all_filenames.lst'
    all_gbk_paths = open_list(listfile)
    
    for path in all_gbk_paths:
	gbk = Gbk_file(path)
	
    
    
    
    

def main(argv):
    """
    @param argv:Obtains the list of all arguments passed. If no argument is passed, it just prints the usage() which is the documentation of the file 
    @function: Creates a class object for each genome path. 
    """

    #Initializing
    
    # Calling database using sqlite3 module
    conn = sqlite3.connect('genome_sequences.db')
    cursor = conn.cursor()
    
    # Calling function parse_gbks to parse all the .gbk files for complete DNA sequence information
    parse_gbks()
    
    
    #Closing database connection
    cursor.close()
    conn.commit()
    conn.close()
    
    
    
    
	




if __name__ == '__main__':
    main(sys.argv[1:])
    









