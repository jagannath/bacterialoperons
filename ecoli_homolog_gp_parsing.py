#! usr/bin/python

# Script opens the blastall results of ecoli vs all orgs (1246 in total). It generates a list of all homolog grps
# This is to avoid or try alternate approach to ortholog groups. So homologs for every gene in Ecoli will be made and tabulated

import pickle
import os
import sys
from warnings import filterwarnings
import sqlite3
import math

def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened : %s "%(name_file,err.strerror)
	sys.exit(0)

    return file

def operate_table(table_name,contents,operation='insert'):
    """
    @function: Creating Table - The details of the table columns are passed on as the string with contents. The contents are passed as within the triple quote format
    """

    table_contents = contents
    
    if operation == 'create': 
	cursor.execute(""" DROP TABLE IF EXISTS %s """%table_name)	#This works well only when the warning sign is turned off
    else: 
	pass 

    cursor.execute(table_contents)

    return 

def parse_blast_result(file_name):
    """
    @param file_name: This is the name of the blast result file (with the path)
    @function: Parses the blast file and updates the table to contain the locus_tag, group_name, -log_eval
    """
    # Initializing
    orgs_acc_nbr = (file_name[:-13])[9:]
    print "Parsing Blast result against %s"%(orgs_acc_nbr)
    result_dir = '/work/jaggu/blast_results_ecoli_vs_all'
    file_path = result_dir + '/' + file_name 

    
    # (1) Open the file and read lines
    ifile = open_file(file_path)
    lines = ifile.readlines()
    ifile.close()
    
    # (2) Iterate through every line and obtain the 1st, 2nd and 11th column (0,1 and 10 element)
    for line in lines:
	# (3) Obtain the items and perform the necessary conversion
	locus_tag, group_name, log_eval = line.split('\t')[1], line.split('\t')[0], line.split('\t')[10]
	locus_tag = locus_tag.split('|')[1]
	group_name = 'ECHOM_' + group_name.split('|')[1]
	try: 
	    log_eval = -1 * math.log10(float(log_eval))	#Can convert the evalue - in the form xe-y to float (not convertable to int)
	except OverflowError:
	    log_eval = 10000
	
	if log_eval > 3:
	    # (4) Update into table
	    insert_contents = '''INSERT INTO ecoli_homolog_group
				(locus_tag, group_name, log_eval)
				VALUES
				('%s','%s','%s')
			    '''%(locus_tag, group_name, log_eval)
	    cursor.execute(insert_contents)
	
    return

def main(argument):
    
    [argv] = argument
    if argv == '1':
	
	# (1) Obtain a list of all ecoli_vs_*.blast_results
	result_dir = '/work/jaggu/blast_results_ecoli_vs_all'
	result_files = os.listdir(result_dir)
	
	# (2) Make a table ecoli_homolog_group containing locus_tag and group_name (Hom_bnumber) and log (evalue)
	table_name = 'ecoli_homolog_group'
	contents = '''CREATE TABLE ecoli_homolog_group
		    (
		    locus_tag VARCHAR (20),
		    group_name VARCHAR (20),
		    log_eval FLOAT  
		    )
		    '''
	#operate_table(table_name, contents, 'create') #Table already created
	
	# (3) For each Blast result, parse it to obtain three columns - (a) bnumber (ecoli) (b) locus tag and (c) evalue. Write them to the local database trial_all_orgs.db
	for result in result_files:
	    parse_blast_result(result)

    if argv == '2':
	# Initializing
	all_locus_tag = []
	unique_locus_tag = []
	
	# (1) Obtail lines in file - ecoli_homolog_data.txt
	ifile = open_file('ecoli_homolog_wgp.txt')
	lines = ifile.readlines()
	ifile.close()
	
	# (2) Run through the lines and eliminate the redundant locus_tags
	for line in lines:
	    locus_tag = line.split('\t')[1]
	    if not locus_tag in all_locus_tag:
		all_locus_tag.append(locus_tag)
		# Drop the group name as it is causing confusion and the tables joined correctly but only the early entries are taken
		line = line[:-1]
		list_line = [item for item in line.split('\t')]
		
		list_line.append('ECHOM_'+locus_tag)
		seq_line = '\t'.join(item for item in list_line) + '\n'
		unique_locus_tag.append(seq_line)
	
	# (3) Write these into a file
	ofile = open_file('ecoli_homolog_data.txt','w')
	for line in unique_locus_tag:
	    ofile.write(line)
	ofile.close()
	
	
	
    return True


if __name__ == '__main__':
    
    argument = sys.argv[1:]
    # Calling database using sqlite module
    conn = sqlite3.connect('trial_all_orgs.db')
    cursor = conn.cursor()
    
    print "Processing %s ..."%(argument)
    main(argument)    
    
    cursor.close()
    conn.commit()
    conn.close()
    
    
    import time
    print "Script - ecoli_homolog_gp_parsing.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    