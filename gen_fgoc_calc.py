#! usr/bin/python

# Script calculates the gen-fgoc score for every organisms. The accession number of the organisms is passed and it calculates the gen-fgoc score. The gen-fgoc score is the 2 * number of conserved adjacent gene pairs (correct direction) / # gene_pairs (ecoli) + # gene_pairs (orgX). This can be normalized with the score obtained between E.coli and E.coli. 

from __future__ import division
import os
import sys
import shelve
import pickle



def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened : %s "%(name_file,err.strerror)
	sys.exit(0)

    return file
    
def get_acc_nbr(pairs):
    """
    @param pairs: Gets the list (makes a set) of all accession_numbers from the gene details of all pairs
    @function: parses the list and gets all the accession_numbers. Then makes a set of non-redundant accession_numbers
    @returns: set of all the accession_numbers
    """
    # Initializing
    list_acc_nbr = []
    acc_nbr = ''
    
    # (1) Change the string list to eval list
    all_pairs = eval(pairs)
    
    # (2) Iterate among the pairs and obtain the 5th element in the 1st pair which is the accession_number
    for pair in all_pairs:
	# Condition to eliminate those cases where there is no adjacent gene pairs. Also if there is no non adjacent gene pairs it is an empty set. 
	if not str(pair).startswith('No'):
	    acc_nbr = (pair[0])[5]
	    
	    #(3) Add to the list of all_accession_numbers
	    list_acc_nbr.append(acc_nbr)
	else:
	    list_acc_nbr = []
    # (4) Make a set of these accession_numbers (it is non redundant) and return this
    return set(list_acc_nbr)
    
def create_dic_group_pair_acc_nbr(lines,status):
    """
    @param lines: This is the set of lines from the walk_in_ecoli_operons_ver2.txt file. 
    @param status: This is to determine whether we have to calculate between operon pairs(2) or within operon pairs (1)
    @function: Makes a dictionary key pair between group_pair : {[list of accession_number of all conserved adjacent genes] | {[list of accession_number of all non adjacent genes]}. 
    @ return: key: value pair
    """
    # Initializing
    pair_acc_nbr_dictionary = {}
    list_acc_adj = list_acc_non_adj = acc_nbrs = pair = adj_pairs = non_adj_pairs = ''
    
    if status == '1':
	# (1) Read through lines and obtain - group_pair, details of adjacent gene pairs, details of non adjacent gene pairs
	for line in lines:
	    if not line.startswith('///'):
		pair, adj_pairs, non_adj_pairs = line.split('\t')[1], line.split('\t')[8], line.split('\t')[9] # For within operon pairs
		
		# (2) Obtain the list of accession_numbers for all adjacent pairs
		list_acc_adj = get_acc_nbr(adj_pairs)
		# (3) Obtain the list of accession_numbers for all non adjacent but conserved gene pairs
		#list_acc_non_adj = get_acc_nbr(non_adj_pairs)
		
		# (4) Convert the pair, list_acc_adj, list_acc_non_adj to strings and join accession_numbers
		#acc_nbrs = str(list_acc_adj) + '|' + str(list_acc_non_adj)
		acc_nbrs = str(list_acc_adj)
		
		# (5) Create dictionary
		pair_acc_nbr_dictionary[pair] = acc_nbrs

    if status == '2':	#Between operon pairs; Very similar but the line split values differ
	# (1) Read through the lines and obtain - group_pair, details of adjacent gene pairs and details of non adjacent gene pairs
	for line in lines:
	    pair, adj_pairs, non_adj_pairs = line.split('\t')[0], line.split('\t')[6], line.split('\t')[7]
	    # (2) Obtain the list of accession_numbers for all adjacent pairs
	    list_acc_adj = get_acc_nbr(adj_pairs)
	    
	    # (3) Obtain the list of accession_numbers for all non adjacent but conserved gene pairs
	    #list_acc_non_adj = get_acc_nbr(non_adj_pairs)
	    
	    # (4) Convert the pair, list_acc_adj, list_acc_non_adj to strings and join accession_numbers
	    #acc_nbrs = str(list_acc_adj) + '|' + str(list_acc_non_adj)
	    acc_nbrs = str(list_acc_adj)
	    
	    # (5) Create dictionary
	    pair_acc_nbr_dictionary[pair] = acc_nbrs

    return pair_acc_nbr_dictionary

def get_nbr_adjacent_pairs(all_group_pairs, org_nbr, pair_acc_nbr_dictionary):
    """
    @param all_group_pairs: This is a list of all group pairs (ref E.coli operonic group pair. But this may change if between walks are considered)
    @param org_nbr: Accession number of the organisms to compare with
    @param pair_acc_nbr_dictionary: Dictionary - key= group_pair (of ecoli): value set(all organism_accession_number)
    @function: Counts the number of times the group pair was conserved in the set of adjacent genes
    """
    
    # Initializing
    nbr_adjacent_pairs = 0
    all_conserved_pairs = []
    
    #(1) Iterate over all group pairs and obtain a group pair
    for pair in all_group_pairs:
	# (2) Obtain the value when this group pair is passed as the key
	set_adj_acc_nbrs = pair_acc_nbr_dictionary[str(pair)]
	# (3) Check if the query_acc_nbr is in the set of adj_acc_nbrs
	if org_nbr in set_adj_acc_nbrs:
	    nbr_adjacent_pairs += 1
	    all_conserved_pairs.append(pair)
    #returns the nbr_adjacent_pairs across all the pairs for that organism
    return nbr_adjacent_pairs, all_conserved_pairs

def create_dic_acc_nbrgenes(file_name):
    """
    @param file_name: File that has accession_nbr \t nbr_genes
    @function: Create dictionary = keys (accession_number): values (number_genes)
    """
    
    # Initializing
    acc_nbrgenes_dictionary = {}
    # (1) Open file and read lines
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    
    for line in lines:
	# (2) Obtain the accession_number and the number of genes
	accession_number, number_genes = line.split('\t')[0], line.split('\t')[1]
	# (3) Create the key: value pair
	acc_nbrgenes_dictionary[accession_number] = number_genes

    return acc_nbrgenes_dictionary


def main(argument):
    """
    @param argument: Passed as the accession number of the organsism 
    @function: Computes the gen-fgoc calculated with respect to E.coli
    """
    
    #Initializing
    line = lines = pair = ''
    all_group_pairs = []
    #assert isinstance(argument,str)	#Checks for correct argument
    [org_nbr,status] = argument[0],argument[1]
    
    if status == '-2':
	
	# (1) Opens file walk_in_ecoli_operons_ver2.tmp
	file_name = 'bw_ecoli_operons_homology_ver6.txt'
	ifile = open_file(file_name)
	lines = ifile.readlines()
	ifile.close()
	
	# (2) Obtain a list of all group_pairs for E.coli (all these are within operon); Applicable only for between operon pairs
	all_group_pairs = [line.split('\t')[0] for line in lines if not line.startswith('///')]

	# (3) Make dictionary pair between pair and all adjacent conserved acc_nbr

	pair_acc_nbr_dictionary = create_dic_group_pair_acc_nbr(lines,'2')	#1 is for status - within (1) and between (2)
	print pair_acc_nbr_dictionary.keys()[0:20]
	
	
	# (4) Shelve the pair: acc_nbr dictionary; This is useful only for between operon pairs
	shelve_file = os.getcwd() + '/shelve_files/bw_operons_pair:acc_nbrs_dictionary'
	s = shelve.open(shelve_file)
	for key, value in pair_acc_nbr_dictionary.items():
	    s[key] = value
	s.close()
	print "Pair : [Accession numbers]  shelved"
	print pair_acc_nbr_dictionary.keys()[0:20]
	# (5) Pickle list of all group_pairs; This is useful only for between operon pairs
	pkl_file = os.getcwd() + '/pkl_files/bw_operons_group_pairs_list'
	ofile = open(pkl_file,'wb')
	pickle.dump(all_group_pairs, ofile)
	ofile.close()
	print "All group pairs (between operons) pickled" 
	
	# (6) Shelve accession number : number of genes dictionary; This is common to both between and within operon pairs
	acc_nbrgenes_dictionary = create_dic_acc_nbrgenes('orgs_acc_nbrgenes.txt')
	
	shelve_file = os.getcwd() + '/shelve_files/acc_nbr:number_genes_dictionary'
	s = shelve.open(shelve_file)
	for key, value in acc_nbrgenes_dictionary.items():
	    s[str(key)] = value
	s.close()
	print acc_nbrgenes_dictionary.keys()[0:10]
	
	print "Accession number : number of genes shelved"

    if status == '-1':
	
	# (1) Opens file walk_in_ecoli_operons_ver2.tmp
	file_name = 'walk_in_ecoli_operons_ver6_homolog.txt'
	ifile = open_file(file_name)
	lines = ifile.readlines()
	ifile.close()
	
	# (2) Obtain a list of all group_pairs for E.coli (all these are within operon); Applicable only for between operon pairs
	all_group_pairs = [line.split('\t')[1] for line in lines if not line.startswith('///')]

	# (3) Make dictionary pair between pair and all adjacent conserved acc_nbr

	pair_acc_nbr_dictionary = create_dic_group_pair_acc_nbr(lines,'1')	#1 is for status - within (1) and between (2)
	print pair_acc_nbr_dictionary.keys()[0:20]
	
	
	# (4) Shelve the pair: acc_nbr dictionary; This is useful only for between operon pairs
	shelve_file = os.getcwd() + '/shelve_files/within_operons_pair:acc_nbrs_dictionary'
	s = shelve.open(shelve_file)
	for key, value in pair_acc_nbr_dictionary.items():
	    s[key] = value
	s.close()
	print "Pair : [Accession numbers]  shelved"
	
	# (5) Pickle list of all group_pairs; This is useful only for between operon pairs
	pkl_file = os.getcwd() + '/pkl_files/within_operons_group_pairs_list'
	ofile = open(pkl_file,'wb')
	pickle.dump(all_group_pairs, ofile)
	ofile.close()
	print "All group pairs (within operons) pickled" 
	
	# (6) Shelve accession number : number of genes dictionary; This is common to both between and within operon pairs
	#acc_nbrgenes_dictionary = create_dic_acc_nbrgenes('orgs_acc_nbrgenes.txt')
	
	#shelve_file = os.getcwd() + '/shelve_files/acc_nbr:number_genes_dictionary'
	#s = shelve.open(shelve_file)
	#for key, value in acc_nbrgenes_dictionary.items():
	    #s[str(key)] = value
	#s.close()
	#print acc_nbrgenes_dictionary.keys()[0:10]
	
	#print "Accession number : number of genes shelved"
    
	
    if status == '1':	# Within operon pair calculation
	
	shelve_file = os.getcwd() + '/shelve_files/acc_nbr:number_genes_dictionary'
	acc_nbrgenes_dictionary = shelve.open(shelve_file)
	
	shelve_file = os.getcwd() + '/shelve_files/within_operons_pair:acc_nbrs_dictionary'
	pair_acc_nbr_dictionary = shelve.open(shelve_file)
	
	pkl_file = os.getcwd() + '/pkl_files/within_operons_group_pairs_list'
	ifile = open(pkl_file)
	all_group_pairs = pickle.load(ifile)
	ifile.close()
	print "All group pairs (within operons) pickled" 
	
	print pair_acc_nbr_dictionary[str(['ECHOM_b0002', 'ECHOM_b0003'])]
	
	## (1) Opens file walk_in_ecoli_operons_ver2.tmp
	#file_name = 'walk_in_ecoli_operons_ver2.txt'
	#ifile = open_file(file_name)
	#lines = ifile.readlines()
	#ifile.close()
	
	## (2) Obtain a list of all group_pairs for E.coli (all these are within operon)
	#all_group_pairs = [line.split('\t')[1] for line in lines if not line.startswith('///')]

	## (3) Make dictionary pair between pair and all adjacent conserved acc_nbr
	#pair_acc_nbr_dictionary = create_dic_group_pair_acc_nbr(lines,status)

	# (4) Get number of conserved adjacent gene pairs
	nbr_conserved_adjacent_pairs, all_conserved_pairs = get_nbr_adjacent_pairs(all_group_pairs,org_nbr, pair_acc_nbr_dictionary)
	
	## (5) Make dictionary pair between accession_number and the number of genes. Made a file from the database using table organisms : accession_nbr \t nbr_genes. Filename is orgs_acc_nbrgenes.txt
	#acc_nbrgenes_dictionary = create_dic_acc_nbrgenes('orgs_acc_nbrgenes.txt')
	
	# (6) Obtain number of genes for E.coli and the organism passed as argument
	nbr_ecoli_genes = acc_nbrgenes_dictionary['NC_000913']
	nbr_query_org_genes = (acc_nbrgenes_dictionary[org_nbr])[:-1]
	
	# (7) Calculate the gen_fgoc_score
	genome_fgoc = (2 * int(nbr_conserved_adjacent_pairs)) / (int(nbr_ecoli_genes) + int(nbr_query_org_genes))
	
	# (8) Calculate normalized gen_fgoc_score
	# Normalizing with the gen_fgoc_score obtained when querying with E.coli accession_number. The value neednt be calculated everytime. 
	# ref_fgoc_score = 0.2129
	ref_fgoc = 0.391012277044
	norm_genome_fgoc = (genome_fgoc / ref_fgoc) * 100
	print genome_fgoc
	print norm_genome_fgoc
	
	# (9) Write the details to file - orgs.genome_fgoc. It is in appending mode and other genome_fgoc can be added when running sequence
	info_to_write = '\t'.join(str(item) for item in [org_nbr, nbr_conserved_adjacent_pairs, nbr_query_org_genes, genome_fgoc, norm_genome_fgoc, all_conserved_pairs, '\n'])
	ofile = open_file('within_orgs_ref_ecoli_ver6.genome_fgoc','a')
	ofile.write(info_to_write)
	ofile.close()

    if status == '2':
	
	# Opening shelves and unpickling
	shelve_file = os.getcwd() + '/shelve_files/bw_operons_pair:acc_nbrs_dictionary'
	pair_acc_nbr_dictionary = shelve.open(shelve_file)
	
	pkl_file = os.getcwd() + '/pkl_files/bw_operons_group_pairs_list'
	ifile = open(pkl_file)
	all_group_pairs = pickle.load(ifile)
	ifile.close()
	
	shelve_file = os.getcwd() + '/shelve_files/acc_nbr:number_genes_dictionary'
	acc_nbrgenes_dictionary = shelve.open(shelve_file)
	
	print pair_acc_nbr_dictionary[str(['ECHOM_b0004', 'ECHOM_b0005'])]

	## (1) Open file - bw_ecoli_operons.ver2.txt; This contains gene pairs between the ecoli operons
	#file_name = 'bw_ecoli_operons_homology_ver6.txt.txt'
	#ifile = open_file(file_name)
	#lines = ifile.readlines()
	#ifile.close()
	
	## (2) Obtain a list of all group_pairs for E.coli (all these are between operon)
	#all_group_pairs = [line.split('\t')[0] for line in lines]

	## (3) Make dictionary pair between pair and all adjacent conserved acc_nbr
	#pair_acc_nbr_dictionary = create_dic_group_pair_acc_nbr(lines,status)

	# (4) Get number of conserved adjacent gene pairs
	nbr_conserved_adjacent_pairs, all_conserved_pairs = get_nbr_adjacent_pairs(all_group_pairs,org_nbr, pair_acc_nbr_dictionary)
	
	## (5) Make dictionary pair between accession_number and the number of genes. Made a file from the database using table organisms : accession_nbr \t nbr_genes. Filename is orgs_acc_nbrgenes.txt
	#acc_nbrgenes_dictionary = create_dic_acc_nbrgenes('orgs_acc_nbrgenes.txt')
	
	# (6) Obtain number of genes for E.coli and the organism passed as argument
	nbr_ecoli_genes = acc_nbrgenes_dictionary['NC_000913']
	nbr_query_org_genes = (acc_nbrgenes_dictionary[org_nbr])[:-1]
	
	# (7) Calculate the gen_fgoc_score
	genome_fgoc = (2 * int(nbr_conserved_adjacent_pairs)) / (int(nbr_ecoli_genes) + int(nbr_query_org_genes))

	# (8) Calculate normalized gen_fgoc_score
	# Normalizing with the gen_fgoc_score obtained when querying with E.coli accession_number. The value neednt be calculated everytime. 
	ref_fgoc = 0.284225156359
	#ref_fgoc = 0.0488765346305
	norm_genome_fgoc = (genome_fgoc / ref_fgoc) * 100

	# (9) Write the details to file - orgs.genome_fgoc. It is in appending mode and other genome_fgoc can be added when running sequence
	info_to_write = '\t'.join(str(item) for item in [org_nbr, nbr_conserved_adjacent_pairs, nbr_query_org_genes,genome_fgoc, norm_genome_fgoc, all_conserved_pairs, '\n'])
	ofile = open_file('bw_operons_orgs_ref_ecoli_ver6.genome_fgoc','a')
	print norm_genome_fgoc
	ofile.write(info_to_write)
	ofile.close()

    return True


if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - gen_fgoc_calc.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    