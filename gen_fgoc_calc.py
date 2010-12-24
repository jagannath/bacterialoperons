#! usr/bin/python

# Script calculates the gen-fgoc score for every organisms. The accession number of the organisms is passed and it calculates the gen-fgoc score. The gen-fgoc score is the 2 * number of conserved adjacent gene pairs (correct direction) / # gene_pairs (ecoli) + # gene_pairs (orgX). This can be normalized with the score obtained between E.coli and E.coli. 

from __future__ import division
import os
import sys



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
    
def create_dic_group_pair_acc_nbr(lines):
    """
    @param lines: This is the set of lines from the walk_in_ecoli_operons_ver2.txt file. 
    @function: Makes a dictionary key pair between group_pair : {[list of accession_number of all conserved adjacent genes] | {[list of accession_number of all non adjacent genes]}. 
    @ return: key: value pair
    """
    # Initializing
    pair_acc_nbr_dictionary = {}
    list_acc_adj = list_acc_non_adj = acc_nbrs = pair = adj_pairs = non_adj_pairs = ''
    
    # (1) Read through lines and obtain - group_pair, details of adjacent gene pairs, details of non adjacent gene pairs
    for line in lines:
	if not line.startswith('///'):
	    pair, adj_pairs, non_adj_pairs = line.split('\t')[1], line.split('\t')[8], line.split('\t')[9]
	    
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
    
    #(1) Iterate over all group pairs and obtain a group pair
    for pair in all_group_pairs:
	# (2) Obtain the value when this group pair is passed as the key
	set_adj_acc_nbrs = pair_acc_nbr_dictionary[pair]
	# (3) Check if the query_acc_nbr is in the set of adj_acc_nbrs
	if org_nbr in set_adj_acc_nbrs:
	    nbr_adjacent_pairs += 1
	
    #returns the nbr_adjacent_pairs across all the pairs for that organism
    return nbr_adjacent_pairs

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
    [org_nbr] = argument
    
    # (1) Opens file walk_in_ecoli_operons_ver2.tmp
    file_name = 'walk_in_ecoli_operons_ver2.txt'
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    
    # (2) Obtain a list of all group_pairs for E.coli (all these are within operon)
    all_group_pairs = [line.split('\t')[1] for line in lines if not line.startswith('///')]

    # (3) Make dictionary pair between pair and all adjacent conserved acc_nbr
    pair_acc_nbr_dictionary = create_dic_group_pair_acc_nbr(lines)

    # (4) Get number of conserved adjacent gene pairs
    nbr_conserved_adjacent_pairs = get_nbr_adjacent_pairs(all_group_pairs,org_nbr, pair_acc_nbr_dictionary)
    
    # (5) Make dictionary pair between accession_number and the number of genes. Made a file from the database using table organisms : accession_nbr \t nbr_genes. Filename is orgs_acc_nbrgenes.txt
    acc_nbrgenes_dictionary = create_dic_acc_nbrgenes('orgs_acc_nbrgenes.txt')
    
    # (6) Obtain number of genes for E.coli and the organism passed as argument
    nbr_ecoli_genes = acc_nbrgenes_dictionary['NC_000913']
    nbr_query_org_genes = (acc_nbrgenes_dictionary[org_nbr])[:-1]
    
    # (7) Calculate the gen_fgoc_score
    genome_fgoc = (2 * int(nbr_conserved_adjacent_pairs)) / (int(nbr_ecoli_genes) + int(nbr_query_org_genes))
    
    # (8) Calculate normalized gen_fgoc_score
    # Normalizing with the gen_fgoc_score obtained when querying with E.coli accession_number. The value neednt be calculated everytime. 
    # ref_fgoc_score = 0.2129
    ref_fgoc = 0.212879314339
    norm_genome_fgoc = (genome_fgoc / ref_fgoc) * 100

    # (9) Write the details to file - orgs.genome_fgoc. It is in appending mode and other genome_fgoc can be added when running sequence
    info_to_write = '\t'.join(item for item in [str(org_nbr), str(nbr_conserved_adjacent_pairs), str(nbr_query_org_genes), str(genome_fgoc), str(norm_genome_fgoc), '\n'])
    ofile = open_file('orgs_ref_ecoli.genome_fgoc','a')
    ofile.write(info_to_write)
    ofile.close()

    return True


if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - gen_fgoc_calc.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    