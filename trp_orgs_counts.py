#! usr/bin/python

# The main question is - In a particular operon design - Trp E->D->C->B->A, how many organisms are in (E->D)..and then how many organisms have (E->D->B) and so on..

import os
import sys



def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	sys.exit(0)

    return file




def get_trp_file():
    """
    @function: Opens the trp_gene_orders.txt file in the directory tryptophan_study and returns the lines
    """
    trp_file = os.getcwd() + '/tryptophan_study/trp_gene_orders.txt'
    ifile = open_file(trp_file)
    lines = ifile.readlines()
    ifile.close()
    
    return lines
    
def get_number_gene_walk(contents):
    """
    @param contents: This is the lines in the trp_gene_orders.txt
    @function: Reads each line and for each line - obtains a set of all organisms conserved in the gene walk. 
    """
    #Initializing
    set_details = []
    walk_number = conserved_details = ''
    line = ''
    all_orgs_in_sets = []
       
    for line in contents:
	walk_number, conserved_details = line.split('\t')[1], line.split('\t')[5]
	
	all_orgs = []
	for gene_pair in eval(conserved_details):
	    gene_1 = gene_pair[0]
	    orgs_acc_no = gene_1[5]
	    all_orgs.append(orgs_acc_no) 
	    
	for orgs in all_orgs:
	    all_orgs_in_sets.append(orgs)
	
	set_details.append([walk_number, all_orgs, len(all_orgs)])
    print len(set(all_orgs_in_sets))
	
	
    return set_details


def get_set_intersection(set1,set2):
    """
    @param set1 and set2: These are set in the form ['walk_number',[set of all_orgs conserved in the walk], number of orgs in the set],..
    @function: Checks the number of times orgs in set 1 occurs in set 2 and what orgs are in both
    @return: ['1_2',[set of all orgs conserved in both sets], number of orgs in this set]
    """
    all_orgs = []
    
    walk_concat = set1[0] + '_' + set2[0]
    set1_orgs = set1[1]
    set2_orgs = set2[1]
    for orgs_1 in set1_orgs:
	for orgs_2 in set2_orgs:
	    if orgs_1 == orgs_2:
		all_orgs.append(orgs_1)

    set_details = [walk_concat, all_orgs, len(set(all_orgs))]
    
    return set_details
    
    

def get_double_walk_orgs(set_details):
    """
    @param set_details: This is a list containing [['walk_number',[set of all_orgs conserved in the walk],number of orgs in the set],..]
    @function: Walks between two such sets - runs from set1 and checks if each orgs in set 1 is present in set2. If present adds this orgs. Obtains the number of orgs conserved between the two sets
    """
    all_set_details = []
    all_orgs_in_sets = []
    set1 = set_details[0]
    i = 1
    
    
    for a_set in set_details:
	if i>1:
	    set2 = a_set
	    set_intersection = get_set_intersection(set1,set2)
	    orgs_in_set = set_intersection[1]
	    for orgs in orgs_in_set:
		all_orgs_in_sets.append(orgs)
	    all_set_details.append(set_intersection)
	    set1 = set2
	i+=1
    
    print len(set(all_orgs_in_sets))
    
    return [all_set_details, all_orgs_in_sets]

#def get_triple_walks_orgs(set_details):
    #"""
    #@param set_details: This is exactly the same as that of the get_double_walk_orgs function. 
    #@function: This computes the intersection between two sets 
    #@return: This returns the set in the same format - 
    #"""
    


if __name__ == '__main__':
    
    # (1) Imports file - trp_gene_orders.txt : This file contains the details of the organisms and the gene pairs that were conserved in the reference (ecoli) tryptophan operon
    file_contents = get_trp_file()
    
    # (2) Get number of organisms in each of the gene walk [nbr in walk1, nbr in walk2 ..]
    set_details_walks = get_number_gene_walk(file_contents)
    
    # (3) Get the set of all organisms conserved in a double gene walk i.e intersection of set (orgs-walk1) and set (orgs-walk2) and so on..
    [set_double_walks,all_orgs_in_sets] = get_double_walk_orgs(set_details_walks)
    
    # (4) Get the set of all organisms conserved in a triple gene walk i.e. intersection of set (orgs walk1-2) and set (orgs walk 2-3) and so on..
    [set_triple_walks,all_orgs_in_sets] = get_double_walk_orgs(set_double_walks)
    
    # (5) Get the set of all organisms conserved in all the four walks, i.e in how many orgs is the entire operon conserved
    [set_four_walks,all_orgs_in_sets] = get_double_walk_orgs(set_triple_walks)
