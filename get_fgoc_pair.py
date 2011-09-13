#!/usr/bin/python

from __future__ import division
import os
import sys
try:
   import cPickle as pickle
except:
   import pickle
   
# The programme takes two group argument; It runs through class group after unpickling the two homolog group pkl files. In the class group, it obtains the number of pairwise distances. Computes two numbers fgoc and directed fgoc which is the number of homologs pair/number of gene1 homologs; Each programme requires three files - two pkl files, each for a homolog group and a acc_nbr_genes_dict

class Groups:
    """ Contains information about the group_order which is passed. Computes the number of conserved walks, list of all possible gene_pairs
    """
    # Shelves imported : cluster_dictionary,group_all_locus_dictionary,locus_cds_dictionary
    def __init__(self, group_pair,list1,list2):
	"""
	@param group_pair: Two homolog groups
	@param list1, list2: The list of all genes with its relevant information
	"""
	self.group_1, self.group_2 = group_pair[0], group_pair[1]
	self.list1 = list1
	self.list2 = list2


    def get_orientation_status(self,orientation_a, orientation_b):
	"""
	@param orientation_a: Orientation of gene A
	@param orientation_b: Orientation of gene B
	@function: Compares the orientation of the two genes. It can be forward, reverse or opposite only. Also for the opposite case - there are two classification - converge if A (forward) -> B (reverse) and diverge if A(reverse) -> B (forward). Each is essential to calculate the different score and assign the directions in the graph
	"""
	orientation_status = ''
	direction_status = ''
	if orientation_a == 'forward' and orientation_b == 'reverse':
	    orientation_status = 'opposite'
	    direction_status = 'converge'
	if orientation_a == 'reverse' and orientation_b == 'forward':
	    orientation_status = 'opposite'
	    direction_status = 'diverge'
	if orientation_a == 'forward' and orientation_b == 'forward':
	    orientation_status = 'forward'
	    direction_status = 'directed'
	if orientation_a == 'reverse' and orientation_b == 'reverse':
	    orientation_status = 'reverse'
	    direction_status = 'directed'

	return orientation_status, direction_status

    def compute_rank_difference(self,orientation_status, gene_a, gene_b):
	"""
	Input:	orientation_status and each of the gene details.
	Function: Checks the orientation status and calculates the difference based on whether it is forward or reverse
	Output:	The difference in the rank and the distance separating them
	"""
	# Initializing
	rank_difference = 0
	distance = 0

	rank_a,orientation_a,left_end_a,right_end_a, acc_no_a = int(gene_a[1]), gene_a[2], int(gene_a[3]), int(gene_a[4]), gene_a[5]
	rank_b,orientation_b,left_end_b,right_end_b, acc_no_b = int(gene_b[1]), gene_b[2], int(gene_b[3]), int(gene_b[4]), gene_b[5]
	
	if orientation_status == 'forward':
	    rank_difference = rank_b - rank_a
	    distance = left_end_b - right_end_a
	if orientation_status == 'reverse':
	    # This is a very very special case where rank_a is 1 and orientation = reverse. So its neighbour rank_b = number_genes in the organisms.  
	    if rank_a == 1:
		print "Special rank"
		sys.exit(1)
	    else:
		rank_difference = rank_a - rank_b
		distance = left_end_a - right_end_b

	return rank_difference, distance

    def find_adjacent_pairs(self):
	"""
	@function: The main function in the class, which computes all the pairwise differences between the members of the list
	"""
	# Initializing
	group_1 = self.group_1
	group_2 = self.group_2
	list1 = self.list1
	list2 = self.list2
	gene_a_occurrences = len(list1)
	gene_b_occurrences = len(list2)
	difference_details = []
	non_adjacent_gene_pairs = []
	count = 0
	converge_count = 0
	diverge_count = 0
	
	# (1) Iterate through all the pairs of the two groups and match their accession numbers. If they match, compute the rank difference
	for gene_a in list1:
	    rank_a,orientation_a,left_end_a,right_end_a, acc_no_a = int(gene_a[1]), gene_a[2], int(gene_a[3]), int(gene_a[4]), gene_a[5]
	    for gene_b in list2:
		rank_b,orientation_b,left_end_b,right_end_b, acc_no_b = int(gene_b[1]), gene_b[2], int(gene_b[3]), int(gene_b[4]), gene_b[5]
		# (2) Check if the accession number of the genes match
		if acc_no_a == acc_no_b:
		    # (3) Obtain the orientation_status (forward, reverse or opposite) and direction_status (converge, diverge or opposite)
		    orientation_status,direction_status = self.get_orientation_status(orientation_a,orientation_b)
		    if orientation_status in ['forward','reverse']:
			# (4) Compute the rank difference and the distance
			rank_difference,distance = self.compute_rank_difference(orientation_status, gene_a, gene_b)
			# (5) Check if the difference is 1.
			if rank_difference == 1:
			    difference_details.append([gene_a,gene_b, distance,direction_status])
			    count+=1
			else:	#All those gene_a that don't have a genepair in gene_b
			    non_adjacent_gene_pairs.append([gene_a,gene_b, rank_difference, distance,direction_status])
		    # (6) This is the complicated book-keeping of gene pairs in opposite strands. 
		    if orientation_status == 'opposite':
			rank_difference = rank_b - rank_a
			distance = left_end_b - right_end_a
			if rank_difference == 1:
			    if direction_status == 'converge':
				converge_count+=1
				difference_details.append([gene_a,gene_b, distance,direction_status])
			    if direction_status == 'diverge':
				diverge_count+=1
				difference_details.append([gene_a,gene_b, distance,direction_status])
    			else:	#All those gene_a that don't have a genepair in gene_b
			    non_adjacent_gene_pairs.append([gene_a,gene_b, rank_difference, distance,direction_status])

	if not (difference_details) : difference_details = ['No Adjacent gene pair']

	return gene_a_occurrences,gene_b_occurrences, count, difference_details, non_adjacent_gene_pairs, converge_count, diverge_count
	

def pickle_file(file_name, contents,path_dir=os.getcwd()):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    #pkl_dir = os.getcwd() + '/pkl_files/'
    pkl_dir = path_dir + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ofile = open(pkl_file,'wb')
    pickle.dump(contents, ofile)
    ofile.close()
    
    return True
    
def unpickle_file(file_name,path_dir=os.getcwd()):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    pkl_dir = path_dir + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ifile = open(pkl_file)
    contents = pickle.load(ifile)
    ifile.close()
    
    return contents

def main(argument):
    
    query_org_id = 'all'.upper()
    upper_path = '/project/marcotte/jagannath/projectfiles/bacterialoperons/graph_genomes/kegg_other_orgs_for_fgoc/' 
    #upper_path = os.getcwd() + '/thiCEFSGH'
    [pair1,pair2] = argument

    # (1) Obtain the pair and the two pickle file paths
    pair = [pair1,pair2]

    #path = os.getcwd() + '/homologs_ecoli/'
    path = upper_path + query_org_id
    file_name = pair1 + '.pkl'
    list1 = unpickle_file(file_name,path)
    
    #path = os.getcwd() + '/homologs_ecoli/'
    #path = os.getcwd() + '/ecoli_kegg/kegg_ecoli_details/'
    file_name = pair2 + '.pkl'
    list2 = unpickle_file(file_name,path)
    
    # (2) Initializing the class
    gene_pair = Groups(pair,list1,list2)

    
    # (3) Obtain the computed values
    gene_a_occurrences, gene_b_occurrences, number_conserved_gene_pair, difference_details, non_adjacent_gene_pairs, converge_count, diverge_count = gene_pair.find_adjacent_pairs()
    #print "yes"
    #print difference_details
    #sys.exit(1)
    
    # (4) Get all the acc_nbrs linked to the non adjacent genes
    non_adjacent_acc_nbrs = set([(item[0])[5] for item in non_adjacent_gene_pairs])
    
    # (5) Obtain Accession numbers of all genes in difference_details
    if len(difference_details) == 1:
	number_conserved_gene_pair = 0
    else:
	acc_nbrs = set([(item[0])[5] for item in difference_details])
	# The number of conserved_gene_pair is always the number of non - redundant acc_nbrs set;
	# Fix -- 10thJun11 : Number of conserved gene pair represents the directed fgoc count.
	# Count the non_adj_acc_nbrs that was not in the conserved_acc_nbrs to get a count of number of non_adjacent_acc_nbrs
	count = 0
	for non_adj_org in non_adjacent_acc_nbrs:
	    if not non_adj_org in acc_nbrs:
		count+=1
	number_non_adjacent_gene_pairs = count
    
    # (6) Calculate the fgoc score (Remember to import _future_ division)
    number_non_adjacent_gene_pairs = len(non_adjacent_acc_nbrs)
    try:
	fgoc_score = number_conserved_gene_pair / (number_conserved_gene_pair + number_non_adjacent_gene_pairs)
	directed_fgoc_score = number_conserved_gene_pair / gene_a_occurrences
	converge_fgoc_score = converge_count / gene_a_occurrences
	diverge_fgoc_score = diverge_count / gene_a_occurrences
    except ZeroDivisionError:
	fgoc_score = 0
	directed_fgoc_score = 0
	converge_fgoc_score = 0
	diverge_fgoc_score = 0

    print "Processing .. %s"%(pair)
    #print number_conserved_gene_pair, number_non_adjacent_gene_pairs, converge_count, diverge_count, gene_a_occurrences, gene_b_occurrences
    #print fgoc_score, directed_fgoc_score, converge_fgoc_score, diverge_fgoc_score
    
    # (7) Writing the information to the file. 
    #information_for_file = '\t'.join(item for item in [str(pair), str(fgoc_score), str(number_conserved_gene_pair), str(number_non_adjacent_gene_pairs), str(gene_a_occurrences), str(gene_b_occurrences), str(difference_details), str(non_adjacent_gene_pairs),'\n'])
    
    #file_name = 'all_pairs_fgoc_pkl_ko.txt'
    #ofile = open(file_name,'a')
    #ofile.write(information_for_file)
    #ofile.close()
    
    # (8) Make the Output file as well
    if fgoc_score > 0:
	file_name = '/project/marcotte/jagannath/projectfiles/regulatory_operons/TPP_riboswitch/' + pair1 + '_' + pair2 + '.txt.' + query_org_id
	ofile = open(file_name,'w')
	info_line = '\t'.join(item for item in [str(pair[0]), str(pair[1]), str(fgoc_score), str(directed_fgoc_score), str(converge_fgoc_score), str(diverge_fgoc_score), str(number_conserved_gene_pair), str(converge_count), str(diverge_count), str(gene_a_occurrences), str(gene_b_occurrences),'\n'])
	ofile.write(info_line)
	ofile.close()

    return True
	
if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - get_fgoc_pair.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
