#! /usr/bin/python

# Script that will compute the operon fgoc score for all the different permutations of the operon. 
from __future__ import division
import sys
import os
import sqlite3
from warnings import filterwarnings
import bsub_fgoc_all_genome as bsub
import itertools
from operator import itemgetter



class Groups:
    """ Contains information about the group_order which is passed. Computes the number of conserved walks, list of all possible gene_pairs
    """
    
    def __init__(self, group_pair, group_dictionary):
	self.group_1, self.group_2 = group_pair[0], group_pair[1]
	self.group_dictionary = group_dictionary
		
    def open_file(self,name_file, open_status = 'r'):
	""" This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """
    
	#Opening and reading/writing the passed file_name """
	try:
	    file = open(name_file,open_status)
	except IOError,err:	# If the file cannot be opened i am exiting
	    raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	    sys.exit(0)
	
	return file
	
    def get_details(self,query_locus_tag, query_acc_no, query_rank):
	""" This function provides the details of the line in question. 
	@param query_locus_tag - The locus_tag of the organism is passed
	@param query_acc_no - The accession_number of the organism is passed
	@param query_rank - The rank of the locus_tag is passed. The latter two can be empty 
	@function - It retrieves the sequence details based on the parameters passed.
	@return - locus_tag, rank, orientation, left_end, right_end, accession_number
	"""
	
	#Initializing

	orientation = left_end = right_end = accession_number = ''
	flag = True
	
	try:
	    assert (query_locus_tag) or (query_acc_no and query_rank)
	except AssertionError:
	    flag = False
	    locus_tag = rank = ''
	    
	if query_acc_no and query_rank:
	    search_command = '''
	    SELECT locus_tag, rank, orientation, left_end, right_end, accession_number
	    FROM coding_sequences
	    WHERE accession_number = '%s' AND rank = %s
	    '''%(query_acc_no, query_rank)
	    
	    cursor.execute(search_command)

	    try: 
		[locus_tag, rank, orientation, left_end, right_end, accession_number] = cursor.fetchone()
	    except TypeError:	#When nothing was returned
		flag = False
		locus_tag = rank = ''
	    

	if query_locus_tag:
	    search_command = '''
	    SELECT locus_tag, rank, orientation, left_end, right_end, accession_number
	    FROM coding_sequences
	    WHERE locus_tag = '%s'
	    '''%(query_locus_tag)
	    cursor.execute(search_command)
	    try: 
		[locus_tag, rank, orientation, left_end, right_end, accession_number] = cursor.fetchone()
	    except TypeError:	#When nothing was returned
		flag = False
		locus_tag = rank = ''
    
	return [locus_tag, rank, orientation, left_end, right_end, accession_number,flag]
			    
    def find_group_members(self, group_name):
	"""
	@param group_name - cluster_accession number  
	@function - This generates the list of all the sequences that belongs to the group. Looks up the table protein_clusters and joins it with the coding_sequences table via the common locus_tag to generate a list of sequences belonging to the group. 
	sequence_details contains = [locus_tag, rank, orientation, left_end, right_end, accession_number,group_name (just for debugging purpose)]
	@param output -group_members - These are the list of all the sequence details in that cluster or group
	"""
	
	#Initializing
	group_members = []
	acc_no = ''
	rank = ''
	
	search_command = '''
	SELECT c.locus_tag, c.rank, c.orientation, c.left_end, c.right_end, c.accession_number, p.cluster_accession
	FROM protein_clusters p
	INNER JOIN coding_sequences c
	ON p.locus_tag = c.locus_tag
	WHERE p.cluster_accession = '%s'
	'''%(group_name)
	
	cursor.execute(search_command)
	group_members = cursor.fetchall()

	return group_members

    def check_orientation(self,orientation_a, orientation_b):
	"""
	Input:	The orientation of the two genes (a and b). They can be 'forward' or 'reverse'. Rank difference is calculated based on that
	Function:Checks the orientation of the two genes. 
	Output:	Returns forward if both forward, reverse if both and opposite if there is a clash
	"""

	if orientation_a == orientation_b:
	    if orientation_a == 'forward':
		orientation_status = 'forward'
	    if orientation_a == 'reverse':
		orientation_status = 'reverse'
	else:
	    orientation_status = 'opposite'
    
	return orientation_status


    def get_special_rank_diff(self,gene_a,gene_b):
	"""
	@param gene_a, gene_b:	information of the two genes
	@ function:		This is the special case calculation of rank difference when the rank_a is 1 and orientation is reverse. I plan to obtain the number of genes and the total size from that particular organism (based on accession number). I will fake the rank difference to be 1 to get this gene_pair into adjacent gene pair which it rightly should. Also the distance will be recalculated based on the genome size. 
	@param rank_difference, distance
	"""
	#Initializing
	rank_a,orientation_a,left_end_a,right_end_a, acc_no_a = int(gene_a[1]), gene_a[2], int(gene_a[3]), int(gene_a[4]), gene_a[5]
	rank_b,orientation_b,left_end_b,right_end_b, acc_no_b = int(gene_b[1]), gene_b[2], int(gene_b[3]), int(gene_b[4]), gene_b[5]
	
	search_command = '''
	SELECT genome_size, number_genes
	FROM organisms
	WHERE accession_number = '%s'
	'''%(acc_no_a)
	cursor.execute(search_command)
	[genome_size, number_genes] = cursor.fetchone()
	
	new_rank_a = rank_a + number_genes
	new_left_end_a = left_end_a + genome_size
	rank_difference = new_rank_a - rank_b
	distance = new_left_end_a - right_end_b

	return rank_difference, distance
	

    def compute_rank_difference(self,orientation_status, gene_a, gene_b):
	"""
	Input:	orientation_status and each of the gene details.
	Function: Checks the orientation status and calculates the difference based on whether it is forward or reverse
	Output:	The difference in the rank and the distance separating them
	"""
	# Initializing
	group_a_members = self.group_a_members
	group_b_members = self.group_b_members
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
		rank_difference, distance = self.get_special_rank_diff(gene_a, gene_b)
	    else:
		rank_difference = rank_a - rank_b
		distance = left_end_a - right_end_b

	return rank_difference, distance


    def get_details_non_adjacent_genes(self):
	"""
	@param Input:	None - Taken from self - list of all non_adjacent_gene_pairs that had a rank difference not equal 1. 
	@function:	1. Gets the first gene - gene_a from each of the member in the list. 
			2. It obviously belongs to the group_a. It checks for the gene next to it in the genome. Gets back details of that gene
			3. Checks if the orientation is the same. If different then - the gene and group is called 'null'
			4. If orientation is correct - it obtains the distance between them and the group which the next gene belonged to.
			5. Sometimes there is no group affliation in which case a keyerror is returned. the group is considered 'null'
			5. Add to list group. Reiterate the steps 1 - 5
	@param out non_adjacent_gene_pair list: this includes [gene_a, gene_b, rank_difference, distance]
	@param out all_non_became_groups: this is a list (non redundant) where the list of all groups (that became when next gene was looked at) 
	"""
	# Initializing
	non_adjacent_gene_pairs = self.non_adjacent_gene_pairs
	locus_tag = rank_a = orientation_a = left_end_a = right_end_a = acc_no_a = ''

	
	next_gene = []
	all_non_became_groups = []
	set_all_non_became_groups = []
	
	try:	# In some cases there are no non_adjacent_gene_pairs. So it is better i send it back without doing any calculation
	    assert non_adjacent_gene_pairs
	except AssertionError:
	    return ['No non_adjacent_gene_pair'],[]
	    
	for gene_pair in non_adjacent_gene_pairs:
	    gene_a = gene_pair[0]
	    rank_a,orientation_a,left_end_a,right_end_a, acc_no_a = gene_a[1], gene_a[2], gene_a[3], gene_a[4], gene_a[5]

	    try:
		assert ((rank_a) and (orientation_a) and (left_end_a) and (right_end_a) and (acc_no_a))
	    
		if orientation_a == 'forward':
		    next_gene = self.get_details(locus_tag, acc_no_a, (rank_a + 1))
		if orientation_a == 'reverse':
		    next_gene = self.get_details(locus_tag, acc_no_a, (rank_a - 1))
		
		flag = next_gene[6]
		# The flag is assigned because sometimes the get_details function cannot get any sequence information from the available coding sequences table. Then in that case i will assign the next group to be zero. 
		
		if flag:
		    locus_tag, rank, orientation, left_end, right_end, accession_number = next_gene[0], next_gene[1], next_gene[2], next_gene[3], next_gene[4], next_gene[5]
		
		    assert (orientation_a) and (orientation)	#Checks whether it retrieved the orientation back
		    
		    orientation_status = self.check_orientation(orientation_a, orientation)
		    
		    if not (orientation_status == 'opposite'):
			rank_difference, distance = self.compute_rank_difference(orientation_status, gene_a, next_gene)

		    if orientation_a == orientation:
			try:
			    next_group = self.group_dictionary[locus_tag]
			except KeyError:
			    next_group = 'Null'
		    else:
			next_group = 'Null'
		else:
		    next_group = 'Null'
	    
	    except AssertionError:
		next_group = 'Null'
		
	    all_non_became_groups.append(next_group)


	return non_adjacent_gene_pairs, all_non_became_groups

    def find_adjacent_pairs(self):
	"""
	@function - This function finds all the adjacent_pairs for the sent group_pair. 
	1. Look into group_1 and obtain the list of all sequence_ids in that group along with necessary details of rank, orientation, position etc
	(Obtains the details from the file - coding_sequences.txt - (sequence_id, rank, orientation, left_end, right_end, sequence_length, accession_nbr, taxon_id))
	2. Do the same with group_2
	3. Do pairwise run through both the groups - check if they have the same accession number, same orientation and obtain (is rank difference 1 (or -1)
	4. If they match that - obtain the distance between the two genes and the gene_pair (via their sequence_id_1:sequence_id_2:distance:taxon_id
	@param output - count, difference_details
	"""
		
	# Initializing
	group_a = self.group_1
	group_b = self.group_2
	count = 0
	difference_details = []
	non_adjacent_gene_pairs = []
	

	self.group_a_members = self.find_group_members(group_a)
	self.group_b_members = self.find_group_members(group_b)
	gene_a_occurrences = len(self.group_a_members)
	gene_b_occurrences = len(self.group_b_members)
	lone_gene_occurrences = abs(gene_a_occurrences - gene_b_occurrences)	#This is the number of genes, which was present in one group for an organism, but for the same organism, the second group didnt have any genes. So the would be group should be considered Null

	for gene_a in self.group_a_members:
	    rank_a,orientation_a,left_end_a,right_end_a, acc_no_a = gene_a[1], gene_a[2], gene_a[3], gene_a[4], gene_a[5]

	    for gene_b in self.group_b_members:
		rank_b,orientation_b,left_end_b,right_end_b, acc_no_b = int(gene_b[1]), gene_b[2], int(gene_b[3]), int(gene_b[4]), gene_b[5]
		
		if acc_no_a == acc_no_b:	#Ensures only the same species is checked
		    orientation_status = self.check_orientation(orientation_a, orientation_b)
		    
		    if not (orientation_status == 'opposite'):
			rank_difference, distance = self.compute_rank_difference(orientation_status, gene_a, gene_b)
			if rank_difference == 1:
			    difference_details.append([gene_a,gene_b, distance])
			    count+=1
			else:	#All those gene_a that don't have a genepair in gene_b
			    non_adjacent_gene_pairs.append([gene_a,gene_b, rank_difference, distance])
	
	self.non_adjacent_gene_pairs = non_adjacent_gene_pairs
	
	if not (difference_details) : difference_details = ['No Adjacent gene pair']

	return gene_a_occurrences,gene_b_occurrences, count, difference_details, non_adjacent_gene_pairs


# Class Groups Ends 


def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened : %s "%(name_file,err.strerror)
	sys.exit(0)

    return file

def compute_within_operon_walks(operon,cluster_dictionary,out_file_name):
    """
    @param operon: Operon name 
    @param cluster_dictionary: The dictionary of locus_tag: cluster_id
    @param out_file_name: The name of the file to which the information computed will be written. 
    @function: One operon name is passed. Computes the FGOC score for each walk within this operon and returns details of all adjacent genes, the genes that weren't adjacent but in the same organism. The total number of gene_pairs conserved etc are calculated. 
    """
    
    # Initializing
    walk_number = 1
    total_groups = []
    total_conserved_gene_pairs = 0
    nbr_group_pairs = 0
    op_fgoc_score = 0
    sum_fgoc_score = 0
    fgoc_score = 0
  
    # (1) Get all genes in the particular operon. genes - information of the gene
    genes, group_order, orientation_status = bsub.get_bsub_genes(operon)
    
    # (2) Check if the number of genes in the operon > 1. This is the correct definition of an operon. If yes continue further. 
    if len(group_order) > 1:
		
	# (3) Reverse group order if the orientation is reverse
	if orientation_status == 'reverse':
	    group_order.reverse()

	# (4) Get group pairs for all the groups present in the operon
	group_pairs = bsub.get_group_pairs(group_order)
	
	# (5) Iterate over each group_pair and obtain the details of total conserved_gene_pairs, number_adjacent_gene_pairs; This will use the class Groups
	for pair in group_pairs:
	    gene_pair = Groups(pair, cluster_dictionary)

	    
	    # Class Groups is called and the details of the find_adjacent_pairs, difference_details etc are obtained. Importantly the became_adjacent gene pair i.e. of the gene_pairs which were not adjacent, the details of what is the gene next to it. It is computationally a little more time consuming but is there in the class Groups and can be modified and used. 
	    gene_a_occurrences, gene_b_occurrences, number_conserved_gene_pair, difference_details, non_adjacent_gene_pairs = gene_pair.find_adjacent_pairs()
	    
	    
	    # (6) Calculate the fgoc score (Remember to import _future_ division)
	    number_non_adjacent_gene_pairs = len(non_adjacent_gene_pairs)
	    fgoc_score = number_conserved_gene_pair / (number_conserved_gene_pair + number_non_adjacent_gene_pairs)
	    
	    # (7) Writing the information to the file. 
	    information_for_file = '\t'.join(item for item in [operon, str(genes), str(walk_number), str(fgoc_score), str(number_conserved_gene_pair), str(number_non_adjacent_gene_pairs), str(gene_a_occurrences), str(gene_b_occurrences), str(difference_details), str(non_adjacent_gene_pairs),'\n'])
	    
	    walk_number += 1

	    #Appending information to file - operon_walk_gene_pairs.txt
	    ofile = open_file(out_file_name,'a')
	    ofile.write(information_for_file)
	    
	    # (8) Summation of fgoc_score
	    sum_fgoc_score += fgoc_score
	
	
	# (9) Calculate the op_fgoc_score
	op_fgoc_score = sum_fgoc_score / len(group_order)
	
	# Keeping a delimiter for every operon.
	ofile.write('Operon FGOC = %s'%(str(op_fgoc_score)))
	ofile.write('\n /// \n')
	ofile.close()
    
    return True

def calculate_op_fgoc(ordering,gene_pair_fgoc_dictionary):
    """
    @param ordering: It is a list of all the genes in a particular order. 
    @function: For each ordering of the operon, the gene pairs are done and the fgoc score for each pair is obtained. This is done using the dictionary generated. Then using the fgoc scores for all gene pairs the op-fgoc_score is computed
    @return: op_fgoc_score
    """
    
    # Initializing
    group_pairs = []
    sum_fgoc = 0
    fgoc = 0
    op_fgoc = 0
    
    # (1) Get group pairs for all the groups present in the operon
    group_pairs = bsub.get_group_pairs(ordering)
    
    # (2) Iterate over all pairs in group_pairs
    for pair in group_pairs:
	joined_pair = pair[0] + '|' + pair[1]
	try:
	    fgoc = gene_pair_fgoc_dictionary[joined_pair]
	except KeyError:
	    print "No key for %s. Exiting..."%(joined_pair)
	    sys.exit(1)
	sum_fgoc += fgoc
    
    op_fgoc = sum_fgoc/len(ordering)
    
    return op_fgoc
    

def compute_shuffled_op_fgoc(operon,cluster_dictionary,out_file_name):
    """
    @param operon: Operon name 
    @param cluster_dictionary: The dictionary of locus_tag: cluster_id
    @param out_file_name: The name of the file to which the information computed will be written. 
    @function: Shuffles the operon. For each operon shuffle a op_fgoc score is calculated. 
    """

    # Initializing
    walk_number = 1
    total_groups = []
    total_conserved_gene_pairs = 0
    nbr_group_pairs = 0
    op_fgoc_score = 0
    sum_fgoc_score = 0
    fgoc_score = 0
    gene_pair_fgoc_dictionary = {}
    sorted_operon_ordering = []
    shuffled_score = []
    
    # (1) Get all genes in the particular operon. genes - information of the gene
    genes, group_order, orientation_status = bsub.get_bsub_genes(operon)
    
    # (2) Check if the number of genes in the operon > 1. This is the correct definition of an operon. If yes continue further. 
    assert len(group_order) > 1

    if len(group_order) < 8:
	# (3) Reverse group order if the orientation is reverse
	if orientation_status == 'reverse':
	    group_order.reverse()

	# (4) Permutation on group_order set to get two gene pairs at a time. [ABCDE] == {AB},{AC},{AD},{AE},{BA}...
	permutation_list = list(itertools.permutations(group_order,2))

	# (5) Iterate over the list to generate group pairs. For each group pair calculate all fgoc_score, gene details etc..
	for pair in permutation_list:

	    gene_pair = Groups(pair,cluster_dictionary)
	    joined_pair = pair[0] + '|' + pair[1]

	    # Class Groups is called and the details of the find_adjacent_pairs, difference_details etc are obtained. Importantly the became_adjacent gene pair i.e. of the gene_pairs which were not adjacent, the details of what is the gene next to it. It is computationally a little more time consuming but is there in the class Groups and can be modified and used. 
	    gene_a_occurrences, gene_b_occurrences, number_conserved_gene_pair, difference_details, non_adjacent_gene_pairs = gene_pair.find_adjacent_pairs()

	    # (6) Calculate the fgoc score (Remember to import _future_ division)
	    number_non_adjacent_gene_pairs = len(non_adjacent_gene_pairs)
	    fgoc_score = number_conserved_gene_pair / (number_conserved_gene_pair + number_non_adjacent_gene_pairs)

	    # (7) Writing the information to the file. 
	    information_for_file = '\t'.join(item for item in [operon, str(pair), str(fgoc_score), str(number_conserved_gene_pair), str(number_non_adjacent_gene_pairs), str(gene_a_occurrences), str(gene_b_occurrences), str(difference_details), str(non_adjacent_gene_pairs),'\n'])

	    # (8) Appending information to file - trp_all_possible_gene_pairs.txt
	    ofile = open_file(os.getcwd() + '/bacillus_shuffled/' + operon + '.shuffled','a')
	    ofile.write(information_for_file)

	    # (9) Make gene_pair : fgoc dictionary
	    gene_pair_fgoc_dictionary[joined_pair] = fgoc_score

	ofile.close()

	# (10) Perform the permutation with the original group_order
	big_permutation_list = list(itertools.permutations(group_order))

	# (11) Iterate through every item on the permutation list. For each item i.e. a particular gene_ordering obtain the op_fgoc_score
	for ordering in big_permutation_list:
	    op_fgoc_score = calculate_op_fgoc(ordering,gene_pair_fgoc_dictionary)
	    shuffled_score.append([ordering, op_fgoc_score])

	# (12) Sort the list - shuffled_score 
	sorted_operon_ordering = sorted(shuffled_score, key = itemgetter(1), reverse=True)


	# (13) Writing the sorted list to the output_file
	for ordering,op_fgoc in sorted_operon_ordering:
	    ofile2 = open_file(out_file_name,'a')
	    ofile2.write(str(ordering) + '\t' + str(op_fgoc) + '\n')

	ofile2.close()
    return True


def main(argument):
    """
    @function: Takes in the argument which consist of two parts - operon_name and shuffle_status (shuffle or normal). It uses this information and accordingly calculates the op_fgoc values
    """
    
    #(1) Get a key:value dictionary of all the clusters. bsub == It is the alias used for calling python script - bsub_fgoc_all_genome
    cluster_dictionary = bsub.get_cluster_locus_tag_dictionary()
    
    # (2) The argument passed contains two items. Item1 == Operon name; Item2 == shuffle or normal. Obtain details of the operon queried and the nature of the task (shuffle operons or normal)
    try:
	assert len(argument) == 2
	assert isinstance(argument[0],str)
	assert isinstance(argument[1],str)
	if (argument[1] == 'shuffle') or (argument[1] == 'normal'):
	    pass
	else:
	    print "Bad second argument. Exiting..."
	    sys.exit(1)
    except (AssertionError):
	print "Invalid Argument. Exiting..."
	sys.exit(1)
	
    operon_name,shuffle_status = argument
    print "Processing %s ..."%(operon_name)

    # Check shuffle status - normal or shuffle. Depending on choice pass operon query differently
    if shuffle_status == 'normal':
	out_file_name = os.getcwd() + '/tryptophan_study/trp_normal.txt'
	compute_within_operon_walks(operon_name,cluster_dictionary,out_file_name)
    else:
	out_file_name = os.getcwd() + '/bacillus_shuffled/' + operon_name + '_shuffled.sort'
	compute_shuffled_op_fgoc(operon_name,cluster_dictionary,out_file_name)
	
    return
    
if __name__ == '__main__':
    # Calling database using sqlite module
    conn = sqlite3.connect('trial_all_orgs.db')
    cursor = conn.cursor()
    
    main(sys.argv[1:])
    
    cursor.close()
    conn.commit()
    conn.close()
    
    import time
    print "Script - operon_shuffling.py %s \t Completed \t %s"%(main, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
