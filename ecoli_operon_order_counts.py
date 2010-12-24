#! usr/bin/python


# This script is near identical to the possible_operons.py. The only minor difference is the looping over of the operon set to generate names of all the operons. Also another minor change is that there is no permutation of the groups done. 
# This is for quantifying gene pairing of the operons - in E.coli. It returns back all the gene pairs that were consistent and those that would be gene pair, thereby indicating another possible gene pair. 

# Step I: I enter the operon_name like 'trpLEDCBA' by looping over the operon list
# Step II: It obtains all the clusters (or groups) associated with each of the gene in this operon. Does a permutation of all the possible order of the operon. 
# Step III: Each of the order set is divided into gene pairs (group_pair) and passed to the class Groups. Also the group_locus_tag dictionary is passed. 
# Step IV: The class (a) looks up the sqlite table for the group and returns back all the sequence_details of all the genes belonging to that group (or cluster) (b) It also uses the first group detail of the gene_pair that did not have adjacent genes and uses the rank and accession number and rank of the first gene to look up and return back the information of the would be adjacent gene if any. The group or cluster information is returned to determine all the group possibility. 

# Sqlite database used - all_trial_orgs.db
# Tables - coding_sequences, protein_clusters

# Files used:	1) 'ecoli_cluster.txt'. - Generated from sqlite command of joining three tables. Operon name, locus tag, sequence_length, rank, orientation, left_end, right_end, cluster_accession, protein_name
# Will not use this	3) 'coding_sequences.txt'- Generated from the sqlite command by obtaining information of the CDS. It has the following columns (tab delimited)
#			(a) locus_tag (but using as sequence_id everywhere to match with the newly made cluster_dictionary, (b) rank, (c) orientation (d) left_end, (e) right_end (f) sequence_length (g) accession_number (i) taxon_id


import numpy as np
import itertools
import os
import sys
import sqlite3
from warnings import filterwarnings
import timeit
import pdb


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
	#print genome_size
	
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
	raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	sys.exit(0)
    
    return file
    
def get_ecoli_genes(query_operon):
    """ 
    @param_input: query_operon: name of the operon in ecoli
    @parm_output: ecoli_genes_information: the name of the genes of ecoli that has the operon. The [[locus_tag, rank, gene_id, orientation, left_end, right_end, group_id],[..],..]
    group_order: [group1, group2, ...]
    @function: Looks up the text file - ecoli_data.txt to obtain information of all ecoli genes that match the operon description
    """
	
    #Initializing
    ecoli_genes_information = []
    group_order = []
    
    ifile = open_file('ecoli_cluster.txt')
    lines = ifile.readlines()
    ifile.close()

    for line in lines:
	if query_operon == line.split('\t')[0]:
	    split_line = line.split('\t')
	    [locus_tag, rank, orientation, left_end, right_end, group_id] = split_line[1], split_line[3], split_line[4], split_line[5], split_line[6], split_line[7]
	    orientation_status = orientation
	    ecoli_genes_information.append([locus_tag, rank, orientation, left_end, right_end, group_id])
	    group_order.append(group_id)
	    
    return ecoli_genes_information, group_order, orientation_status


def prob_group_order(groups, operon,group_dictionary,ofile_handle):
    """
    @param groups: This is the list of the groups in the given order
    @param group_dictionary: this is the key:value pairs which helps identify the group given the sequence_id
    @param operon: name of the operon
    @function:	Obtains the pairwise probability of each pairwise groups in the group list. 
    @param output: number of conserved gene pairs and total number of possibilities of gene pairs
    """
    # Initializing
    i = 1
    group_pair = []
    group_pairs = []
    walk_number = 1
    total_groups = []
    total_conserved_gene_pairs = 0
    nbr_group_pairs = 0
    #out_file_name = 'walk_in_ecoli_operons.txt'
    out_file_name = os.getcwd() + '/shuffled_operons/' + operon + '_shuffled.txt'
    #if os.path.exists(out_file_name):
	#os.remove(out_file_name)    
    
    group_1 = groups[0]
    
    for group_2 in groups:
	if i>1:
	    group_pair = [group_1, group_2]
	    group_pairs.append(group_pair)
	    group_1 = group_2
	i+=1

    for group_pair in group_pairs:
	gene_pair = Groups(group_pair,group_dictionary)
	#gene_a_occurrences, gene_b_occurrences, number_conserved_gene_pair, difference_details = gene_pair.find_adjacent_pairs()	#Old better one
	# Changing to ignore the number of all_non_became_groups, so getting additional non_adjacent_gene_pairs from the find_adjacent_pairs function
	gene_a_occurrences, gene_b_occurrences, number_conserved_gene_pair, difference_details, non_adjacent_gene_pairs = gene_pair.find_adjacent_pairs()
	#non_adjacent_gene_pairs, all_non_became_groups = gene_pair.get_details_non_adjacent_genes()
	
	# non_adjacent_gene_pairs = These refer to those in which gene_a is in group_a and there is gene_b in group_b but the gene_a and gene_b aren't next to one another  
	# For all_non_became_groups = gene_a doesn't have an adjacent gene_b in group_b. So I look to see if the next gene of gene_a in the same organism belongs to which group. So I get a list of all possible groups this gene_a was a neighbour of apart from group_b
	#all_non_became_groups.append(group_pair[1])
	number_non_adjacent_gene_pairs = len(non_adjacent_gene_pairs)
	#total_groups = all_non_became_groups
	#number_total_groups = len(set(total_groups))
	total_conserved_gene_pairs += number_conserved_gene_pair
	
	ofile_handle.write("Operon %s  walk : %d \n"%(operon,walk_number))	#Walk number and operon name - tells the step number in that operon
	ofile_handle.write("# Conserved gene pair \t %d \n"%(number_conserved_gene_pair))
	ofile_handle.write("# Occurrences of gene 1 \t %d \n"%(gene_a_occurrences))
	ofile_handle.write("# Occurrences of gene 2 \t %d \n"%(gene_b_occurrences))
	ofile_handle.write("# Non adjacent gene pairs \t %d \n"%(number_non_adjacent_gene_pairs))
	#ofile_handle.write("# Total possible adjacent genes \t %d \n"%(number_total_groups))

	#Old information_for_file: has number_total_groups
	#information_for_file = '\t'.join(item for item in [operon, str(walk_number), str(number_conserved_gene_pair),str(number_total_groups), str(gene_a_occurrences), str(gene_b_occurrences), str(difference_details), str(non_adjacent_gene_pairs),'\n'])
	#eval(string) will convert the string into a list
	
	information_for_file = '\t'.join(item for item in [operon, str(walk_number), str(number_conserved_gene_pair), str(gene_a_occurrences), str(gene_b_occurrences), str(difference_details), str(non_adjacent_gene_pairs),'\n'])
	
	walk_number += 1
	
	#Appending information to file - operon_walk_gene_pairs.txt
	ofile = open_file(out_file_name,'a')
	ofile.write(information_for_file)

    ofile.close()
    number_total_groups = 10	#This is the new part just to pass a value. This is never used. 
    return total_conserved_gene_pairs, number_total_groups

def make_dir(dir_name,path = os.getcwd()):
    """This function makes a directory in the path passed to the function. If no path is passed then it takes the current working directory as default. If the directory already exists then it doesnt do anything and returns back the directory path. In both cases the new directory path is returned back """
    
    #Initializing
    new_dir_path = ''
    
    new_dir_path = path + '/' + dir_name
    try:
	os.mkdir(new_dir_path)
    except OSError:
	    #print "Directory already exists. Overwriting on files in directory if any!"
	pass
	
    return new_dir_path
	


def operon_walks(operon, group_dictionary):
    """
    @param operon - This is the query operon
    @param group_dictionary - This is the dictionary created - values - cluster_accession: key - locus_tag
    @function:		Calls the function get_ecoli_genes(operon) to get the details of the sequence of the ecoli_genes belonging to the operon.
			It then checks the group orientation and then calls the function prob_group_order(group_order, operon, group_dictionary) and obtains the number_gene_pairs (that were adjacent) and the total number of gene pairs that are possible. 
    @param output - None
    """

    #Initializing
    group_order = []
    all_gene_pairs = []
    total_conserved_gene_pairs = 0
    total_possible_groups = []
    ecoli_genes, group_order, operon_orientation_status = get_ecoli_genes(operon)	
    dir_name = 'shuffled_operons'
    new_dir_path = make_dir(dir_name)
    ofile_name_out = new_dir_path + '/' + operon + '.shuffled'
    ofile = open_file(ofile_name_out,'a')
    
    if operon_orientation_status == 'reverse':
	group_order.reverse()
    
    permutation_list = list(itertools.permutations(group_order))	# This is only for the special case where I want to analyse all the different permutations of the group order. This mustnt be applied in the case of regualar within ecoli_operons. 
    
    for group_order in permutation_list:
	try:
	    assert len(group_order) > 1
	    ofile.write("Group Order : %s \n"%(str(group_order)))
	    total_conserved_gene_pairs, total_possible_groups = prob_group_order(group_order, operon, group_dictionary,ofile)
	    ofile.write("Total conserved gene_pairs in arrangement = %d \n "%(total_conserved_gene_pairs))
	    ofile.write("/// \n")
	    print "%s done"%(str(group_order))
	except AssertionError:
	    os.remove(ofile_name_out)
	    pass
    
    
    ofile.close()
    
    return
    
def get_cluster_locus_tag_dictionary():
    """
    @param input:	None
    @function:		This function obtains the key:value dictionary pair for cluster_accession as the value and locus_tag being the key. 
    It uses the file - cluster_locus_tag.txt. This has two columns - cluster accession and locus_tag. Just converting it to key:value pair for faster identification
    @param output:	cluster_dictionary -It returns back the dictionary list. I am keeping the name cluster_dictionary just for this function as I dont want to mess up the class which everywhere uses group_dictionary
    """
    
    # Initializing
    file_name = 'cluster_locus_tag.txt'
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    cluster_dictionary = {}
    
    for line in lines:
	cluster_accession, locus_tag = line.split('\t')[0],line.split('\t')[1]
	cluster_dictionary[locus_tag] = cluster_accession
	
    return cluster_dictionary
    
def between_operon_walks(group_dictionary):
    """
    @param group_dictionary:	This is the key:value dictionary pair for cluster_accession as the value and locus_tag being the key
    @function: This function reads the 'ecoli_cluster.txt' file and watches when the operon name changes and passes the groups of these two as the group_order to the class
    @parm output:		None
    """
    #Initializing
    ifile = open_file('ecoli_cluster.txt')
    lines = ifile.readlines()
    ifile.close()
    i = 1
    all_group_pairs = []
    gene_pair = []
    all_non_became_groups = []
    walk_number = 1
    all_information = ''
    
    group_1 = (lines[0].split('\t')[7])
    operon_1 = lines[0].split('\t')[0]

    for line in lines:
	if i>1:
	    split_line = line.split('\t')
	    operon_2 = split_line[0]
	    group_2 = (split_line[7])
	    if not operon_1 == operon_2:
		operon_pair = [operon_1,operon_2]
		group_pair = [group_1,group_2]
		all_group_pairs.append([operon_pair,group_pair])	#Note this difference, there is an [operon_pair],[group_pair]
		
	    group_1 = group_2
	    operon_1 = operon_2
	i+=1


    for all_info_pair in all_group_pairs:
	
	operon_pair = all_info_pair[0]
	group_pair = all_info_pair[1]

	gene_pair = Groups(group_pair,group_dictionary)
	gene_a_occurrences, gene_b_occurrences, number_conserved_gene_pair, difference_details = gene_pair.find_adjacent_pairs()

	non_adjacent_gene_pairs, all_non_became_groups = gene_pair.get_details_non_adjacent_genes()
	# non_adjacent_gene_pairs = These refer to those in which gene_a is in group_a and there is gene_b in group_b but the gene_a and gene_b aren't next to one another  
	# For all_non_became_groups = gene_a doesn't have an adjacent gene_b in group_b. So I look to see if the next gene of gene_a in the same organism belongs to which group. So I get a list of all possible groups this gene_a was a neighbour of apart from group_b
	all_non_became_groups.append(group_pair[1])
	number_non_adjacent_gene_pairs = len(non_adjacent_gene_pairs)
	total_groups = all_non_became_groups
	number_total_groups = len(set(total_groups))
	
	print "Operon Pair %s  walk : %d "%(str(operon_pair),walk_number)	#Walk number and operon name - tells the step number in that operon
	print "# Conserved gene pair \t %d"%(number_conserved_gene_pair)
	print "# Occurrences of gene 1 \t %d"%(gene_a_occurrences)
	print "# Occurrences of gene 2 \t %d"%(gene_b_occurrences)
	print "# Non adjacent gene pairs \t %d"%(number_non_adjacent_gene_pairs)
	print "# Total possible adjacent genes \t %d"%(number_total_groups)
	print "///"
	
	information_for_file = '\t'.join(item for item in [str(operon_pair), str(walk_number), str(number_conserved_gene_pair), str(number_total_groups),  str(gene_a_occurrences), str(gene_b_occurrences) ,str(difference_details), str(non_adjacent_gene_pairs),'\n'])	
	#eval(string) will convert the string into a list

	all_information += information_for_file
	walk_number += 1
	
    #Appending information to file - walk_between_ecoli_operons.txt
    out_file_name = 'walk_between_ecoli_operons.txt' 
    if os.path.exists(out_file_name):
	os.remove(out_file_name)
    
    ifile = open_file(out_file_name,'w')
    ifile.write(all_information)
    
    ifile.close()

    return


if __name__ == '__main__':
    
    # Calling database using sqlite3 module
    conn = sqlite3.connect('trial_all_orgs.db')
    cursor = conn.cursor()

    ifile = open_file('ecoli_cluster.txt')
    lines = ifile.readlines()
    ifile.close()

    all_operon_list =([line.split('\t')[0] for line in lines]) 
    operon_list = set(all_operon_list)
    
    cluster_dictionary = get_cluster_locus_tag_dictionary()
    total_group_pairs = 0
    #for operon in operon_list:
	#operon_walks(operon,cluster_dictionary)
	
	
    #between_operon_walks(cluster_dictionary)

    #for operon in operon_list:
	#operon_walks(operon,cluster_dictionary)
	
    operon = sys.argv[1]

    operon_walks(operon,cluster_dictionary)
    
    import time
    print "Script - ecoli_operon_order_counts.py %s \t Completed \t %s"%(operon, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    