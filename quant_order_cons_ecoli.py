#! /usr/bin/bash

# This is for the Trp operon prediction. The query will be an operon name. I will look into the ecoli_gene_operonset or rather allwalks_output.txt to determine which of the walks must be considered. Then run the other txt file to pull out all the other kind of pairings. I hoped to lay out something like a similar operons (based on groups) more graphically


# Rephrasing question - I first query Trp operon: It tells me the genes involved, the walks it is involved in, the group_ids of all the member genes

# Step II: I now enter the group_ids in some order - the order being determined is the order of the operon I want to check. Return back the probability of all the walks again. i.e. query = (gp_1,gp_2,gp_3). It now considers these three to be the genes in that group and tells if there is a probability that two genes between the two adjacent groups have members that are adjacent. 

# Files used:	1) 'ecoli_data.txt' - Generated from sqlite command of joining tables. Operon name, sequence_id, sequence_length, rank, gene_id, orientation (forward, reverse), left_end, right_end, group_id
# 		2) 'groups.txt' - Generated from orthomcl. Group id (JEM_#): all sequence ids. Each sequence_id is separated by a space
# 		3) 'coding_sequences.txt' - Generated from the sqlite command by obtaining information of the CDS. It has the following columns (tab delimited)
#			(a) sequence_id, (b) rank, (c) orientation (d) left_end, (e) right_end (f) sequence_length (g) accession_number (i) taxon_id

import numpy as np

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
	
    def get_details(self,sequence_id, acc_no, rank):
	""" This function provides the details of the line in question. 
	Parameters passed - sequence_id, accession_number, rank. The latter two can be empty but either sequence_id must be passed or acc_no and rank be passed
			    File opened: 'coding_sequences.txt'
	Returned - sequence_id, rank, orientation, left_end, right_end, accession_number
	"""
	
	#Initializing
	flag = True
	ifile = self.open_file('coding_sequences.txt')
	lines = ifile.readlines()
	ifile.close()
	orientation = left_end = right_end = accession_number = ''
	
	for line in lines:
	    if (flag):
		split_line = line.split('\t')
		if (sequence_id): 
		    if (sequence_id == split_line[0]):
			sequence_id, rank, orientation, left_end, right_end, accession_number = split_line[0], int(split_line[1]), split_line[2], int(split_line[3]), int(split_line[4]), split_line[6]
			flag = False
		    
		if (rank) and (acc_no):
		    if (str(rank) == split_line[1] and acc_no == split_line[6]):
			sequence_id, rank, orientation, left_end, right_end, accession_number = split_line[0], int(split_line[1]), split_line[2], int(split_line[3]), int(split_line[4]), split_line[6]
			flag = False
		
		
	#print sequence_id, rank, orientation, left_end, right_end, accession_number
	return [sequence_id, rank, orientation, left_end, right_end, accession_number]
			    

    def find_group_members(self,group_name):
	""" Input:	group_name = example: 'JEM_5'
	    Function:	1. Adds a ':' to the group_name
			2. Look for the group name in 'groups.txt' and return back the sequence_id of the members (they are separated by ' ')
	    Output:	List of all the sequence_ids of the members in that group
	"""
	#Initializing
	group_members = []
	acc_no = ''
	rank = ''
	
	ifile = open_file('groups.txt')
	lines = ifile.readlines()
	ifile.close()
	
	for line in lines:
	    if line.startswith(group_name+':'):
		split_line = line.split(' ')	#Space is all that separates the different sequenceids
		(split_line.pop(0))[:-1]	#It has a ':' in the end
		
		for sequence in split_line:
		    sequence_id = sequence.replace('\n','')
		    sequence_details = self.get_details(sequence_id,acc_no,rank)
		    group_members.append(sequence_details)

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
	    rank_difference = rank_a - rank_b
	    distance = left_end_a - right_end_b

	return rank_difference, distance


    def get_details_non_adjacent_genes(self):
	"""
	Input:	Taken from self - list of all non_adjacent_gene_pairs that had a rank difference not equal 1. 
	Function:	1. Gets the first gene - gene_a from each of the member in the list. 
			2. It obviously belongs to the group_a. It checks for the gene next to it in the genome. Gets back details of the gene
			3. Checks if the orientation is the same. If different then - the gene and group is called 'null'
			4. If orientation is correct - it obtains the distance between them and the group which the next gene belonged to.
			5. Sometimes there is no group affliation in which case a keyerror is returned. the group is considered 'null'
			5. Add to list group. Reiterate the steps 1 - 5
	@out non_adjacent_gene_pair list: this includes [gene_a, gene_b, rank_difference, distance]
	@out all_non_became_groups: this is a list (non redundant) where the list of all groups (that became when next gene was looked at) 
	"""
	# Initializing
	non_adjacent_gene_pairs = self.non_adjacent_gene_pairs
	sequence_id = rank_a = orientation_a = left_end_a = right_end_a = acc_no_a = ''
	
	next_gene = []
	all_non_became_groups = []
	set_all_non_became_groups = []
	
	try:	# In some cases there are no non_adjacent_gene_pairs. So it is better i send it back without doing any calculation
	    assert non_adjacent_gene_pairs
	except AssertionError:
	    return ['No non_adjacent_gene_pair'],[]
	    
	
	for gene_pair in non_adjacent_gene_pairs:
	    gene_a = gene_pair[0]
	    rank_a,orientation_a,left_end_a,right_end_a, acc_no_a = int(gene_a[1]), gene_a[2], int(gene_a[3]), int(gene_a[4]), gene_a[5]
	    try:
		print rank_a, orientation_a, left_end_a, right_end_a
		assert ((rank_a) and (orientation_a) and (left_end_a) and (right_end_a) and (acc_no_a))
	    
		if orientation_a == 'forward':
		    next_gene = self.get_details(sequence_id, acc_no_a, (rank_a + 1))
		if orientation_a == 'reverse':
		    next_gene = self.get_details(sequence_id, acc_no_a, (rank_a - 1))
		
		sequence_id, rank, orientation, left_end, right_end, accession_number = next_gene[0], next_gene[1], next_gene[2], next_gene[3], next_gene[4], next_gene[5]
		
		assert (orientation_a) and (orientation)
		
		orientation_status = self.check_orientation(orientation_a, orientation)
		
		if not (orientation_status == 'opposite'):
		    rank_difference, distance = self.compute_rank_difference(orientation_status, gene_a, next_gene)
		    
		if orientation_a == orientation:
		    try:
			next_group = self.group_dictionary[sequence_id]
		    except KeyError:
			next_group = 'Null'
		else:
		    next_group = 'Null'
	
	    except AssertionError:
		next_group = 'Null'
		
	    all_non_became_groups.append(next_group)
	    
	return non_adjacent_gene_pairs, all_non_became_groups

    def find_adjacent_pairs(self):
	""" This function finds all the adjacent_pairs for the sent group_pair. 
	1. Look into group_1 and obtain the list of all sequence_ids in that group along with necessary details of rank, orientation, position etc
	(Obtains the details from the file - coding_sequences.txt - (sequence_id, rank, orientation, left_end, right_end, sequence_length, accession_nbr, taxon_id))
	2. Do the same with group_2
	3. Do pairwise run through both the groups - check if they have the same accession number, same orientation and obtain (is rank difference 1 (or -1)
	4. If they match that - obtain the distance between the two genes and the gene_pair (via their sequence_id_1:sequence_id_2:distance:taxon_id
	"""
		
	# Initializing
	group_a = self.group_1
	group_b = self.group_2
	count = 0
	difference_details = []
	non_adjacent_gene_pairs = []
	
	
	self.group_a_members = self.find_group_members(group_a)
	self.group_b_members = self.find_group_members(group_b)
	
	for gene_a in self.group_a_members:
	    rank_a,orientation_a,left_end_a,right_end_a, acc_no_a = int(gene_a[1]), gene_a[2], int(gene_a[3]), int(gene_a[4]), gene_a[5]
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
	
	return count, difference_details


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
    @parm_output: ecoli_genes_information: the name of the genes of ecoli that has the operon. The [[sequence_id, rank, gene_id, orientation, left_end, right_end, group_id],[..],..]
    group_order: [group1, group2, ...]
    @function: Looks up the text file - ecoli_data.txt to obtain information of all ecoli genes that match the operon description
    """
	
    #Initializing
    ecoli_genes_information = []
    group_order = []
    
    ifile = open_file('ecoli_data.txt')
    lines = ifile.readlines()
    ifile.close()
    
    for line in lines:
	if query_operon == line.split('\t')[0]:
	    split_line = line.split('\t')
	    [sequence_id, rank, orientation, left_end, right_end, group_id] = split_line[1], split_line[3], split_line[5], split_line[6], split_line[7], split_line[8][:-1]
	    orientation_status = orientation
	    ecoli_genes_information.append([sequence_id, rank, orientation, left_end, right_end, group_id])
	    group_order.append(group_id)
	    
    return ecoli_genes_information, group_order, orientation_status
		
def sequence_group_dictionary():
    """ This function creates a dictionary of sequence_id and group_id with sequence_id (key): group_id (value). The details are taken from the groups.txt file created by orthomcl. 
    Returns: sequenceid_group_dictionary
    """
    #Initializing
    group_dictionary = {}
    split_line = line = sequence_id = sequence = ''
    
    ifile = open_file('groups.txt')
    lines = ifile.readlines()
    ifile.close()
    
    for line in lines:
	split_line = line.split(' ')	#Space is all that separates the different sequenceids
	group_name = (split_line.pop(0))[:-1]	#It has a ':' in the end
	
	for sequence in split_line:
	    sequence_id = sequence.replace('\n','')	#Some seq id has this annoying carriage return
	    group_dictionary[sequence_id] = group_name
	    
    return group_dictionary
			    
def prob_group_order(groups, operon,group_dictionary):
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
    
    group_1 = groups[0]
    
	
    for group_2 in groups:
	if i>1:
	    group_pair = [group_1, group_2]
	    group_pairs.append(group_pair)
	    group_1 = group_2
	i+=1

    for group_pair in group_pairs:

	gene_pair = Groups(group_pair,group_dictionary)
	number_conserved_gene_pair, difference_details = gene_pair.find_adjacent_pairs()

	non_adjacent_gene_pairs, all_non_became_groups = gene_pair.get_details_non_adjacent_genes()
	print group_pair, non_adjacent_gene_pairs, all_non_became_groups
	all_non_became_groups.append(group_pair[1])
	total_groups = all_non_became_groups
	print total_groups
	number_total_groups = len(set(total_groups))
	
	print "Operon : %s  walk : %d "%(operon,walk_number)	#Walk number and operon name - tells the step number in that operon
	print "# Conserved gene pair: %d"%(number_conserved_gene_pair)
	print "# Total number of possible adjacent genes: %d"%(number_total_groups)
	

	information_for_file = '\t'.join(item for item in [operon, str(walk_number), str(number_conserved_gene_pair),str(number_total_groups), str(difference_details),  str(non_adjacent_gene_pairs),'\n'])	#eval(string) will convert the string into a list
	
	walk_number += 1
	print information_for_file
	
	# Appending information to file - operon_walk_gene_pairs.txt
	ifile = open_file('walk_in_operon.txt','a')
	ifile.write(information_for_file)

    ifile.close()
    return number_conserved_gene_pair, number_total_groups
    
def operon_walks(group_dictionary):
    """
    Input:	None
    Function:	1.Opens the file - 'ecoli_data.txt' and makes a list of operons
    """

    #Initializing
    ecoli_genes_information = []
    group_order = []
    all_gene_pairs = []

    ifile = open_file('ecoli_data.txt')
    lines = ifile.readlines()
    ifile.close()

    operon_list = set([line.split('\t')[0] for line in lines])
    
    for operon in operon_list:
		
	ecoli_genes, group_order, operon_orientation_status = get_ecoli_genes(operon)	

	if operon_orientation_status == 'reverse':
	    group_order.reverse()
	
	try:
	    assert len(group_order) > 1
	    number_gene_pairs, total_gene_pairs = prob_group_order(group_order, operon, group_dictionary)
	    all_gene_pairs.append(number_gene_pairs)
	except AssertionError:
	    pass
	

    
    return 

def between_operon_walks(group_dictionary):
    """
    @function: This function reads the 'ecoli_data.txt' file and watches when the operon name changes and passes the groups of these two as the group_order to the class
    """
    #Initializing
    ifile = open_file('ecoli_data.txt')
    lines = ifile.readlines()
    ifile.close()
    i = 1
    all_group_pairs = []
    gene_pair = []
    all_non_became_groups = []
    walk_number = 1
    all_information = ''
    
    group_1 = (lines[0].split('\t')[8])[:-1]	#The stupid return character
    operon_1 = lines[0].split('\t')[0]
    for line in lines:
	if i>1:
	    split_line = line.split('\t')
	    operon_2 = split_line[0]
	    group_2 = (split_line[8])[:-1]
	    if not operon_1 == operon_2:
		operon_pair = [operon_1,operon_2]
		group_pair = [group_1,group_2]
		all_group_pairs.append([operon_pair,group_pair])
		
	    group_1 = group_2
	    operon_1 = operon_2
	i+=1

    for all_info_pair in all_group_pairs:
	
	operon_pair = all_info_pair[0]
	group_pair = all_info_pair[1]
	gene_pair = Groups(group_pair,group_dictionary)
	number_conserved_gene_pair, difference_details = gene_pair.find_adjacent_pairs()
	
	non_adjacent_gene_pairs, all_non_became_groups = gene_pair.get_details_non_adjacent_genes()
	all_non_became_groups.append(group_pair[1])

	total_groups = all_non_became_groups
	number_total_groups = len(set(total_groups))
	
	information_for_file = '\t'.join(item for item in [str(operon_pair), str(walk_number), str(number_conserved_gene_pair), str(number_total_groups),  str(difference_details),  str(non_adjacent_gene_pairs),'\n'])	#eval(string) will convert the string into a list
	
	print information_for_file
	all_information += information_for_file
	walk_number += 1
	
    # Appending information to file - operon_walk_gene_pairs.txt
    ifile = open_file('walk_between_operon.txt','w')
    ifile.write(all_information)
    
    ifile.close()

    return
    

if __name__ == '__main__':
    
    
    ##query_operon = 'trpLEDCBA'
    ##query_operon = 'leuLABCD'
    #query_operon = 'mraZ-rsmH-ftsLI-murEF-mraY-murD-ftsW-murGC-ddlB-ftsQAZ-lpxC'
    group_dictionary = sequence_group_dictionary()
    #operon_walks(group_dictionary)
    between_operon_walks(group_dictionary)
    