#! usr/bin/bash

# Making a dictionary of sequence_id and groups
# Notes on dictionary - Key:value; Key is unique. Can get value for a key by having dictionary[key]. 
# For this case the group_id will be the value and the sequence_id the key (which will be unique)

import re

#Initializing
group_dictionary = {}
count = 0



def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	sys.exit(0)

    return file

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
   
def divide_gene_pair_walks(file_name):
    """ This function takes in the column of all tried gene pairs. It differentiates the two into adjacent gene pair groups and non adjacent gene pairs. The adjacent gene pairs are also in the column preceding it in the file. As an identifier it returns back the name of the walk. 
    """
    #Initializing
    parsed_gene_pairs = []
    
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    
    for line in lines:
	# Initializing
	adjacent_gene_pairs = []
	non_adjacent_gene_pairs = []
	walk_name = []
	intermediate_gene_pair = []
	number_walk = []
	
	split_line = line[:-1].split('\t')
	walk_name, number_walk, adjacent_gene_set, all_gene_set = split_line[1][:-3], split_line[2], split_line[4], split_line[5]
	adjacent_genes = adjacent_gene_set.split(',')

	for gene_pairs in adjacent_genes:
	    if not ((gene_pairs == '') or (gene_pairs == 'Null')):
		gene_pairs = gene_pairs.split(':')
		start_walk, stop_walk = gene_pairs[0], gene_pairs[1]
		start_walk = start_walk.split('(')[0]
		stop_walk = stop_walk.split('(')[0]
		gene_paired = start_walk + ':' + stop_walk
		adjacent_gene_pairs.append(gene_paired)

	#print all_gene_set
	all_genes = all_gene_set.split(',')
	for gene_pairs in all_genes:
	    if not ((gene_pairs == '') or (gene_pairs == 'Null')):
		gene_pairs = gene_pairs.split(':')
		start_walk = gene_pairs[0]
		stop_walk = gene_pairs[1]
		gene_paired = start_walk + ':' + stop_walk
		if not (gene_paired in adjacent_gene_pairs):	#Ensures that the gene_pair isnt already in the other list
		    non_adjacent_gene_pairs.append(gene_paired)
	
	if not (adjacent_gene_pairs): adjacent_gene_pairs.append('Null')
	if not (non_adjacent_gene_pairs): non_adjacent_gene_pairs.append('Null')
	
	intermediate_gene_pair.append(walk_name)
	intermediate_gene_pair.append(number_walk)
	intermediate_gene_pair.append(adjacent_gene_pairs)
	intermediate_gene_pair.append(non_adjacent_gene_pairs)
	
	parsed_gene_pairs.append(intermediate_gene_pair)
	
    
    return parsed_gene_pairs
    
    
def get_details(lines, sequence_id, acc_no, rank):
    """ This function provides the details of the line in question. 
    Parameters passed - lines, sequence_id, accession_number, rank. The latter three are empty but either sequence_id must be passed or acc_no and rank be passed
    Returned - rank, orientation, left_end, right_end
    """
    
    #Initializing
    flag = True

    
    for line in lines:
	if (flag):
	    split_line = line.split('\t')
	    if (sequence_id): 
		if (sequence_id == split_line[0]):
		    sequence_id, rank, orientation, left_end, right_end = split_line[0], int(split_line[1]), split_line[2], split_line[3], split_line[4]
		    flag = False
	    
	    
	    if (rank) and (acc_no):
		if (str(rank) == split_line[1] and acc_no == split_line[6]):
		    sequence_id, rank, orientation, left_end, right_end = split_line[0], int(split_line[1]), split_line[2], split_line[3], split_line[4]
		    flag = False

    
    return sequence_id, rank, orientation, left_end, right_end

def get_distance(orientation_1, left_end_1, right_end_1, orientation_2, left_end_2, right_end_2):
    """ This function computes the distance between two genes
    Parameters passed: 1st gene (orientation, left_end, right_end) and 2nd gene (orientation, left_end, right_end)
    Returns: distance (as integer)
    """
    #Initializing
    
    if not (orientation_1 == orientation_2) : 
	#print "Opposite orientation"
	return 0,'opposite'
    
    if orientation_1 == 'forward' : distance = int(left_end_2)- int(right_end_1)
    if orientation_1 == 'reverse' : distance = int(left_end_1) - int(right_end_2)
    
    return distance,'correct'

def get_group(gene_pair, lines, sequence_group_dictionary):
    """ This function takes in one gene - pair of a non adjacent genepair column and returns back the group of the next gene (rankwise) that follows the first gene
    This also takes into account the orientation. If there is a mismatch with the next gene a Null is returned
    """
    
    #Initializing
    start_walk_seq_id = acc_no = seq_id = split_line = ''
    
    
    start_walk_seq_id = gene_pair.split(':')[0]
    acc_no = start_walk_seq_id.split('|')[0]
    
    start_sequence_id, start_rank, start_orientation, start_left_end, start_right_end = get_details(lines, start_walk_seq_id,None,None)
    if start_orientation == 'forward' :
	next_rank = start_rank + 1
    if start_orientation == 'reverse' :
	next_rank = start_rank - 1

    next_sequence_id, next_rank, next_orientation, next_left_end, next_right_end = get_details(lines,None,acc_no, next_rank)

    try:
	group_id = sequence_group_dictionary[next_sequence_id]
    except KeyError:
	group_id = 'Null'
    
    distance, orientation_status = get_distance(start_orientation, start_left_end,start_right_end, next_orientation, next_left_end, next_right_end)

#    print start_sequence_id,start_orientation, start_left_end,start_right_end, next_sequence_id, next_orientation, next_left_end, next_right_end

    return group_id, distance, orientation_status


def get_other_groups(walk_list, sequence_group_dictionary):
    """ This function determines a set of all the orthogroups the next gene was in the list of genes that was the starting gene for the walk. It makes a set and returns the number of distinct orthogroups. In that ways one can obtain a count of all A-X pairings. Note that the orientation, gene distance must be checked. If orientation is reverse then decrease the rank by one and ask which orthogroup did the next gene belong to. 
    """
    
    #Initializing
    all_groups = []
    
    ifile = open_file('coding_sequences.txt')
    lines = ifile.readlines()
    ifile.close()
    
    
    for walk in walk_list:
	#Initializing
	group_one_walk = []
	# This is one gene walk from now
	non_adjacent_gene_pairs = walk[2]
	if not (non_adjacent_gene_pairs == ['Null']):
	    for gene_pair in non_adjacent_gene_pairs:
		# This is one gene-pair in the set of non_adjacent_gene_pairs
		group_next, distance_next, orientation_status = get_group(gene_pair,lines, sequence_group_dictionary)
		if not (orientation_status == 'opposite'):
		    group_one_walk.append(group_next + ':' + orientation_status)
	else:
	    group_one_walk = ['Null']
	
	if not (group_one_walk) : group_one_walk = ['Null']
	group_one_walk_set = set(group_one_walk)
	len_one_walk_set = len(group_one_walk_set) + 1	#The 1 includes the fact that the actual walk is accounted
	group_walk_information = [[item, group_one_walk.count(item)] for item in group_one_walk_set]
	
	walk_information = [walk[0], walk[1], group_walk_information, len_one_walk_set]
	all_groups.append(walk_information)
	print walk_information
	
    return all_groups
    
    

if __name__ == '__main__':
    
    file_name_allwalks = 'allwalks_output.txt'
    
    sequenceid_groupid_dict = sequence_group_dictionary()
    
    
    #walk_name, adjacent_gene_pairs, non_adjacent_gene_pairs = divide_gene_pair_walks(file_name_allwalks)
    parsed_walk_list = divide_gene_pair_walks(file_name_allwalks)
    
    all_groups_set = get_other_groups(parsed_walk_list,sequenceid_groupid_dict)
    print all_groups_set[1]
    
    print parsed_walk_list[1][1]

    


## Checking non walk genes from 'allwalks_output.txt' and returning their group_ids

#ifile = open('allwalks_output.txt','r')
#lines = ifile.readlines()
#ifile.close()

#for line in lines:
    #split_line = line[:-1].split('\t')
    #walk_name = split_line[1]


