#! /usr/bin/bash

# This is for the Trp operon prediction. The query will be an operon name. I will look into the ecoli_gene_operonset or rather allwalks_output.txt to determine which of the walks must be considered. Then run the other txt file to pull out all the other kind of pairings. I hoped to lay out something like a similar operons (based on groups) more graphically


# Rephrasing question - I first query Trp operon: It tells me the genes involved, the walks it is involved in, the group_ids of all the member genes

# Step II: I now enter the group_ids in some order - the order being determined is the order of the operon I want to check. Return back the probability of all the walks again. i.e. query = (gp_1,gp_2,gp_3). It now considers these three to be the genes in that group and tells if there is a probability that two genes between the two adjacent groups have members that are adjacent. 

# Files used:	1) 'ecoli_data.txt' - Generated from sqlite command of joining tables. Operon name, sequence_id, sequence_length, rank, gene_id, orientation (forward, reverse), left_end, right_end, group_id
# 		2) 'groups.txt' - Generated from orthomcl. Group id (JEM_#): all sequence ids. Each sequence_id is separated by a space








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
	    ecoli_genes_information.append([sequence_id, rank, orientation, left_end, right_end, group_id])
	    group_order.append(group_id)

    return ecoli_genes_information, group_order

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

def prob_group_order(groups, group_dictionary):
    """
    @parameter input:	groups - This is the list of the groups in the given order
			group_dictionary - this is the key:value pairs which helps identify the group given the sequence_id
    @function:		Determines the probability of the group_order. Obtains the pairwise probability of each pairwise groups
    @parameter output:	probability_gene_order
    """
    # Initializing
    i = 1
    group_pair = []
    group_pairs = []
    
    group_1 = groups[0]
    
    for group_2 in groups:
	if i>1:
	    group_pair = [group_1, group_2]
	    group_pairs.append(group_pair)
	    group_1 = group_2
	i+=1
    
    print group_pairs

    return
    
#class Groups:
    """ Contains information about the group_order which is passed. Performs 
    """

    



if __name__ == '__main__':
    
    
    query_operon = 'trpLEDCBA'
    
    ecoli_genes, group_order = get_ecoli_genes(query_operon)	
    group_dictionary = sequence_group_dictionary()
    print group_order
    print ecoli_genes
    
    prob_group_order(group_order, group_dictionary)
    
    
    