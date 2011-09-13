#! usr/bin/python

from __future__ import division
import os
import sys
try:
   import cPickle as pickle
except:
   import pickle
import shelve
import re

# I am making a big matrix of all operon genes (including non-operonic); So 1-2 would mean a score of Homolog(1) -> Homolog(2) -- the FGOC score. I will first make the matrix as a whole and then divide it into 2, representing the forward and the reverse strand. 

pkl_dir = os.getcwd() + '/shelve_files/'
cluster_dictionary = shelve.open(pkl_dir + 'ecoli_homolog_locus_tag:group_dictionary')
#group_all_locus_dictionary = shelve.open(pkl_dir + 'ecoli_homolog_group:all_locus_tag_dictionary')	# Needed
#locus_cds_dictionary = shelve.open(pkl_dir + 'locus:cds_information_dictionary')


# Unpickle the dictionary file
file_name = 'acc_nbr_genes_dict.pkl'
pkl_dir = os.getcwd() + '/pkl_files/'
pkl_file = pkl_dir + file_name
ifile = open(pkl_file)
acc_nbr_genes = pickle.load(ifile)
ifile.close()

class Groups:
    """ Contains information about the group_order which is passed. Computes the number of conserved walks, list of all possible gene_pairs
    """
    # Shelves imported : cluster_dictionary,group_all_locus_dictionary,locus_cds_dictionary
    def __init__(self, group_pair):
	self.group_1, self.group_2 = group_pair[0], group_pair[1]
	self.group_dictionary = cluster_dictionary
	self.group_all_locus_dictionary = group_all_locus_dictionary
	self.locus_cds_dictionary = locus_cds_dictionary
		
    def open_file(self,name_file, open_status = 'r'):
	""" This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """
    
    #Opening and reading/writing tof_he passed file_name """
	try:
	    file = open(name_file,open_status)
	except IOError,err:	# If the file cannot be opened i am exiting
	    raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	    sys.exit(0)
	
	return file
   
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
	
	[genome_size, number_genes] = cursor.fetchone()
	[genome_size, number_genes] = acc_nbr_genes_dict[acc_no_a]
	
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
	# NOT USING THIS FUNCTION IN THIS PROGRAMME
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

    def get_group_members(self,group):
	"""
	@param group: This is the name of the group
	@function: Open pickle dictionary of group: [all locus_tag]; For each of the locus_tag obtain information. Append the information of the group_name to it as well. 
	@return [[locus_tag, rank, orientation, left_end, right_end, accession_number, group_name],[..]..]
	"""
	# Initializing
	all_cds_info = []
	count = 0
	# (1) Opens the pickle_dictionary of group:[all locus_tag]
	try:
	    all_locus_tag = self.group_all_locus_dictionary[group]
	except KeyError:
	    return all_cds_info
	    
	# (2) Iterate through each locus_tag and obtain the information from the locus_cds_dict
	for locus in all_locus_tag:
	    cds_info = self.locus_cds_dictionary[locus]
	    cds_info.append(group)
	    if len(cds_info) == 7:
		all_cds_info.append(cds_info)
	    else:
		print all_locus_tag
		print locus
		sys.exit(1)
		count+=1

	return all_cds_info


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
	
	#self.group_a_members = self.find_group_members(group_a)
	#self.group_b_members = self.find_group_members(group_b)
	self.group_a_members = self.get_group_members(group_a)	 # This uses dictionaries to retrieve the data
	self.group_b_members = self.get_group_members(group_b)
	
	gene_a_occurrences = len(self.group_a_members)
	gene_b_occurrences = len(self.group_b_members)
	lone_gene_occurrences = abs(gene_a_occurrences - gene_b_occurrences)	#This is the number of genes, which was present in one group for an organism, but for the same organism, the second group didnt have any genes. So the would be group should be considered Null

	for gene_a in self.group_a_members:
	    rank_a,orientation_a,left_end_a,right_end_a, acc_no_a = gene_a[1], gene_a[2], gene_a[3], gene_a[4], gene_a[5]
	    

	    for gene_b in self.group_b_members:
		rank_b,orientation_b,left_end_b,right_end_b, acc_no_b = int(gene_b[1]), gene_b[2], int(gene_b[3]), int(gene_b[4]), gene_b[5]
		if acc_no_a == acc_no_b:	#Ensures only the same species is checked
		    #orientation_status = self.check_orientation(orientation_a, orientation_b)
		    #pdb.set_trace()
		    orientation_status = 'forward'	# This is just a fix; I want to just see the difference between two genes.
		    rank_difference, distance = self.compute_rank_difference(orientation_status, gene_a, gene_b)
		    if rank_difference == 1:
			difference_details.append([gene_a,gene_b, distance])
			count+=1
		    else:	#All those gene_a that don't have a genepair in gene_b
			non_adjacent_gene_pairs.append([gene_a,gene_b, rank_difference, distance])
	
	self.non_adjacent_gene_pairs = non_adjacent_gene_pairs
	
	if not (difference_details) : difference_details = ['No Adjacent gene pair']

	return gene_a_occurrences,gene_b_occurrences, count, difference_details, non_adjacent_gene_pairs
#**********************************#
# Class Groups Ends 




def pickle_file(file_name, contents):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    pkl_dir = os.getcwd() + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ofile = open(pkl_file,'wb')
    pickle.dump(contents, ofile)
    ofile.close()
    
    return True
    
def unpickle_file(file_name):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    pkl_dir = os.getcwd() + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ifile = open(pkl_file)
    contents = pickle.load(ifile)
    ifile.close()
    
    return contents


def calculate_fgoc_score(pair):
    """
    @param pair: This is the group pair that is passed either from the function or from the shell line.
    @function: It calculates the FGOC score and outputs the relevant details into file (appendable)
    """
    
    # (1) Initialize the Class groups
    gene_pair = Groups(pair)
    # Class Groups is called and the details of the find_adjacent_pairs, difference_details etc are obtained. Importantly the became_adjacent gene pair i.e. of the gene_pairs which were not adjacent, the details of what is the gene next to it. It is computationally a little more time consuming but is there in the class Groups and can be modified and used. 
    # (2) Obtain all the computed values from the function in class
    gene_a_occurrences, gene_b_occurrences, number_conserved_gene_pair, difference_details, non_adjacent_gene_pairs = gene_pair.find_adjacent_pairs()
    #*********Modification for Ecoli - homolog group *****
    non_adjacent_acc_nbrs = set([(item[0])[5] for item in non_adjacent_gene_pairs])
    print "gene occurrences : %d %d"%(gene_a_occurrences, gene_b_occurrences)
    # (3) Obtain Accession numbers of all genes in difference_details
    if len(difference_details) == 1:
	number_conserved_gene_pair = 0
    else:
	acc_nbrs = set([(item[0])[5] for item in difference_details])
	# The number of conserved_gene_pair is always the number of non - redundant acc_nbrs set
	number_conserved_gene_pair = len(acc_nbrs)
	# Count the non_adj_acc_nbrs that was not in the conserved_acc_nbrs to get a count of number of non_adjacent_acc_nbrs
	count = 0
	for non_adj_org in non_adjacent_acc_nbrs:
	    if not non_adj_org in acc_nbrs:
		count+=1
	number_non_adjacent_gene_pairs = count
    
    # (4) Calculate the fgoc score (Remember to import _future_ division)
    number_non_adjacent_gene_pairs = len(non_adjacent_acc_nbrs)
    try:
	fgoc_score = number_conserved_gene_pair / (number_conserved_gene_pair + number_non_adjacent_gene_pairs)
    except ZeroDivisionError:
	fgoc_score = 0
    print "Processing .. %s"%(pair)
    print number_conserved_gene_pair, number_non_adjacent_gene_pairs
    print fgoc_score

    # (5) Writing the information to the file. 
    information_for_file = '\t'.join(item for item in [str(pair), str(fgoc_score), str(number_conserved_gene_pair), str(number_non_adjacent_gene_pairs), str(gene_a_occurrences), str(gene_b_occurrences), str(difference_details), str(non_adjacent_gene_pairs),'\n'])
    
    file_name = 'all_pairs_fgoc.txt'
    ofile = open(file_name,'a')
    ofile.write(information_for_file)
    ofile.close()

    # (6) Make the mclInput file as well
    file_name = 'mclInput_pairs_fgoc.txt'
    ofile = open(file_name,'a')
    info_line = '\t'.join(item for item in [str(pair[0]), str(pair[1]), str(fgoc_score), '\n'])
    ofile.write(info_line)
    ofile.close()

    return True


def get_group_pairs(groups):
    """
    @param groups: This is a list of groups in order. 
    @function: Provides a pair of these groups (order maintained) to get a list of all group pairs
    @return: List of pairs
    """
    
    #Initializing
    group_1 = group_2 = ''
    group_pairs = []
    pair = []
    i = 1
    
    group_1 = groups[0]
    
    for group_2 in groups:
	if i>1:
	    group_pairs.append([group_1, group_2])
	    group_1 = group_2
	i+=1
	    
    return group_pairs


def make_pairwise_groups(group_list):
    """
    @param group_list: a list of all groups in E.coli (homolog of each gene)
    @function: Makes all pairwise groups in the E.coli genome
    """
    
    # Initializing
    all_pairwise_groups = []
    
    for group_1 in group_list:
	for group_2 in group_list:
	    print group_1, group_2
	    all_pairwise_groups.append([group_1,group_2])
	    
    print len(all_pairwise_groups)
    return all_pairwise_groups

def main(argument):
    
    if not argument[0] == '2':
	pair1 = pair2 = 'Null'
	[argv] = argument
    else:
	[argv,pair1,pair2] = argument
    print "In main"
    if argv == '1':
	# (1) Open the file - ecoli_homolog_data.txt and obtain ecoli_homolog pairs
	file_name = 'ecoli_homolog_data.txt'
	ifile = open(file_name)
	lines = ifile.readlines()
	ifile.close()
	group_list = [line.split('\t')[7][:-1] for line in lines]	#[:-1] is to delete the carriage return
	
	# (2) Make pairwise list of all these genes
	pairwise_group_list = make_pairwise_groups(group_list)
	
	# (3) Pickle the pairwise group list
	file_name = 'pairwise_group_list.pkl'
	pickle_file(file_name,pairwise_group_list)
	print "file pickled"
	
    if argv == '2':
	
	group_1 = pair1[1:-1]
	group_2 = pair2[:-1]
	group_pair = [group_1,group_2]
	print group_pair
	calculate_fgoc_score(group_pair)

    if argv == '3':
	# (1) Use the pairwise list to obtain the FGOC score. The orientation is not considered.
	# Unpickle pairwise group list
	file_name = 'pairwise_group_list.pkl'
	print "Unpickling ..."
	pairwise_group_list = unpickle_file(file_name)
	
	# (2) Write into output files (4 of them)
	index = 0
	lst_dir = os.getcwd() + '/sh_files/'
	for pair in pairwise_group_list:
	    index += 1
	    if index == 5: index = 1
	    
	    file_name = lst_dir + 'fgoc_all_pairs_pkl_' + str(index) + '.sh'
	    info_line = 'python get_fgoc_pair.py ' + str(pair[0]) + ' ' + str(pair[1]) + '\n'
	    ofile = open(file_name,'a')
	    ofile.write(info_line)
	    ofile.close()


    if argv == '4':
	# Making a dictionary of organism acc_nbr: nbr_genes and genomesize
	file_name = 'orgs_accnbr_nbrgenes.txt'
	ifile = open(file_name)
	lines = ifile.readlines()
	ifile.close()
	acc_nbr_genes_dict = {}
	for line in lines:
	    [acc_nbr,genome_size,nbr_genes] = (line[:-1]).split('\t')
	    acc_nbr_genes_dict[acc_nbr] = [genome_size, nbr_genes]
	
	# Pickle dictionary
	file_name = 'acc_nbr_genes_dict.pkl'
	pickle_file(file_name,acc_nbr_genes_dict)

    if argv == '5':
	file_name = os.getcwd() + '/sh_files/' + 'fgoc_all_pairs_pkl_ko_8.sh'
	ifile = open(file_name)
	lines = ifile.readlines()
	ifile.close()
	index = 0
	flag = False
	for line in lines:
	    
	    if flag:
		index += 1
		if index == 9: index = 1
		file_name = os.getcwd() + '/sh_files/' + 'fgoc_break_8_' + str(index) + '.sh'
		ofile = open(file_name,'a')
		ofile.write(line)
		ofile.close()
		print line[:-1]
		# 'b1442', 'b2199'
	    if line.split(' ')[2] == 'b1442' and (line.split(' ')[3])[:-1] == 'b2199':
		print "Hurray Found"
		print line
		flag = True



if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - make_genevsgene_mat.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    
