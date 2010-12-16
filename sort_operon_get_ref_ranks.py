#******************************************************************************************* #
#I had to block the above section because this is a modification to the programme. The script will now be able to open all the .shuffled file in the directory shuffled_operon, parse and sort the list to and write operon_name.sort It will also retrieve what the rank of the gene order (ref ecoli gene order) is among all the various possibilities

from __future__ import division
import os
import sys
from operator import itemgetter

def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """
    
    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	sys.exit(0)
    
    return file

def parse_one_shuffle(one_group):
    """
    @param one_group:	This is string passed from parse_shuffled_operon which split the gene order into one shuffle. 
    @function:		Parses this information of one shuffled operon and obtains (a) Order (b) Total conserved gene pair and (c) #conserved gene_pair(occurrence of gene1,gene2):#conserved gene_pair(..) 
    @param output:	Returns back a line of this information, each column tab delimited
    """
    # Initializing
    total_conserved = 0
    operon_order = ''
    conserved_gene_pair = ''
    conserved_gene = ''
    one_conserved_gene_pair = ''
    gene_1_occ = gene_2_occ = ''
    frac_conserved_gene_pair = 0
    total_frac = 0
    split_group = one_group.split('\n')
    
    
    
    for item in split_group:
	if item.startswith('Group Order'):
	    operon_order = item.split(':')[1]
	if item.startswith('# Conserved gene pair'):
	    one_conserved_gene_pair = int((item.split('\t'))[1])
	if item.startswith('# Occurrences of gene 1'):
	    gene_1_occ = (item.split('\t'))[1]
	if item.startswith('# Occurrences of gene 2'):
	    gene_2_occ = (item.split('\t'))[1]
	if item.startswith('# Non adjacent gene pairs'):
	    non_adjacent_gene_pair = int((item.split('\t'))[1])
	    try:
		frac_conserved_gene_pair = one_conserved_gene_pair / (one_conserved_gene_pair + non_adjacent_gene_pair)
	    except ZeroDivisionError: 
		print one_conserved_gene_pair, non_adjacent_gene_pair
	    conserved_gene += str(one_conserved_gene_pair) + '(' + gene_1_occ + ',' + gene_2_occ + '):' + str(frac_conserved_gene_pair) + '||'
	    total_frac += frac_conserved_gene_pair
	if item.startswith('Total conserved gene_pairs'):
	    total_conserved = item.split('=')[1]

    number_gene_pairs = conserved_gene.count('||')
    norm_frac = total_frac/ number_gene_pairs
    
    return [operon_order, norm_frac, conserved_gene]	
    
    
def parse_shuffled_operon(file_name):
    """
    @param file_name:	This is the file_name containing the output information of the operon_name.shuffled. Each of the operon shuffle is delimited by '///'
    @function:		Parses the file generating a list containing the following information, each column separated by a \t
    @param output:	Returns the created list
    """

    #Initializing
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    
    one_group = ''
    one_group_list = []
    all_pair_frac_list = []
    operon_count = 0
    
    flag = False
    
    for line in lines:
	if line.startswith('Group Order'):
	    flag = True
	if line.startswith(' ///'):
	    flag = False
	    one_group_list = parse_one_shuffle(one_group)
	    all_pair_frac_list.append(one_group_list)
	    one_group = ''
	if flag:
	    one_group += line
    
		    
    return all_pair_frac_list


def write_sorted_list(file_name,sorted_list,cut_off=0):
    """
    @param file_name:	The file name for writing the output of the calculated and sorted list
    @param sorted_list:	The list contains three items - the first can be the operon_pair, group_pair or operon_order: the second is always the fraction : the third is the gene_occurrences
    @function:		Writes a file with the output of the sorted list. It also numbers the first column as the serial_nbr
    @output:		None
    """
    
    #Initializing
    num_list = []
    above_cut_off_group = []
    number_gene_pairs = 0
    serial_nbr = 0

    ofile = open_file(file_name,'w')
    for one_row in sorted_list:
	serial_nbr += 1
	num_list.append(one_row[1])
	group_order = '[' + (str(one_row[0]))[2:-2] + ']'

	row = str(serial_nbr) + '\t' + str(group_order) + '\t' + str(one_row[1]) + '\t' + one_row[2] + '\n'
   
	ofile.write(row)
	if one_row[1] > cut_off:
	    above_cut_off_group.append([one_row[0], one_row[1], one_row[2]])
		
    number_gene_pairs = len(num_list)
    #print "Number_gene_pairs : %d"%(number_gene_pairs)
    #print "Number_gene_pairs above %f cut off of fraction conserved gene pair: %d  or %d percentage "%(cut_off, len(above_cut_off_group), (((len(above_cut_off_group))/number_gene_pairs)*100))
    #print "***"
    ofile.close()
    
    
    return num_list

def get_list_shuffled_files():
    """
    @function:	This function obtains the path of all the .shuffled files in the directory - shuffled_operons
    @param path_list:	Returns back the path of all the files ending with .shuffled 
    """
    
    dir_path = os.getcwd() + '/shuffled_operons'
    path_list = [dir_path + '/' + file_name for file_name in os.listdir(dir_path) if file_name.endswith('shuffled')]
    
    
    return path_list

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

def get_operon_names():
    """
    @function:	This function searches the directory shuffled_operons and from that obtains only the written operon file names (with a .shuffled extension)
    @param output:	Returns the list of all the operon names
    """
    
    dir_path = os.getcwd() + '/shuffled_operons'
    operon_names = [file_name[:-9] for file_name in os.listdir(dir_path) if file_name.endswith('shuffled')]
    
    return operon_names


def check_rank(operon_name, group_order):
    """
    @param operon_name:	This is the name of the operon for which the rank of the group order is obtained
    @param group_order:	This is the actual group order of the operon
    @function:		The function first opens the file operon_name.shuffled_sort. Converts the group order to a string and gets the serial_nbr where the group order which was passed matches with the 2nd column (the group_order column) in the file.
    @param output:	The serial number that matched
    """
    
    #Initializing
    ref_group_order = str(group_order)
    file_name = os.getcwd() + '/shuffled_operons/' + operon_name + '.shuffled_sort'
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    
    for line in lines:
	if ref_group_order == line.split('\t')[1]:
	    rank = line.split('\t')[0]
	    found_group_order = line.split('\t')[1]
	    return rank
    
    
def get_incomplete_operons():
    """
    @function:	This function obtains operons which are incomplete - i.e. one of its gene isnt a member of any cluster
    @return:	Number of such operons
    """
    ifile = open_file('ecoli_cluster.txt')
    lines = ifile.readlines()
    ifile.close()
    all_incomplete_operons = []
    
    row_0 = lines[0].split('\t')
    i = 1
    operon_0, rank_0 = row_0[0], row_0[3]
    for row in lines:
	status = 'complete'
	if i>1:
	    operon_next, rank_next = row.split('\t')[0],row.split('\t')[3]
	    
	    if operon_0 == operon_next:
		rank_diff = int(rank_next) - int(rank_0)
		if rank_diff == 1:
		    pass
		else:
#		    print operon_0, rank_0, operon_next, rank_next
		    incomplete_operon = operon_0
		    all_incomplete_operons.append(incomplete_operon)
	    else:
		pass
	
	    operon_0 = operon_next
	    rank_0 = rank_next
	i+=1

    nbr_incomplete = len(set(all_incomplete_operons))
    
    return nbr_incomplete
    


if __name__ == '__main__':
    
    all_files_paths = get_list_shuffled_files()		#Obtains a path list of all the files 
    #for file_path in all_files_paths:
	#all_shuffles_list = parse_shuffled_operon(file_path)
	#sorted_list_operon = sorted(all_shuffles_list, key = itemgetter(1), reverse=True)
	#new_file_name = file_path + '_sort'
	#num_list = write_sorted_list(new_file_name,sorted_list_operon)

    #cluster_dictionary = get_cluster_locus_tag_dictionary()
    operon_list = get_operon_names()
    list_ranks = []
    
    for operon in operon_list:
	ecoli_genes_information, group_order, orientation_status = get_ecoli_genes(operon)
	if orientation_status == 'reverse':
	    group_order.reverse()
	get_rank = check_rank(operon, group_order)
	if get_rank != '1':
	    print operon, group_order
	list_ranks.append(get_rank)

    nbr_incomplete_operons = get_incomplete_operons()
