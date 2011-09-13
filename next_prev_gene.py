#!usr/bin/python

# Last modified: 10th Feb 2010
# Script for retrieving a next gene neighborhood upstream and next gene downstream for all the genes in all organisms. 

import os
import sys
import pickle
import pdb
import string
from operator import itemgetter
import shelve
import fnmatch
#fav_gene = 'SGR_5936'	# This is the favorite gene that one can enter to see which context it is in
#fav_gene = 'BAG84255.1'
fav_gene = 'BAG84248.1'

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

def get_files(path_dir,pattern):
    """
    @param path_dir: This is the directory that must be searched (os.walk) to retrieve a list of all files matching the pattern
    @param pattern: This is the pattern that must be looked for in the file
    @function: This function walks through the path_dir looking for the file with the matching pattern and returning a list
    """
    file_list = []

    for root, dirs,files in os.walk(path_dir):
	for file in files:
	    if (file.find(pattern)>-1):
		file_list.append(os.path.abspath(os.path.join(root,file)))
	
    return file_list

def unpickle_file(file_name):
    """
    @function: Unpickles the file and returns the contents
    """
    ifile = open(file_name)
    contents = pickle.load(ifile)
    ifile.close()
    return contents

def pickle_file(contents, file_name):
    """
    @function: Pickles the file with the passed contents
    """
    ofile = open(file_name,'wb')
    pickle.dump(contents,ofile)
    ofile.close()

def get_adj_gene(locus_cds_dict,status):
    """
    @param locus_cds_dict: This is locus tag: cds information
    @param status : 'next' or 'prev'(previous)
    @function: Produce a dictionary of next genes. i.e next_gene[gene1] should be gene2. And so on.
    @notes : I must ensure the circularity of genome and the next gene of the last gene must be the first gene
    """
    # Initializing
    nbr_genes = len(locus_cds_dict)
    adj_gene_dict = {}
    
    all_cds = locus_cds_dict.values()
    for cds in all_cds:
	cds[1] = int(cds[1])
    sort_all_cds = sorted(all_cds,key = itemgetter(1))
    
    if status == 'next':
	for i in xrange(0,nbr_genes):
	    if not i == nbr_genes-1:
		try:
		    adj_gene = sort_all_cds[i+1]
		except:
		    print i, sort_all_cds[i],locus_cds_dict[sort_all_cds[i]]
		    sys.exit(1)
	    else:
		adj_gene = sort_all_cds[0]
	    adj_gene_dict[sort_all_cds[i][0]] = adj_gene
    
    if status == 'prev':
	for i in xrange(0,nbr_genes):
	    if not i == 0:	#['RSal33209_0002', 1, 'forward', '1864', '3090', 'NC_010168']
	    	adj_gene = sort_all_cds[i-1]
	    else:
		adj_gene = sort_all_cds[nbr_genes-1]
	    adj_gene_dict[sort_all_cds[i][0]] = adj_gene
	
    return adj_gene_dict

def get_gene_context(locus_tag, acc_nbr, acc_nbr_dir, window_limit):
    """
    @param locus_tag: Locus tag which is the blast hit of the favorite gene
    @param acc_nbr_dir: The directory in which the locus_tag is present. Get back the locus_tag of all the genes both upstrem as well as downstream
    @function: Obtain next genes (upto the window limit) for the locus_tag queried. Go both downstream (next) and upstream (up)
    """
    gene_context_list = [locus_tag]
    next_genes_list = [locus_tag]
    prev_genes_list = [locus_tag]
    next_gene_dict_name = acc_nbr_dir + '/' + acc_nbr + '_next_gene_dict.pkl'
    next_gene_dict = unpickle_file(next_gene_dict_name)
    prev_gene_dict_name = acc_nbr_dir + '/' + acc_nbr + '_prev_gene_dict.pkl'
    prev_gene_dict = unpickle_file(prev_gene_dict_name)
    
    for i in xrange(0,window_limit):
	cur_gene = next_genes_list[i]
	next_gene = next_gene_dict[cur_gene][0]
	next_genes_list.append(next_gene)
	gene_context_list.append(next_gene)
	cur_gene = prev_genes_list[i]
	prev_gene = prev_gene_dict[cur_gene][0]
	prev_genes_list.append(prev_gene)
	gene_context_list.append(prev_gene)
	
    return next_genes_list, prev_genes_list, gene_context_list
    
def get_seq_genes(acc_nbr, acc_nbr_dir, list_genes):
    """
    @param acc_nbr: This is the accession number from which a locus tag: sequence information can be obtained
    @param acc_nbr_dir: Directory containing the dictionary 
    @param list_genes: A list of all genes whose sequence must be obtained
    @function: Run through the list of genes and retrieve the locus tag and sequence. 
    """
    
    locus_seq_dict_name = acc_nbr_dir + '/' + acc_nbr + '_locus_sequence_dictionary'
    locus_seq_dict = unpickle_file(locus_seq_dict_name)
    gene_seq_list = []
    
    for gene in list_genes:
	locus_tag, seq = locus_seq_dict[gene]
	gene_seq_list.append([acc_nbr + '|' + locus_tag, seq])

    return gene_seq_list
    

def get_neighbourhood(blast_hits):
    """
    @param blast_hits: These are the acc_nbr|locus_tag of all the genes that were the blast hits for the favorite gene
    @function: Obtain a 10 gene neighborhood for each of the genes and give back a list of all the genes along with their acc_nbr
    """
    # Initializing
    window_limit = 5	#This is the number of genes one must get 5 genes upstream and 5 genes downstream
    path_dir = os.getcwd() + '/parsed_gbk_files_genomes'
    acc_nbrandlocus_context_dict = {}
    all_genes = []
    # The file to store all the gene sequences. The identifier is the acc_nbr|locus_tag 
    fav_gene_all_seq_name = os.getcwd() + '/fav_genes/' + fav_gene + '/' + fav_gene + '_all_gene_context_seq.fasta'
    # Iterate through the gene_seq_list and append the locus tag and sequence to the all-seq-file
    ofile = open(fav_gene_all_seq_name,'w')

    for hit in blast_hits:
	# Obtain acc_nbr to track the next and prev genes dict 
	acc_nbr, locus_tag = hit.split('|')
	acc_nbr_dir = path_dir + '/' + acc_nbr + '/'
	next_genes_list, prev_genes_list, gene_context_list = get_gene_context(locus_tag, acc_nbr, acc_nbr_dir, window_limit)
	# Make a acc_nbr|locus_tag : [locus_tag, next_genes, prev_genes, set of these 11 genes surrounding the locus_tag
	acc_nbrandlocus_context_dict[hit] = [acc_nbr, locus_tag, next_genes_list,prev_genes_list,gene_context_list]
	# Also obtain sequences from all the locus_tag in the gene_context_list. Append to a common file 
	gene_seq_list = get_seq_genes(acc_nbr, acc_nbr_dir, gene_context_list)
	for locus_seq in gene_seq_list:
	    locus_tag, seq = locus_seq
	    ofile.write('>' + locus_tag + '\n' + seq + '\n')
    ofile.close()
    
    # Pickle the dictionary
    pkl_file = os.getcwd() + '/fav_genes/' + fav_gene + '/' + 'locus_id_gene_context_dict.pkl'
    pickle_file(acc_nbrandlocus_context_dict,pkl_file)
    print pkl_file
    
    return 

def get_locus_gp_dict(gp_file):
    """
    @param gp_file: The path to the groups.txt file - generated after orthomcl clustering. 
    @function: Parse the file so that each locus_tag gets assigned a group_number
    """
    locus_gp_dict = {}

    ifile = open(gp_file)
    lines = ifile.readlines()
    ifile.close()
    
    for line in lines:
	line = line.strip('\n')
	gp_nbr, all_locus = line.split(':')
	acc_locus_list = all_locus.split(' ')
	locus_list = [item.split('|')[1] for item in acc_locus_list if not (item == '')]
	for locus in locus_list:
	    locus_gp_dict[locus] = gp_nbr

    return locus_gp_dict

def order_list(next_list, prev_list,locus_tag):
    """
    @function: There is overlapping locus_tag between the two lists. 
    """
    ordered_list = []
    prev_list.reverse()
    for item in prev_list:
	ordered_list.append(item)
    for item in next_list:
	if not item == locus_tag:
	    ordered_list.append(item)
    return ordered_list


def main(argument):
    [argv] = argument
  
    
    if argv == '0':
	
	# Initializing
	base_dir = os.getcwd()
	path_dir = base_dir + '/parsed_gbk_files_genomes'
	locus_cds_pkl = []

	# (1) Obtain the list of all files - locus:cds_information dictionary
	pattern = 'locus_cds_information_dictionary.pkl'
	locus_cds_pkl_files = get_files(path_dir,pattern)
	
	for file in locus_cds_pkl_files:
	    locus_cds_dict = {}
	    dir_name = os.path.dirname(file)
	    org_name = dir_name.split('/')[-1]
	    print " Processing %s ..."%(org_name)
	    # (2) Unpickle file for locus:cds information dictionary
	    locus_cds_dict = unpickle_file(file)
	    # (3) Obtain a dictionary of next gene
	    next_gene_dict = get_adj_gene(locus_cds_dict,'next')
	    # (4) Obtain a dictionary for the previous gene
	    prev_gene_dict = get_adj_gene(locus_cds_dict,'prev')
	    # (5) Pickle the two dictionaries in the folder containing all other details
	    pkl_file = dir_name + '/' + org_name + '_next_gene_dict.pkl'
	    pickle_file(next_gene_dict,pkl_file)
	    pkl_file = dir_name + '/' + org_name + '_prev_gene_dict.pkl'
	    pickle_file(prev_gene_dict,pkl_file)
	    

    if argv == '1':
	# This part will take in the query protein - locus_tag. It will then search the shelve file to obtain information and then go to the particular organism locus tag: sequence information to obtain the sequence and make a file. Then write a script which will blast against all genomes and output to a common blast file.
	# Initializing
	only_seq_there = True
	base_dir = os.getcwd()
	path_dir = base_dir + '/parsed_gbk_files_genomes'
	
	if not (only_seq_there):
	    pkl_dir = os.getcwd() + '/shelve_files/'
	    locus_cds_information_dictionary = shelve.open(pkl_dir + 'locus:cds_information_dictionary','c')

	    
	    # (1) Obtain the information of the favorite gene - especially the org_nbr
	    org_nbr = locus_cds_information_dictionary[fav_gene][5]
	    # (2) Based on the org_nbr unpickle the locus: sequence file
	    seq_file = path_dir + '/' + org_nbr + '/' + org_nbr + '_locus_sequence_dictionary'
	    locus_seq_dict = unpickle_file(seq_file)
	    fav_gene_seq = locus_seq_dict[fav_gene][1]
	    
	    # (3) Make a fasta file of this gene which will be the query
	    dir_path = os.getcwd() + '/fav_genes'
	    fav_gene_dir_path = make_dir(fav_gene,dir_path)
	    query_file_name = fav_gene_dir_path + '/' + fav_gene + '_query.fasta'
	    ofile = open(query_file_name,'w')
	    ofile.write('>' + fav_gene + '\n' + fav_gene_seq + '\n')
	    ofile.close()
	
	else:
	    dir_path = os.getcwd() + '/fav_genes'
	    fav_gene_dir_path = make_dir(fav_gene,dir_path)
	    query_file_name = fav_gene_dir_path + '/' + fav_gene + '_query.fasta'

	# (4) Obtain a list of all .fasta files 
	matches = []
	for root, dirs,files in os.walk(path_dir):
	    for file in fnmatch.filter(files,'*.fasta'):
		matches.append(os.path.join(root,file))
	# (5) Make a sh file with path of query file blasted against the list of files in the matches. 
	sh_file = fav_gene_dir_path + '/' + fav_gene + '_blast_all.sh'
 	evalue = '0.1'
	ofile = open(sh_file,'w')
	for db_file in matches:
	    command_line = 'blastall -p blastp -d ' + db_file + ' -i ' + query_file_name  + ' -e ' + evalue + ' -m 8 ' + '\n' 
	    # Will have to run the command from shell and redirect output. The result file cannot be appended
	    ofile.write(command_line)

    if argv == '2':
	# Parse the blast output which is generally fav_gene name _blast_all.result
	# Initializing
	dir_path = os.getcwd() + '/fav_genes'
	fav_gene_dir_path = dir_path + '/' + fav_gene 
	
	# (1) Parse the blast result and get a list of locus_tag and the acc_nbr for blast hits that are above a particular identity threshold
	blast_result = fav_gene_dir_path + '/' + fav_gene + '_blast_all.result'
	ifile = open(blast_result)
	lines = ifile.readlines()
	ifile.close()
	identity_cutoff = 10
	blast_hits = [line.split('\t')[1] for line in lines if float(line.split('\t')[2]) > identity_cutoff]
	# (2) Get information of the next 5 and previous 5 genes for each of the blast hit
	get_neighbourhood(blast_hits)
	
    if argv == '3':
	pkl_dir = os.getcwd() + '/shelve_files/'
	locus_cds_information_dictionary = shelve.open(pkl_dir + 'locus:cds_information_dictionary','c')
	diff_context = {}
	
	# This part analyses the groups.txt result (after clustering the genes using orthomcl) and retrieves all possible grouping pattern this gene was in.
	# (1) Parse groups.txt and make a dictonary of locus_tag: group_number
	gp_file = os.getcwd() + '/fav_genes/' + fav_gene + '/groups.txt'
	locus_gp_dict = get_locus_gp_dict(gp_file)
	# (2) Use the locus_id_gene_context_dict to retrieve data about the different genes. 
	#acc_nbrandlocus_context_dict[hit] = [acc_nbr, locus_tag, next_genes_list,prev_genes_list,gene_context_list]
	pkl_file = os.getcwd() + '/fav_genes/' + fav_gene + '/' + 'locus_id_gene_context_dict.pkl'
	locus_id_gene_context_dict = unpickle_file(pkl_file)
	for key, value in locus_id_gene_context_dict.items():
	    locus_id = key
	    gp_order = []
	    gp_order_locus = []
	    acc_nbr, locus_tag, next_genes_list, prev_genes_list = value[0], value[1], value[2], value[3]
	    ordered_list = order_list(next_genes_list, prev_genes_list,locus_tag)
	    for item in ordered_list:
		try:
		    gp_order_locus.append([item, locus_gp_dict[item]])
		    gp_order.append(locus_gp_dict[item])
		except KeyError:
		    gp_order_locus.append([item, 'None'])
		    gp_order.append('None')
	    
	    if gp_order.count('None') < 8 and gp_order.count('JEM_1'):
		gp_order_seq = '->'.join(item for item in gp_order)
		diff_context[locus_id] = [ordered_list,gp_order,gp_order_seq, gp_order_locus]
		
	all_group_order_seq = []
	for key, value in diff_context.items():
	    locus_order, group_order, group_order_seq, group_order_locus = value
	    print "%s : %s : %s"%(key,group_order_seq, group_order_locus)
	    all_group_order_seq.append(group_order_seq)
	
	for item in set(all_group_order_seq):
	    print item, all_group_order_seq.count(item)

    if argv == '4':
	pkl_dir = os.getcwd() + '/shelve_files/'
	locus_cds_information_dictionary = shelve.open(pkl_dir + 'locus:cds_information_dictionary','c')
	
	print locus_cds_information_dictionary['BAG84255']

if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - next_prev_gene.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    


