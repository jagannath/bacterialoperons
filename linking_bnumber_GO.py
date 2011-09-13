#! usr/bin/python

# Understanding genes/operons involved with respect to GO terminologies. Have downloaded the GO (bp=Biological process, mf = Molecular Function and cc = Cellular Component). Will join tables and pickle with bnumber : [[bp, mf, cc], operon]. Sometimes the operon name will not be available as so defaulted.
# Downloaded from bacteria Ensembl information about the GO information (in download directory - biomart/ecoli_GO_info.csv and there they dont have a bnumber. So have to use ecogenedata.txt to retrieve information of gene_name: bnumber. 
from __future__ import division
import os
import sys
import pickle
from operator import itemgetter
import pdb

from math import factorial


# Code copied from http://stackoverflow.com/questions/2096573/counting-combinations-and-permutations-efficiently
def product(iterable):
    prod = 1
    for n in iterable:
        prod *= n
    return prod

def npr(n, r):
    """
    Calculate the number of ordered permutations of r items taken from a
    population of size n.

    >>> npr(3, 2)
    6
    >>> npr(100, 20)
    1303995018204712451095685346159820800000
    """
    assert 0 <= r <= n
    return product(range(n - r + 1, n + 1))

def ncr(n, r):
    """
    Calculate the number of unordered combinations of r items taken from a
    population of size n.

    >>> ncr(3, 2)
    3
    >>> ncr(100, 20)
    535983370403809682970
    >>> ncr(100000, 1000) == ncr(100000, 99000)
    True
    """
    assert 0 <= r <= n
    if r > n // 2:
        r = n - r
    return npr(n, r) // factorial(r)


def get_go_terms(bnbr_list, bnbr_go_dict):
    """
    @param bnbr_list: This is a list of bnumbers
    @param bnbr_go_dict: This is a bnumbers:GO terms - bnumber: [bp,mf,cc]
    @function: Iterates through the list and sums up the different bp, mf, cc and gives back the highest occuring terms for each as well as the list of all
    """
    # Initializing
    bp_list = []
    mf_list = []
    cc_list = []
    
    for group_nbr in eval(bnbr_list):
	bnbr = group_nbr[6:]
	try: 
	    bp, mf, cc = bnbr_go_dict[bnbr]
	except KeyError:
	    bp = mf = cc = 'Not Annotated'
	bp_list.append(bp)
	mf_list.append(mf)
	cc_list.append(cc)

    return set(bp_list), set(mf_list), set(cc_list)

def go_enrich(top_go,total_genes_set, go_count_dict):
    """
    @param top_go: This is the top GO enrichments in the list in format (GO term, counts)
    @function: Calculate the hyper geometric distribution; Think this is correct
    """
    total_go_terms = 4323
    go_enrich_list = []
    
    for one_go, count in top_go:
	try:
	    total_go = go_count_dict[one_go]
	except KeyError:
	    total_go = go_count_dict['']
	observed_go = count
	# Total go terms = 3926 (N); Total number of genes in set (i.e number of draws) = n; Number of GO terms observed in the set (draw) = k; Total number of the particular GO term in all GO terms = m; I am not sure if it is correct but I am calculating some probability value by using the hyper geometric test; Formula from wikipedia - hypergeometric
	N = total_go_terms
	n = total_genes_set	#The total number of genes observed in this set 
	m = total_go	# The total number of times this go term was observed in the all the genes - big GO term set 
	k = observed_go	# Number of occurrences of this GO term in this set
	prob = (ncr(m,k) * ncr((N-m), (n-k)))/ncr(N,n)
	go_enrich_list.append([one_go,count,prob,(N,n,m,k)])

    return go_enrich_list

def main(argument):
    
    [argv] = argument
    
    if argv == '1':
	# Initializing
	name_bnumber_dict = {}
	bnbr_go_dict = {}
	all_bp_list = []
	all_mf_list = []
	all_cc_list = []
	bp_count_dict = {}
	mf_count_dict = {}
	cc_count_dict = {}

	# (1) Open the ecogenedata.txt file to retrieve gene_name: bnumber information
	ifile = open('ecogenedata.txt','r')
	lines = ifile.readlines()
	ifile.close()
	name_bnumber_tuple = [(line.split('\t')[1], line.split('\t')[3]) for line in lines]

	# (2) Make a dictionary of gene name : bnumber (locus_tag for ecoli)
	for name,bnumber in name_bnumber_tuple:
	    name_bnumber_dict[name] = bnumber

	# (3) Open the GO information file and read lines
	file_name = '/project/marcotte/jagannath/projectfiles/Downloads/Biomart/ecoli_GO_info.csv'
	ifile = open(file_name)
	lines = ifile.readlines()
	ifile.close()

	# (4) Iterate through lines and ensure no duplicity first line is Description
	for line in lines:
	    if not line.startswith('Gene'):
		gene_name, bp, mf, cc = line.split(',')[1], line.split(',')[2], line.split(',')[3], line.split(',')[4][:-1]
		# (5) Make dictionary of bnumber : [bp, mf, cc] by calling the name_bnumber_dict
		try: 
		    bnbr = name_bnumber_dict[gene_name]
		    bnbr_go_dict[bnbr] = [bp, mf, cc]
		# (6) Append the bp, mf, cc to a list
		except KeyError:
		    pass
	
	all_bp_list = [item[0] for item in bnbr_go_dict.values()]
	all_mf_list = [item[1] for item in bnbr_go_dict.values()]
	all_cc_list = [item[2] for item in bnbr_go_dict.values()]
	
	# (7) Make a pickle of bnumber : GO terms
	pkl_file = os.getcwd() + '/pkl_files/bnumber:GOterms'
	ofile = open(pkl_file,'wb')
	pickle.dump(bnbr_go_dict, ofile)
	ofile.close()
	print "bnbr_go_dict pickled" 
	print len(bnbr_go_dict)
	print len(name_bnumber_dict)
	
	# (8) Make a shortened list containing all the bp etc
	bp_count_list = [(item, all_bp_list.count(item)) for item in set(all_bp_list)]
	mf_count_list = [(item, all_mf_list.count(item)) for item in set(all_mf_list)]
	cc_count_list = [(item, all_cc_list.count(item)) for item in set(all_cc_list)]
	# (9) Make a dictionary; For '' call it Not Annotated
	for key, value in bp_count_list: bp_count_dict[key] = value
	for key, value in mf_count_list: mf_count_dict[key] = value
	for key, value in cc_count_list: cc_count_dict[key] = value
	
	# (10) Pickle the files
	pkl_file = os.getcwd() + '/pkl_files/ecoli_bp:counts'
	ofile = open(pkl_file,'wb')
	pickle.dump(bp_count_dict, ofile)
	ofile.close()
	pkl_file = os.getcwd() + '/pkl_files/ecoli_mf:counts'
	ofile = open(pkl_file,'wb')
	pickle.dump(mf_count_dict, ofile)
	ofile.close()
	pkl_file = os.getcwd() + '/pkl_files/ecoli_cc:counts'
	ofile = open(pkl_file,'wb')
	pickle.dump(cc_count_dict, ofile)
	ofile.close()

	print "Files pickled"
	print bnbr_go_dict['b3252']
	print bp_count_dict['oxidation reduction']

    if argv == '2':
	# This analyses the operons that did not come back as rank 1. There were like 45 operons in E.coli. Will get a list of all GO terms included in the list, sort them in decreasing order and see if it makes any sense
	# File arranged as rank \t [bnumber1, bnumber2, ..] \t operon
	# Initializing
	all_bp_list = []
	all_mf_list = []
	all_cc_list = []
	total_go_terms = 3926
	total_gene_pairs = 0
	
	# (0) Unpickle the bnbr_go_dict
	pkl_file = os.getcwd() + '/pkl_files/bnumber:GOterms'
	ifile = open(pkl_file)
	bnbr_go_dict = pickle.load(ifile)
	ifile.close()
	print "bnbr_go_dict unpickled" 
	
	pkl_file = os.getcwd() + '/pkl_files/ecoli_bp:counts'
	ifile = open(pkl_file)
	bp_count_dict = pickle.load(ifile)
	ifile.close()

	pkl_file = os.getcwd() + '/pkl_files/ecoli_mf:counts'
	ifile = open(pkl_file)
	mf_count_dict = pickle.load(ifile)
	ifile.close()
	
	pkl_file = os.getcwd() + '/pkl_files/ecoli_cc:counts'
	ifile = open(pkl_file)
	cc_count_dict = pickle.load(ifile)
	ifile.close()
	

	# (1) Open and readlines in file
	file_name = 'ecoli_homolog_greater_rankers.txt'
	ifile = open(file_name,'r')
	lines = ifile.readlines()
	ifile.close()

	for line in lines:
	    rank, bnbr_list, operon = int(line.split('\t')[0]), line.split('\t')[1], line.split('\t')[2]

	    total_gene_pairs += len(eval(bnbr_list))
	    set_bp_list, set_mf_list, set_cc_list = get_go_terms(bnbr_list, bnbr_go_dict)
	    # Need to flatten the set
	    bp_list = [item for item in set_bp_list]
	    mf_list = [item for item in set_mf_list]
	    cc_list = [item for item in set_cc_list]
	    for item in bp_list: all_bp_list.append(item)
	    for item in mf_list: all_mf_list.append(item)
	    for item in cc_list: all_cc_list.append(item)

	# (2) Sorting the list after obtaining the counts
	cc_counts = [(item, all_cc_list.count(item)) for item in set(all_cc_list) if (item)]
	bp_counts = [(item, all_bp_list.count(item)) for item in set(all_bp_list) if (item)]
	mf_counts = [(item, all_mf_list.count(item)) for item in set(all_mf_list) if (item)]
	sorted_cc_counts = sorted(cc_counts, key = itemgetter(1), reverse=True)
	sorted_bp_counts = sorted(bp_counts, key = itemgetter(1), reverse=True)
	sorted_mf_counts = sorted(mf_counts, key = itemgetter(1), reverse=True)
	
	top_bp = sorted_bp_counts[0:10]
	top_mf = sorted_mf_counts[0:10]
	top_cc = sorted_cc_counts[0:10]
	total_genes_set = total_gene_pairs * 2
    
	bp_enrich_list = go_enrich(top_bp,total_genes_set, bp_count_dict)
	mf_enrich_list = go_enrich(top_mf,total_genes_set, mf_count_dict)
	cc_enrich_list = go_enrich(top_cc,total_genes_set, cc_count_dict)
	print bp_enrich_list
	print mf_enrich_list
	print cc_enrich_list
	
	
    if argv == '3':
	# This part is for obtaining GO terms for gene pairs within operons that were conserved as only < 10% of the genome fgoc score across all species. Seeing if there is any type of enrichment
	# Utilizes the file - within_orgs_ref_ecoli_ver6.genome_fgoc where the org_acc_nbr[0], norm_fgoc[4], gene_pair list[5]
	
	pkl_file = os.getcwd() + '/pkl_files/bnumber:GOterms'
	ifile = open(pkl_file)
	bnbr_go_dict = pickle.load(ifile)
	ifile.close()
	print "bnbr_go_dict unpickled" 
	
	# Initializing
	all_bp_list = []
	all_mf_list = []
	all_cc_list = []
	all_sp_gene_pairs = []
	
	# (1) Opening the file containing the information of genome fgoc
	ifile = open('within_orgs_ref_ecoli_ver6.genome_fgoc','r')
	lines = ifile.readlines()
	ifile.close()
	
	# (2) Filter the orgs to obtain only those with a cutoff of norm fgoc score of less than 10% on the norm_fgoc score
	cutoff = 5
	count = 0
	
	for line in lines:
	    if float(line.split('\t')[4]) < cutoff:
		count += 1
		acc_nbr, norm_fgoc, gene_pair_list = line.split('\t')[0], float(line.split('\t')[4]), eval(line.split('\t')[5])
		print len(gene_pair_list), norm_fgoc
		# (3) Iterate through the list to obtain all the gene pairs and get the set of GO terms 
		for pair in gene_pair_list:
		    # (4) Make a list of gene pairs from all the species within this cut off. 
		    all_sp_gene_pairs.append(pair)
	print len(set(all_sp_gene_pairs))
		    
	# (5) Make a non redundant set of all gene pairs
	nr_all_spp_gene_pairs = set(all_sp_gene_pairs)
	for group_pair in nr_all_spp_gene_pairs:
	    set_bp_list, set_mf_list, set_cc_list = get_go_terms(group_pair, bnbr_go_dict)
	    # Need to flatten the set
	    bp_list = [item for item in set_bp_list]
	    mf_list = [item for item in set_mf_list]
	    cc_list = [item for item in set_cc_list]
	    for item in bp_list: all_bp_list.append(item)
	    for item in mf_list: all_mf_list.append(item)
	    for item in cc_list: all_cc_list.append(item)

	# (4) Sorting the list after obtaining the counts
	cc_counts = [(item, all_cc_list.count(item)) for item in set(all_cc_list) if (item)]
	bp_counts = [(item, all_bp_list.count(item)) for item in set(all_bp_list) if (item)]
	mf_counts = [(item, all_mf_list.count(item)) for item in set(all_mf_list) if (item)]
	sorted_cc_counts = sorted(cc_counts, key = itemgetter(1), reverse=True)
	sorted_bp_counts = sorted(bp_counts, key = itemgetter(1), reverse=True)
	sorted_mf_counts = sorted(mf_counts, key = itemgetter(1), reverse=True)
	
	print count
	print sorted_bp_counts[0:10]
	print sorted_mf_counts[0:10]
	print sorted_cc_counts[0:10]




if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    
    print "Processing %s ..."%(argument)
    main(argument)    
    
    
    
    import time
    print "Script - linking_bnumber_GO.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    

