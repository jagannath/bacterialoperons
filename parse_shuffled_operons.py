#!/usr/bin/env python

# Script parses the output of the shuffled tryptophan file or any of the operon that was shuffled. It obtains the total number of conserved gene pairs for a given operon order. It then sorts the list. 
# I dont want to use probability score as there may be one gene pair with 0 conserved gene pair which will make the total score fall to 0. I could add the individual probability and normalize it based on the number of walks.

# Mainly parses the file - possible_trp.out (in this case) and obtains ordering information and total conserved gene pair information for that ordering
# Output will be in columns (1) Ordered group list (2) Total conserved pairs (3) Number conserved gene pair across each walk delimited by ':'
from __future__ import division
import numpy as np
import pylab as p
import matplotlib
import os
import sys
from operator import itemgetter
import csv


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
    @param one_group:	This is string passed from parse_shuffled_trp which split the file into one shuffle. 
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
    
    
def parse_shuffled_trp(file_name,operate_on):
    """
    @param file_name:	This is the file_name containing the output information of the shuffled_operon. Each of the operon shuffle is delimited by '///'
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
    
    if operate_on == 'operon':
	for line in lines:
	    if line.startswith('Group Order'):
		flag = True
	    if line.startswith('///'):
		flag = False
		one_group_list = parse_one_shuffle(one_group)
		all_pair_frac_list.append(one_group_list)
		one_group = ''
	    if flag:
		one_group += line
    
    # This is for the all_ecoli_genome_walk in an operon
    if operate_on == 'organism':
	
	for line in lines:
	    if line.startswith('Group Order'):
		operon_order = line.split(':')[1]
		operon_count += 1
		operon_order = eval(operon_order)	#This converts the string list to a regular list
	    if line.startswith('Operon'):
		walk_number = int(line.split(':')[1])
		group_pair = [operon_order[walk_number-1],operon_order[walk_number]]
	    if line.startswith('# Conserved gene pair'):
		one_conserved_gene_pair = int((line.split('\t'))[1])
	    if line.startswith('# Occurrences of gene 1'):
		gene_1_occ = (line.split('\t'))[1][:-1]	#There is a \n character
	    if line.startswith('# Occurrences of gene 2'):
		gene_2_occ = (line.split('\t'))[1][:-1]
	    if line.startswith('# Non adjacent gene pairs'):
		non_adjacent_gene_pair = int((line.split('\t'))[1])
		try:
		    frac_conserved_gene_pair = one_conserved_gene_pair / (one_conserved_gene_pair + non_adjacent_gene_pair)
		except ZeroDivisionError: 
		    print one_conserved_gene_pair, non_adjacent_gene_pair
	    
		all_pair_frac_list.append([group_pair, frac_conserved_gene_pair, '(' + gene_1_occ + ',' + gene_2_occ + ')'])

	print "Number of operons : %d"%(operon_count)

    # Parsing the file for genome walk between operon 
    if operate_on == 'organism_between_operon':
	for line in lines:
	    if line.startswith('Operon'):
		operon_pair, walk_number = line.split(':')[0],int(line.split(':')[1])
	    if line.startswith('# Conserved gene pair'):
		one_conserved_gene_pair = int((line.split('\t'))[1])
	    if line.startswith('# Occurrences of gene 1'):
		gene_1_occ = (line.split('\t'))[1][:-1]	#There is a \n character
	    if line.startswith('# Occurrences of gene 2'):
		gene_2_occ = (line.split('\t'))[1][:-1]
	    if line.startswith('# Non adjacent gene pairs'):
		non_adjacent_gene_pair = int((line.split('\t'))[1])
		try:
		    frac_conserved_gene_pair = one_conserved_gene_pair / (one_conserved_gene_pair + non_adjacent_gene_pair)
		except ZeroDivisionError: 
		    print one_conserved_gene_pair, non_adjacent_gene_pair
		    
		all_pair_frac_list.append([operon_pair, frac_conserved_gene_pair, '(' + gene_1_occ + ',' + gene_2_occ + ')'])
		    
    return all_pair_frac_list

def draw_histogram(number_list):
    """ 
    @param number_list:	This is the list of all numbers and have to make a histogram from this
    @function:		Makes a histogram based on the list of numbers. 
    @param output:	None
    """
    
    bins = 500
    p.figure()
    n,bins,patches = p.hist(number_list,bins,facecolor='g')
    p.xlabel('fraction of conserved gene_pair to the total number of combined gene occurrences')
    p.ylabel('frequency')
    p.title('Histogram of conservation of gene pairs within Ecoli operons')

    
    return
    
def write_sorted_list(file_name,sorted_list,cut_off=0,operate_on = 'organism'):
    """
    @param file_name:	The file name for writing the output of the calculated and sorted list
    @param sorted_list:	The list contains three items - the first can be the operon_pair, group_pair or operon_order: the second is always the fraction : the third is the gene_occurrences
    @function:		Writes a file with the output of the sorted list. It also makes a histogram of the conserved_gene_pair vs frequency of observation
    @output:		None
    """
    #Initializing
    num_list = []
    above_cut_off_group = []
    number_gene_pairs = 0

    
    if operate_on == 'organism':
	ofile = open_file(file_name,'w')
	for one_row in sorted_list:
	    row = str(one_row[0]) + '\t' + str(one_row[1]) + '\t' + one_row[2] + '\n'
	    ofile.write(row)
	    if one_row[1] > cut_off:
		above_cut_off_group.append([one_row[0], one_row[1], one_row[2]])
		num_list.append(one_row[1])

		
    if operate_on == 'operon':
	ofile = open_file(file_name,'w')
	for one_row in sorted_list:
	    #conserved_genes = one_row[2].split('||')[:-1]
	    #nbr_conserved_gene_pairs_list = [int(item.split('(')[0]) for item in conserved_genes]
	    #if 0 in nbr_conserved_gene_pairs_list:
		#pass
	    #else:
	    num_list.append(one_row[1])
	    row = str(one_row[0]) + '\t' + str(one_row[1]) + '\t' + one_row[2] + '\n'
	    
	    ofile.write(row)
	    if one_row[1] > cut_off:
		above_cut_off_group.append([one_row[0], one_row[1], one_row[2]])
		
	    
    number_gene_pairs = len(num_list)
    print "Number_gene_pairs : %d"%(number_gene_pairs)
    print "Number_gene_pairs above %f cut off of fraction conserved gene pair: %d  or %d percentage "%(cut_off, len(above_cut_off_group), (((len(above_cut_off_group))/number_gene_pairs)*100))
    print "***"
    ofile.close()
    
    
    return num_list
    



if __name__ == '__main__':
    file_name = 'ecoli_operon_walk.out'
    num_list = []
    cut_off = 0.9
    above_cut_off_group = []
    
    
    parsed_list_organism = parse_shuffled_trp(file_name,'organism')
    sorted_list_organism = sorted(parsed_list_organism, key = itemgetter(1), reverse=True)
    
    #parsed_list_between_operon = parse_shuffled_trp('ecoli_between_operon.out','organism_between_operon')
    #sorted_list_between_operon = sorted(parsed_list_between_operon, key = itemgetter(1), reverse=True)
    
    #parsed_list_operon = parse_shuffled_trp('possible_trp.out','operon')
    #sorted_list_operon = sorted(parsed_list_operon, key = itemgetter(1), reverse=True)
    
    #num_list_operon = write_sorted_list('trp_frac_all_possibilities.txt',sorted_list_operon,-1,'operon')
    
    #num_list_within_operon = write_sorted_list('ecoli_walk_frac.txt',sorted_list_organism,0)
    #num_list_between_operon = write_sorted_list('ecoli_between_walk_frac.txt',sorted_list_between_operon,-1)
    #num_list_between_operon_nozero = write_sorted_list('ecoli_between_walk_frac_nozero.txt',sorted_list_between_operon,0.1)
    
    #bins = 500
    #fig = p.figure(figsize=(15,15))
    #ax1 = fig.add_subplot(1,1,1)
    #n,bins,patches = p.hist(num_list_operon,bins,color='g')
    #p.xlabel('fraction of conserved gene_pair to the total number of combined gene occurrences')
    #p.ylabel('frequency')
    ##p.title('Histogram of conservation of gene pairs of Tryptophan operon (includes only gene pairs that were observed in atleast one organism')
    ##fig_name = 'tryptophan_operon_shuffled_zerocutoff.png'
    #p.title('Histogram of fraction of conserved gene pairs of Tryptophan operon with all possible shuffles')
    
    
    #bins = 500
    #fig = p.figure(figsize=(15,15))
    #ax1 = fig.add_subplot(2,1,1)
    #n,bins,patches = p.hist(num_list_between_operon,bins,color='g')
    #p.xlabel('fraction of conserved gene_pair to the total number of combined gene occurrences')
    #p.ylabel('frequency')
    #p.title('Histogram of conservation of gene pairs between Ecoli operons')
    
    ## This is the inset on between operon gene pair graph
    #ax3 = p.axes([.65, .6, .2, .2])
    #n,bins,patches = p.hist(num_list_between_operon_nozero,bins,color='g')
    #p.title('Histogram with fraction score cutoff > 0')
    #p.setp(ax3, xticks=[0.2,0.5,1.0])
    
        
    ##ax2 = p.twinx()
    #ax2 = fig.add_subplot(2,1,2)

    #n,bins,patches = p.hist(num_list_within_operon,bins,color='r')
    #p.xlabel('fraction of conserved gene_pair to the total number of combined gene occurrences')
    #p.ylabel('frequency')
    #p.title('Histogram of conservation of gene pairs within Ecoli operons')
    
    
    
    
    fig_dir = os.getcwd() + '/figures/'
    fig_name = fig_dir + 'trp_all_shuffles.png'
    
    #p.show()
    p.savefig(fig_name,format = "png", orientation='landscape')

