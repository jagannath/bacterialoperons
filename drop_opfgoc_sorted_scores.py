#! usr/bin/python

from __future__ import division
import os
import sys
import pickle
import numpy as np
import pylab as p
import matplotlib
from operator import itemgetter

# This programme calculates the percentage drop in the operonic fgoc score between the first ranked and the next ranked operon order. A histogram of the drop may indicate which operons have the greatest drop % and those operon orders that have lowest. 


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
    
def calculate_score_drop(file):
    """
    @file: The operon-name_shuffled.sort file. Contains operon_order op-fgoc score
    @function: Calculate the percentage score drop in the operon order. 
    """
    file_name = os.getcwd() + '/ecoli_homolog_shuffled/' + file
    ifile = open(file_name)
    lines = ifile.readlines()
    ifile.close()
    flag = True
    order_list = []
    count = 0

    first_order, first_opfgoc = lines[0].split('\t')
    for line in lines:
	if flag:
	    second_order, second_opfgoc = line.split('\t')
	    if first_opfgoc > second_opfgoc:
		flag = False
    
    first_opfgoc = float(first_opfgoc)
    second_opfgoc = float(second_opfgoc)
    operon_name = file.split('_')[0]
    nbr_genes = len(eval(first_order))
    try:
	drop_perc = ((first_opfgoc - second_opfgoc)/first_opfgoc)*100
    except ZeroDivisionError:
	print operon_name
	count+=1
	drop_perc = 0
    drop_wrt_genenbr = drop_perc * nbr_genes
    order_list.append([first_order, second_order])
    
    return [operon_name, drop_perc, drop_wrt_genenbr, order_list, nbr_genes]

def main(argument):
    
    [argv] = argument
    
    if argv == '1':
	# This operates on the ecoli_shuffled.sort files and obtains the difference in rank. The operon name, and the order is remembered. 
	# (1) Obtain the list of sort files
	shuffled_dir = os.getcwd() + '/ecoli_homolog_shuffled'
	sort_files = [file for file in os.listdir(shuffled_dir) if file.endswith('sort')]

	# (2) Iterate through the list and obtain the operon_name, score drop (in %), [order(rank1),order(rank2)], number of genes in the operon
	# Initializing
	all_operon_sort_info = []
	operon_drop_details_dict = {}
	for file in sort_files:
	    [operon_name, drop_perc, drop_wrt_genenbr, order_list, nbr_genes] = calculate_score_drop(file)
	    all_operon_sort_info.append([operon_name, drop_perc, drop_wrt_genenbr, order_list, nbr_genes])
	    # Make a dictionary of operon_name : all drop details
	    operon_drop_details_dict[operon_name] = [operon_name, drop_perc, drop_wrt_genenbr, order_list, nbr_genes]
	    
	# (3) Pickle the dictionary
	file_name = 'operon_drop_details_dict.pkl'
	pickle_file(file_name, operon_drop_details_dict)
	
    if argv == '2':
	# This part analyses the drop perc and makes histogram and associates operons that behave differently
	# (1) Obtain list of only drop_perc
	file_name = 'operon_drop_details_dict.pkl'
	operon_drop_details_dict = unpickle_file(file_name)
	
	count = 0
	all_operon_sort_info = operon_drop_details_dict.values()
	all_drop_perc = [(item[0], item[1]) for item in all_operon_sort_info if (item[1] < 100 and item[1] > 0 and item[4] > 2)]
	
	# Ranking them 
	sorted_drop_perc = sorted(all_drop_perc, key = itemgetter(1), reverse=True)
	print sorted_drop_perc[0:10]
	
	sys.exit(1)
	
	print all_drop_perc
	sys.exit(1)
	#all_drop_perc = [item[1] for item in all_operon_sort_info if item[4] == 3]
    
	for operon_name, drop_details in operon_drop_details_dict.items():
	    if drop_details[1] < 100 and drop_details[1] > 0 and drop_details[4] > 3:
		print operon_name, drop_details
		count+=1 
	
	print count
	# (2) Draw Histogram of drop_perc
	bins = 500
	fig = p.figure(figsize=(15,15))
	ax1 = fig.add_subplot(1,1,1)
	n,bins,patches = p.hist(all_drop_perc,bins,color='g')
	p.xlabel('Percentage drop of Operonic FGOC score (from highest to the next highest ranked operon order)')
	p.ylabel('frequency')
	p.title('Histogram of Percentage drop of Operonic FGOC score across Ecoli operons (excluding 100 and 0)')
	#ax2 = p.axes([.5, .5, .3, .3])
	#p.title('Histogram excluding zero')
	##ax2 = fig.add_subplot(0.65,0.6,0.2)
	#n,bins,patches = p.hist(cutofflist,bins,color='b')
	#p.setp(ax2)
	fig_dir = os.getcwd() + '/figures/'
	fig_name = fig_dir + 'ecoli_operon_drop_perc_ver2.png'
	#p.savefig(fig_name,format = "png", orientation='landscape')
	#p.show()
    




if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - drop_opfgoc_sorted_scores.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    