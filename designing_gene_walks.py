#! usr/bin/python

# Designing the walk from a given gene. For a given gene, it will obtain all the possible connections it has, what kind and the weight for each. It will then calculate the maximum of all the connections and use the next bnumber as the next gene

from __future__ import division
import os
import sys
try:
   import cPickle as pickle
except:
   import pickle
import shelve
import re
import itertools
import math
import numpy as np
import pylab as p
import matplotlib
from operator import itemgetter

def get_maxfgoc_details(line):
    """
    @param line: Line from the .fgoc file. 
    @function: Parses the line to obtain the source gene, target gene, maximum of the fgoc scores and the type of connection corresponding to that - i.e. directed, convergent or divergent
    """
   
    source, target, dir_fgoc, conv_fgoc, div_fgoc = line.split('\t')[0], line.split('\t')[1], line.split('\t')[3], line.split('\t')[4], line.split('\t')[5]
    max_fgoc = max([dir_fgoc,conv_fgoc,div_fgoc])
    max_index = [dir_fgoc, conv_fgoc, div_fgoc].index(max_fgoc)
    status = (['directed','convergent','divergent'])[max_index]
    
    return source, target, max_fgoc, status

def get_connection_dict(file_name):
    """
    @param file_name: File containing the various fgoc scores and the connection status
    @function: Makes a dictionary of b_nbr: [b_nbr_connected, fgoc, connection_status]; The connection status and fgoc is determined by the maximum of the fgoc values 
    """
    
    ifile = open(file_name)
    lines = ifile.readlines()
    ifile.close()
    
    dict_connects = {}
    for line in lines:
	status = False
	source, target, max_fgoc, connection_status = get_maxfgoc_details(line)
	#print source, target, max_fgoc, connection_status
	try:
	    prev_list = dict_connects[source]
	    status = True
	except KeyError:
	    dict_connects[source] = [[target, max_fgoc, connection_status]]

	if status:
	    prev_list = dict_connects[source]
	    prev_list.append([target, max_fgoc, connection_status])
	    dict_connects[source] = prev_list
    
    return dict_connects

def main(argument):
    [argv] = argument
    
    if argv == '1':
	# (1) Make a dictionary of bnumber : [adj gene1, adj gene2]. Each adj gene is defined as [bnumber, kind of connection (directed, convergent or divergent) and weight of it]
	file_name = os.getcwd() + '/directed_fgoc_ecoli_ko/' + 'ecoli_kegg_dir_convdiv.fgoc'
	dict_connections = get_connection_dict(file_name)
	
	# (2) For a given bnumber, which is the most prefered connection
	bnbr = 'b4237'
	all_connections = dict_connections[bnbr]
	all_fgoc_scores = [float(item[1]) for item in all_connections]
	print all_fgoc_scores
	print all_connections
	sort_scores = sorted(all_fgoc_scores, key = itemgetter, reverse=True)
	print sort_scores
	first_index = all_fgoc_scores.index(sort_scores[0])
	second_index = all_fgoc_scores.index(sort_scores[3])
	third_index = all_fgoc_scores.index(sort_scores[4])
	print all_connections[first_index]
	print all_connections[second_index]
	print all_connections[third_index]
	
	    


if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - designing_gene_walks.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))