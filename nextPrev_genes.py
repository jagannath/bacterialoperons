#! usr/bin/python

# Uses the pickle file for next genes and prev genes in all organisms and remakes all them into 2 shelve files

from __future__ import division
import os
import sys
try:
    import cPickle as pickle
except:
    import pickle
import shelve
import re
import random
import pdb
from operator import itemgetter
import fnmatch







def pickle_file(file_name, contents, path_dir=os.getcwd()):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    #pkl_dir = os.getcwd() + '/pkl_files/'
    pkl_dir = path_dir + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ofile = open(pkl_file,'wb')
    pickle.dump(contents, ofile)
    ofile.close()
    
    return True

def unpickle_file(file_name, path_dir=os.getcwd()):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    pkl_dir = path_dir + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ifile = open(pkl_file)
    contents = pickle.load(ifile)
    ifile.close()
    
    return contents

def update_shelve_next(nextGene_dict):
    """ Updates the nextShelve file to contain all the details of the nextGene_dict """
    lastGeneRank = 0
    lastGene = []
    firstGene = []
    for key, cdsDetail in nextGene_dict.items():
	if int(cdsDetail[1]) > lastGeneRank:
	    lastGeneRank = int(cdsDetail[1])
	    lastGene = cdsDetail
	if int(cdsDetail[1]) == 1:firstGene = cdsDetail
	nextGeneShelve[key] = cdsDetail
    nextGeneShelve[lastGene[0]] = firstGene
    return 

def update_shelve_prev(prevGene_dict):
    """ Updates the nextShelve file to contain all the details of the nextGene_dict """
    lastGeneRank = 0
    lastGene = []
    firstGene = []
    for key, cdsDetail in prevGene_dict.items():
	if int(cdsDetail[1]) > lastGeneRank:
	    lastGeneRank = int(cdsDetail[1])
	    lastGene = cdsDetail
	if int(cdsDetail[1]) == 1:firstGene = cdsDetail
	prevGeneShelve[key] = cdsDetail
    prevGeneShelve[firstGene[0]] = lastGene
    return 


def main(argument):
    [argv] = argument
    
    if argv == '1':
	# Using os.walk to find the files *.fasta
	rootPath = '/project/marcotte/jagannath/projectfiles/bacterialoperons/parsed_gbk_files_genomes/'
	pattern = '*_next_gene_dict.pkl' # Can include any UNIX shell-style wildcards
	count = 0
	for root, dirs, files in os.walk(rootPath):
	    for filename in fnmatch.filter(files, pattern):
		count+=1
		nextGeneFile = os.path.join(root, filename)
		print "Processing %d : %s ..."%(count,filename)
		nextGene_dict = pickle.load(open(nextGeneFile))
		update_shelve_next(nextGene_dict)
		print len(nextGeneShelve)
    
    if argv == '2':
	rootPath = '/project/marcotte/jagannath/projectfiles/bacterialoperons/parsed_gbk_files_genomes/'
	pattern = '*_prev_gene_dict.pkl' # Can include any UNIX shell-style wildcards
	count = 0
	for root, dirs, files in os.walk(rootPath):
	    for filename in fnmatch.filter(files, pattern):
		count+=1
		prevGeneFile = os.path.join(root, filename)
		print "Processing %d : %s ..."%(count, filename)
		prevGene_dict = pickle.load(open(prevGeneFile))
		update_shelve_prev(prevGene_dict)
		
    if argv == '3':
	# Redraws the tree by mapping the taxonID with three Letter code
	path_dir = '/project/marcotte/jagannath/projectfiles/Downloads/Strings_database'
	file_name = 'species.tree.v9.0.txt'
	ifile = open(path_dir + '/' + file_name)
	treelines = ifile.readlines()
	ifile.close()
	
	path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/graph_genomes/all_orgs_dictionaries'
	file_name = 'accNbr_threeCode_dict.pkl'
	accNbr_threeCode_dict = unpickle_file(file_name, path_dir)
	file_name = 'TaxonID_AccNbr_dict.pkl'
	taxonID_AccNbr_dict = unpickle_file(file_name, path_dir)
	
	allLines = []
	file_name = 'template.tree'
	ofile = open(file_name,'w')
	pat = re.compile('\s+\d+')
	for line in treelines:
	    status = True
	    if pat.match(line):
		nbr = pat.match(line).group()
		nbr = nbr.lstrip()
		try:
		    accNbr = taxonID_AccNbr_dict[nbr]
		    threeCode = accNbr_threeCode_dict[accNbr]
		    aline = line.replace(nbr,threeCode)
		    line = aline
		except KeyError: line = line
	    allLines.append(line)
	    ofile.write(line)
	ofile.close()
	file_name = 'template.tree'
	pickle_file(file_name,allLines, path_dir)





if __name__ == '__main__':
    
    argument = sys.argv[1:]
    #pdb.set_trace()
    print "Processing %s ..."%(argument)
    # Setting some global dictionaries
    # Making a nextShelve file
    filename = 'nextGene.shelve'
    nextGeneShelve = shelve.open(filename)
    filename = 'prevGene.shelve'
    prevGeneShelve = shelve.open(filename)
    
    main(argument)
    nextGeneShelve.close()
    prevGeneShelve.close()
    import time
    print "Script - orthogrouping_newgenes.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))