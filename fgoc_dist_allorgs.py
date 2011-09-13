#! usr/bin/python

# This script creates a graph of the dissociation of fgoc score (gen-fgoc scores) with distance. (First - within operon pairs)

import os
import sys
import numpy as np
import pylab as p
import matplotlib


def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened : %s "%(name_file,err.strerror)
	sys.exit(0)

    return file

def create_dic_acc_gen_fgoc(file_name):
    """
    @param file_name: The file = orgs_ref_ecoli.genome_fgoc or other file name
    @function: Creates a key : value pair of accession number : genome_fgoc score that was calculated previously
    """
    
    # (1) Open file - orgs_ref_ecoli.genome_fgoc; This file is for genome_fgoc within operon
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    
    # (2) Initializing
    acc_genome_fgoc_dictionary = {}
    
    for line in lines:
	# (3) Makes a dictionary of acc_no:genome_fgoc
	acc_nbr, genome_fgoc = line.split('\t')[0], line.split('\t')[4]
	acc_genome_fgoc_dictionary[acc_nbr] = genome_fgoc

    return acc_genome_fgoc_dictionary

def make_distfile_ecoli(file_name):
    """
    @param file_name: The distance file. This is acc_nbr1 acc_nbr2 distance; This file can be enormous
    @function: Parse through the file and write lines containing ecoli acc_nbr in either the 1st or 2nd column. 
    """
    # (1) Initializing
    out_lines = []
    
    # (2) Iterate through lines
    for line in open(file_name,'r'):	# This is an iterator. It keeps only one line in memory at a time
	# (3) Obtain the three columns. 
	try:
	    acc_nbr1, acc_nbr2, distance = line.split(' ')[0], line.split(' ')[1], line.split(' ')[2]
	except IndexError:
	    print line
	    raise

	# (4) Check if either of the two accession numbers is Ecoli = NC_000913
	if (acc_nbr1 == 'NC_000913') or (acc_nbr2 == 'NC_000913'):
	    out_lines.append(line)
	
    # (5) Write the lines into output_file
    output_file = 'ecoli_orgs.dist'
    ofile = open_file(output_file,'w')
    for line in out_lines:
	ofile.write(line)
    ofile.close()
    
    return True

def create_dict_acc_dist(distance_file):
    """
    @param distance_file: This is the distance file containing only E.coli  someother orgs distance
    @function: Makes a dictionary other orgs: distance
    """
    # Initializing
    org_dist_dictionary = {}
    
    # (1) Open the distance file
    for line in open(distance_file,'r'):
	org1,org2,distance = line.split(' ')
	distance = distance[:-1]
	# (2) Either of the organism could be NC_000913; So an appropriate check and dictionary update
	if org1 == 'NC_000913':
	    org_dist_dictionary[org2] = distance
	if org2 == 'NC_000913':
	    org_dist_dictionary[org1] = distance

    return org_dist_dictionary

def create_dict_org_normgenfgoc(file_name):
    """
    @param file_name: This is the file *.genome_fgoc where the first column is the accession number and the 4th column is the normalised gen_fgoc_score
    @function: Parses the file and makes a dictionary acc_nbr: norm_genome_fgoc
    """
    # Initializing
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    org_normgenfgoc_dictionary = {}
    
    for line in lines:
	acc_nbr, norm_genome_fgoc = line.split('\t')[0], (line.split('\t')[4])
	org_normgenfgoc_dictionary[acc_nbr] = norm_genome_fgoc
    
    return org_normgenfgoc_dictionary

def get_xy_list(org_dist_dict, org_fgoc_dict):
    """
    @param org_dist_dict: This is the dictionary organism acc_nbr: distance
    @param org_fgoc_dict: This is the dictionary organism acc_nbr: normalised genome_fgoc_score
    @function: Iterate through the org_fgoc_dict and map the org acc_nbr to obtain distance. Append the x and the y list simultaneously
    # Note of caution: Many of the dictionary may not be complete. There may be no mapping at all. Key error will be introduced which will be passed
    """
    # Initializing
    distance_array = []
    gen_fgoc_array = []
    
    # (1) Iterate through the org_fgoc_dict
    for org, fgoc in org_fgoc_dict.items():
	# (2) Obtain distance from the org_dist_dictionary
	try:
	    distance = org_dist_dict[org]
	    # (3) Append the x and the y array; Convert to float
	    distance_array.append(float(distance))
	    gen_fgoc_array.append(float(fgoc))
	except KeyError:
	    pass
	

    return distance_array, gen_fgoc_array
    


def main(argument):
    """
    """
    [argv] = argument
    if argv == '1':
	
	# (1) Create a dictionary of acc_no: genome_fgoc; File used : orgs_ref_ecoli.genome_fgoc. This is within operon calculation
	file_name = 'orgs_ref_ecoli_.genome_fgoc'	# This file is within operon
	acc_genome_fgoc_dictionary = create_dic_acc_gen_fgoc(file_name)
	print "Dictionary of accession_nbr: genome_fgoc created"
	
	# (2) Make a distance file containing only the E.coli organism distance with others. 
	distance_file = 'prok_reduced_16SrRNA.filter.dist'
	print "Making distance file"
	make_distfile_ecoli(distance_file)

    if argv == '2':
	
	distance_file = 'ecoli_orgs.dist'
	
	# (3) Create a dictionary of acc_no: distance (only between ecoli)
	org_dist_dictionary = create_dict_acc_dist(distance_file)
	
	# (4) Create a dictionary of acc_no: norm_genome_fgoc for within operon
	file_name = 'within_orgs_ref_ecoli_ver6.genome_fgoc'
	within_org_norm_genfgoc_dictionary = create_dict_org_normgenfgoc(file_name)
	print within_org_norm_genfgoc_dictionary['NC_000913']
	dist_list1, gen_fgoc_list1 = get_xy_list(org_dist_dictionary, within_org_norm_genfgoc_dictionary)
	
	# (5) Create a dictionary of acc_no: norm_genome_fgoc for within operon
	file_name = 'bw_operons_orgs_ref_ecoli_ver6.genome_fgoc'
	bw_org_norm_genfgoc_dictionary = create_dict_org_normgenfgoc(file_name)
	dist_list2, gen_fgoc_list2 = get_xy_list(org_dist_dictionary, bw_org_norm_genfgoc_dictionary)
	print len(dist_list1)

	# Plot graph
	fig = p.figure(figsize=(15,15))
	ax1 = fig.add_subplot(1,1,1)
	p.xlabel('# Distance (from 16SrRNA) ')
	p.ylabel('# Normalised Genomic FGOC score')
	p.title('Decay of Gene order with genetic distance')
	scatter1 = p.scatter(dist_list1,gen_fgoc_list1,color='g',marker='^')
	scatter2 = p.scatter(dist_list2,gen_fgoc_list2,color='r',marker='o')
	#p.legend((line1,line2,line3, line4),('E cut; K,R labeled', 'E cut; K, R, C labeled', 'P cut; K, R labeled','P cut; K, R, C labeled'), loc=4)
	p.legend((scatter1,scatter2),('pairs within operons', 'pairs between operons'), loc=4)
	fig_dir = os.getcwd() + '/figures/'
	fig_name = fig_dir + 'gen_fgoc_with_distance(homology).png'
	p.savefig(fig_name,format = "png", orientation='landscape')
	#p.show()
	
    

if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - fgoc_dist_allorgs.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    