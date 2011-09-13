#!usr/bin/python

# This script finds the gene neighborhood of the Streptomyces genes and tries to bridge a gap to the nearest and most abundant very homologous ecoli genes; 
# Idea is to blast this unusual gene with all the orgs and get a homolog group of that particular gene; Will keep cut off later. But will keep a track of identity score. Also will need to make the ecoli homolog groups again but in the list have %identity, length and evalue as well. 
import os
import sys
import pickle
import shelve
shelve_dir = os.getcwd() + '/shelve_files/'
locus_cds_dictionary = shelve.open(shelve_dir + 'locus:cds_information_dictionary')










def main(argument):
    
    [argv] = argument
    goi_locus_tag = 'SGR_5936' # From Kegg database look up the gene name for the unique enzyme
    goi_info = locus_cds_dictionary[goi_locus_tag]
    org_nbr = goi_info[5]
    rank_goi = goi_info[1]
    window = 5
    for key,value in locus_cds_dictionary.items():
	if value[5] == org_nbr and rank in xrange(rank_goi-window, rank_goi+window):
	    print key
    
    #if argv == '1':
	# (1) 








if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    main(argument)

    import time
    print "Script - finding_gene_neighborhood.py - %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    
