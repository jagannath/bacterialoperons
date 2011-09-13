#! usr/bin/python

from __future__ import division
import os
import sys
try:
   import cPickle as pickle
except:
   import pickle
import shelve

# This programme is to redo the ecoli homologues into 4000 odd pkl files. 

print "yes"

pkl_dir = os.getcwd() + '/shelve_files/'
#cluster_dictionary = shelve.open(pkl_dir + 'ecoli_homolog_locus_tag:group_dictionary')
group_all_locus_dictionary = shelve.open(pkl_dir + 'ecoli_homolog_group:all_locus_tag_dictionary')	# Needed
locus_cds_dictionary = shelve.open(pkl_dir + 'locus:cds_information_dictionary')

def pickle_file(file_name, contents,path_dir=os.getcwd()):
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
    
def unpickle_file(file_name,path_dir=os.getcwd()):
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

def filter_out_orgs(genus_name,file_name):
    """
    @param genus_name: This is the genus that must be filtered out from the list of genes
    @param file_name: File name of all_filenames.lst that contains the .gbk information
    @function: Obtains all the acc_nbr of organisms in genus (say Escherichia)
    """
    acc_nbr_list = []
    ifile = open(file_name)
    lines = ifile.readlines()
    ifile.close()

    acc_nbr_list = [(line.split('/')[6])[:-5] for line in lines if (line.split('/')[5]).startswith(genus_name)]

    return acc_nbr_list


def main(argument):
    
    [argv] = argument

    if argv == '1':
	#group = 'ECHOM_b0002'
	# Obtain details of all genes belonging to this homolog group
	
	# (1) Get a list of all groups in E.coli
	# Open the file - ecoli_homolog_data.txt and obtain ecoli_homolog pairs
	file_name = 'ecoli_homolog_data.txt'
	ifile = open(file_name)
	lines = ifile.readlines()
	ifile.close()
	group_list = [line.split('\t')[7][:-1] for line in lines]	#[:-1] is to delete the carriage return
	
	# (2) Get list of all acc_nbr of genus Escherichia
	genus_name = 'Escherichia'
	file_name = 'all_filenames.lst'
	filtered_out_list = filter_out_orgs(genus_name,file_name)
	
	count = 0
	for group in group_list:
	    count+=1
	    # (3) Obtain list of all genes in this homolog group
	    try:
		list_genes = group_all_locus_dictionary[group]
	    except KeyError:
		list_genes = []
	    print "Processing %d %s ..."%(count,group)
	    
	    # (4) For each gene in the list obtain a list containing the CDS information
	    index = 0
	    locus_info = []
	    locus_info_filtered = []
	    for gene in list_genes:
		cds_info = locus_cds_dictionary[gene]
		if cds_info[5] in filtered_out_list:
		    locus_info.append(cds_info)
		else:
		    locus_info.append(cds_info)
		    locus_info_filtered.append(cds_info)

	    # (4) Pickle the two files
	    homolog_dir = os.getcwd() + '/homologs_ecoli/' 
	    file_name = group + '.pkl'
	    pickle_file(file_name,locus_info,homolog_dir)
	    homolog_filtered_dir = os.getcwd() + '/homologs_ecoli_filtered/' 
	    file_name = group + '-filtered_ecoli.pkl'
	    pickle_file(file_name,locus_info_filtered,homolog_filtered_dir)




if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - make_genevsgene_mat.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    