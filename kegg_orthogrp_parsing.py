#!/usr/bin/python

# This programme parses the Kegg Orthogroup database (downloaded Apr 21st 2011). Makes a dictionary pair of group_name : list of all locus_tags. Then another part will search ones with ECO id or b# and link b# to list of all locus_tags. This will filter out the ECO gene from the orthogroup

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


def get_gene_detail(line):
    """
    @param line: This is the line containing the org_id: gene (gene_name sometimes) or number of genes separated by space
    @function: Obtain the list of genes in this line.
    @returns : A list
    """
    line = line.lstrip()
    org_id,genelist = line.split(':')
    # It still leaves a trailing space genelist
    genelist = genelist.lstrip()
    genes = genelist.split(' ')	#Separates if there are a number of genes
    all_genes = []
    for gene in genes:
	# Check if there is a ( in the gene. If yes, take the part that is before (; (gene name) is provided
	if re.search('\(',gene):
	    gene = gene.split('(')[0]
	all_genes.append(gene)

    return org_id, all_genes

def get_genes(group):
    """
    @param group: This is the set of lines belonging to a group
    @function: Parses the group to return back the group_name, list_genes and ecoli b_number
    """
    #Initializing
    lines = group.split('\n')
    all_locus_tags = []
    b_number = ['Null']	#Default if no Ecoli gene
    flag = False
    group_name = 'Null'
    include_orgs=['ECJ','ECD','EBW','ECE','ECS','ECF','ETW','EOJ','EOI','EOH','ECG','EOK','ECC','ECP','ECI','ECV','ECX','ECW','ECM','ECY','ECR','ECQ','ECK','ECT','EUM','ECZ','ECL','EBR','EBD','EFE']
    exclude_orgs = []
    
    for line in lines:
	# ENTRY is the key word for Kegg orthogroup number
		if line.startswith('ENTRY'):
			group_name = line.split('       ')[1]
			print "Processing %s ..."%(group_name)
		# GENES is the key word for the list of locus_tags. Unfortunately the first one is in line with the GENES word. Replacing the word GENES with space. Then removing the space to obtain the relevant details
		if line.startswith('GENES'):
			flag = True	#Flag to tell it encountered the beginning of gene list
			line = line.replace('GENES','     ')
		if flag and line.startswith('            ') and not line.count(':')>1:
			# Obtain details of the gene
			org_id, genelist = get_gene_detail(line)
			# If the aim is to make just ecoli orthogroups then run the b_number etc 
			if (make_ecoli_orthogrps):
			    if org_id == 'ECO':
				    b_number = genelist
			    elif org_id in include_orgs and not org_id in exclude_orgs:
				    print org_id, genelist
				    for gene in genelist:
					    all_locus_tags.append(gene)
			    # Obtain the Ecoli b_number if exist // Must also include (maybe later) to neglect all particular org_id
			    # There are two lists - include orgs list and exclude orgs list. 
  
    return group_name, all_locus_tags, b_number

def get_adjacent_genepairs(ecoli_gene_list):
    """
    """
    # Initializing
    adjacent_genepairs = []
    
    for i in xrange(1,len(ecoli_gene_list)):
		gene_1 = ecoli_gene_list[(i-1)]
		gene_2 = ecoli_gene_list[i]
		adjacent_genepairs.append([gene_1,gene_2])

    # Since it is circular, the gene pair of last gene and 1st gene must be considered. 
    gene_1 = ecoli_gene_list[0]
    gene_last = ecoli_gene_list[len(ecoli_gene_list)-1]
    adjacent_genepairs.append([gene_last,gene_1])

    return adjacent_genepairs

def main(argument):
    [argv] = argument
    
    if argv == '1':
	# This part will parse the file - kegg_orthogroups
	kegg_ecoli_cluster = {}
	# (1) Read the kegg_orthogroups file
	#path_dir = '/project/marcotte/jagannath/projectfiles/Downloads/'
	path_dir = os.getcwd()
	file_name = '/kegg_orthogroups'
	ifile = open(path_dir + file_name)
	
	# (2) Separator for each group in the file is ///. Using this to separate the group
	groups = (ifile.read()).split('///')
	ifile.close()
	
	# (3) Iterate through each group and get a list of all genes belonging to group, the ecoli b_number if exist and group_name; Note that Ecoli gene is not included in the orthogroup
	red_bnbr = []
	all_bnbr_groups = []
	group_count = 0
	# I will randomise the group number so that adjacent_genepairs arent next to one another in group numbers
	random_range = [i for i in xrange(1,3000)]
	random.shuffle(random_range)
	file_name = 'kegg_ecoli_bnbr.clust'
	ofile = open(file_name,'w')
	for group in groups:
		group_name, list_genes, b_number = get_genes(group)
		all_bnbr_groups.append(b_number)
		if len(b_number) > 1:
			red_bnbr.append(b_number)
		if not b_number[0] == 'Null':
			group_count += 1
			print group_count

			for ecoligene in b_number:
				try:
					cluster_name = ecoligene + '_' + group_name
				except:
					print b_number, group_name, ecoligene
					sys.exit(1)
				
				#info_towrite = ecoligene + '\t' + str(random_range[group_count]) + '\n'
				#ofile.write(info_towrite)
				kegg_ecoli_cluster[cluster_name] = list_genes
				path_dir = os.getcwd() + '/ecoli_kegg'
				file_name = ecoligene + '_onlyecoli_ko.pkl'
				pickle_file(file_name, list_genes, path_dir)
				ofile.close()
				
				print red_bnbr
				
				# (4) Pickle the kegg_ecoli_cluster
				file_name = 'kegg_ecoli_cluster_list_genes_dict.pkl'
				path_dir = os.getcwd() + '/ecoli_kegg/'
				pickle_file(file_name, kegg_ecoli_cluster, path_dir)
				
				file_name = 'kegg_ecoli_within_orthologs.pkl'
				path_dir = os.getcwd() + '/ecoli_kegg/'
				pickle_file(file_name, red_bnbr, path_dir)
				
		else: # This part is because I opted to get all the locus ids for the kegg orthogroups
		    kegg_clusters[group_name] = list_genes 
		
		    # (4) Pickle the kegg_clusters dictionary
		    file_name = 'kegg_clusters_dict.pkl'
		    path_dir = os.getcwd() + ''

	

    if argv == '2':
	# This makes list of all 2993 vs 2993 gene pairs. But two types of genes (of E.coli) are flagged. Ones where a genes of Ecoli are missing and another in the 2993 gene list in which the gene upstream or downstream is missing. 
	# (1) Getting Ecoli gene list
	file_name = 'ecoli_homolog_data.txt'
	ifile = open(file_name)
	lines = ifile.readlines()
	ifile.close()
	ecoli_gene_list = [line.split('\t')[1] for line in lines]
	
	# (2) Getting the list of Kegg ortholog ecoli genes
	path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/ecoli_kegg/pkl_files'
	file_list = os.listdir(path_dir)
	ko_gene_list = [item[:-4] for item in file_list]
	
	# (3) Get the list of missing genes i.e those that are present in the ecoli_gene_list but not in ko_gene_list
	missing_genes = [gene for gene in ecoli_gene_list if not gene in ko_gene_list]
	check_missing = [gene for gene in missing_genes if gene in ko_gene_list]
	
	# (4) Flag those genes in the KO_genelist if they are missing in known upstream genes or downstream genes
	upstream_flag = []
	downstream_flag = []
	adjacent_genepairs = get_adjacent_genepairs(ecoli_gene_list)

	count = 0
	missing_genepairs = []
	for genepair in adjacent_genepairs:
	    [gene_1, gene_2] = genepair
	    if gene_1 in missing_genes and gene_2 in missing_genes:
			missing_genepairs.append([gene_1,gene_2])
	    else:
		if gene_1 in ko_gene_list and gene_2 in missing_genes:
		    downstream_flag.append(gene_1)
		if gene_1 in missing_genes and gene_2 in ko_gene_list:
		    upstream_flag.append(gene_2)
		if gene_1 in ko_gene_list and gene_2 in ko_gene_list:
		    count+=1
		    print genepair
	
	
	print len(upstream_flag)
	print len(downstream_flag)
	print count
	print len(ko_gene_list)
	print len(missing_genepairs)
	print len(ecoli_gene_list)

    if argv == '3':
	# Make b_nbr.pkl files now in kegg_ecoli_details directory; It contains list of CDS details of all genes belonging to the homologue
	# (1) Obtain list of 2991 ecoli genes
	path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/ecoli_kegg/onlyEcoli_details/bnbr_onlyecoli_list/pkl_files'
	file_list = os.listdir(path_dir)
	ko_gene_list = [item[:-17] for item in file_list]
	# (2) For each gene in list, use the CDS shelve to get details of the gene
	shelve_dir = os.getcwd() + '/shelve_files/'
	locus_cds_dictionary = shelve.open(shelve_dir + 'locus:cds_information_dictionary')
	path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/ecoli_kegg/onlyEcoli_details/bnbr_onlyecoli_list'
	write_path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/ecoli_kegg/onlyEcoli_details'
	for gene in ko_gene_list:
	    print "Processing %s ..."%(gene)
	    # (3) Unpickle the file
	    file_name = gene + '_onlyecoli_ko.pkl'
	    ortholog_list = unpickle_file(file_name,path_dir)
	    all_locus_details = []
	    for ortholog in ortholog_list:
		try:
		    locus_details = locus_cds_dictionary[ortholog]
		    all_locus_details.append(locus_details)
		except KeyError:
		    pass
	    # (4) Pickle the file bnumber.pkl to contain all the locus_details
	    file_name = gene + '.pkl'
	    pickle_file(file_name,all_locus_details,write_path_dir)
	
    if argv == '4':
	# This will create 8 .sh files to write down all the pairwise combinations of these 2991 genes. 
	# (1) Obtain all the 2991 genes
	path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/ecoli_kegg/pkl_files'
	file_list = os.listdir(path_dir)
	ko_gene_list = [item[:-4] for item in file_list]
	
	# (2) Write into output files (4 of them)
	index = 0
	count = 0
	lst_dir = os.getcwd() + '/sh_files/'
	for gene_1 in ko_gene_list:
	    for gene_2 in ko_gene_list:
		count+=1
		print "Processing %d %s"%(count,(gene_1 + ' ' + gene_2))
		index += 1
		if index == 9: index = 1
		file_name = lst_dir + 'fgoc_all_pairs_pkl_ko_' + str(index) + '.sh'
		info_line = 'python get_fgoc_pair.py ' + gene_1 + ' ' + gene_2 + '\n'
		ofile = open(file_name,'a')
		ofile.write(info_line)
		ofile.close()


if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - kegg_orthogrp_parsing.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    
