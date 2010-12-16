#!/usr/bin/python
# The script parses the file - "PRK_Clusters.bcp" and creates a table in sqlite3. It also creates a text file with the table (just in case).The table contains the following columns - 
# cluster accession - gene_id - locus_tag - protein_gi - protein_name - organism_name
# The file is in similar format. Each cluster description is separated by a ////. Conveniently there is good delimitation and the cluster names starts off with ENTRY
# Also the clusters have all the protein information that starts off with Proteins. 

import os
import sys
import sqlite3
from warnings import filterwarnings
import timeit

def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """
    
    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	sys.exit(0)
    
    return file

def write_table(protein_info):
    """
    @param input:	protein_info - This contains the information of one single protein information with the cluster accession
    @function:		This just inserts the protein_information into the table
    @param output:	None
    """
    table_name = 'protein_clusters'
    [cluster_accession, gene_id, locus_tag, protein_gi, protein_name, organism_name] = protein_info
    
    insert_table ='''
    INSERT INTO protein_clusters
    (cluster_accession, gene_id, locus_tag, protein_gi, protein_name, organism_name)
    VALUES
    ('%s','%s','%s','%s','%s','%s')
    '''%(cluster_accession, gene_id, locus_tag, protein_gi, protein_name, organism_name)
    cursor.execute(insert_table)
    #try:
	#cursor.execute(insert_table)
    #except:
	#print insert_table
    return

def get_cluster_details(cluster):
    """
    @param input: 	cluster: this is a large string with \n and contains all the information of one cluster
    @function:		The function parses this cluster to retrieve the cluster number and all the proteins with their details in that particular cluster
    @param output:	A list of all the protein details with the extra column of the cluster_number
    """
    #Initializing
    flag_proteins = False
    all_clustered_proteins = []
    gene_id = locus_tag = protein_gi = protein_orgs_info = protein_name = organism_name = ''
    split_cluster = cluster.split('\n')[:-1]	#This last list is an empty string and is removed
    
    for line in split_cluster:
	if line.startswith('ENTRY'):
	    cluster_accession = line.split('\t')[1]
	if line.startswith('PROTEINS'):		#Since this is the last identifier, only now will the flag_proteins be flagged True.
	    flag_proteins = True
	if (flag_proteins):
	    if not line.startswith('PROTEINS'):
		split_line = line.split('\t')
		gene_id, locus_tag, protein_gi, protein_orgs_info = split_line[1],split_line[2],split_line[3],split_line[4]
		locus_tag = locus_tag.lstrip()	#Somehow there are spaces in the locus_tag and this is to remove them. 
		protein_name, organism_name = protein_orgs_info.split('[')[0], protein_orgs_info.split('[')[1][:-1]
		organism_name = organism_name.replace('\'','_')	#Remember these changes. Will need to do it in other places when the organism name needs to be matched up. 
		protein_name = protein_name.replace('\'','_')
		one_protein_info = [cluster_accession, gene_id, locus_tag, protein_gi, protein_name, organism_name]
		write_table(one_protein_info)
		
		all_clustered_proteins.append(one_protein_info)
		
    return all_clustered_proteins
    

def parse_clusters(file_name):
    """
    @param input:	 file_name: The file downloaded from NCBI - 'PRK_Clusters.bcp' It contains each cluster that ends by '///' ; ENTRY begins the line that tells the cluster number. A line that begins with PROTEINS starts the list of all proteins and its details that belong to the cluster. 
    @function: 		Breaks up the file into Clusters and then makes a list containing all the protein details along with a column for a cluster. 
    @param output:	A gigantic list containing all the details of each cluster
    """
    
    #Initializing
    ifile = open_file(file_name)
    all_info = ifile.read()
    ifile.close()
    flag_proteins = False
    clusters = ''
    line = ''
    
    all_clusters = all_info.split('////')
    for cluster in all_clusters:
	clustered_proteins = get_cluster_details(cluster)

    return
def operate_table(table_name,contents,operation='insert'):
    """Creating Table - The details of the table columns are passed on as the string with contents. The contents are passed as within the triple quote format
    """
	
    table_contents = contents
    if operation == 'create': 
	cursor.execute(""" DROP TABLE IF EXISTS %s """%table_name)	#This works well only when the warning sign is turned off
    else: 
	pass 
    
    print table_name
    cursor.execute(table_contents)
    
    return 
    

def create_table(table_name):
    """
    @param input:	table_name - Inputs the name of the table to be created
    @function:		This function specifically creates the protein_clusters table with the following columns - 
    cluster_accession VARCHAR(10), gene_id INT (10), locus_tag VARCHAR (25), protein_gi INT (10), protein_name MEDIUMINT, organism_name MEDIUM INT
    @param output:	None
    """
    
    contents =''' CREATE TABLE protein_clusters
    (cluster_accession VARCHAR(10),
    gene_id INT (10),
    locus_tag VARCHAR (25), 
    protein_gi INT (10),
    protein_name MEDIUMINT,
    organism_name MEDIUMINT
    )
    '''
    
    operate_table(table_name, contents,'create')
    
    return


if __name__ == '__main__':
    
    # Calling database using sqlite3 module
    conn = sqlite3.connect('bsub.db')
    cursor = conn.cursor()
    
    table_name = 'protein_clusters'
    file_name = '/project/marcotte/jagannath/projectfiles/Downloads/PRK/PRK_Clusters.bcp'
    
    create_table(table_name)
    print "Processing PRK_Clusters.bcp and Creating Table .."
    parse_clusters(file_name)
    
    cursor.close()
    conn.commit()
    conn.close()
    
  