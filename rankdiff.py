#! bin/usr/python

# This script file opens up the ecoli_data.txt file which is again ecoli gene data with all identifiers including operon name and group id. 
# It - (1) Takes the group_id from ecoli genes in an operon (2) Finds the set of all members in the group and fetches their details (3) Computes their rank difference and distance between them. If their rank difference is one, then it makes a count. It prints a count against a walk


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

def sql_execute_two_sets(group_name_a, group_name_b):
    """ This function queries the sqlite3 database 'orthomcl.db' for all the sequence_ids in the group_id. It then computes the difference in ranking between the two sequence_ids between one specie (unique accno). If it is 1 (absolute 1 for now) then it increments a count. It returns back the frequency of counts obtained
    """
    
    #Initializing
    count = 0
    set1 = ''
    set2 = ''
    output_contents = ''
    orthologous_gene_pairs = ''
    distance = 0
    number_tried = 0
    tried_gene_pairs = ''
    
    output_contents = """SELECT c.sequence_id, c.rank, c.orientation, c.left_end, c.right_end, c.sequence_length, c.taxon_id, g.group_id
    FROM
    (SELECT group_id, sequence_id
    FROM groups
    WHERE group_id = '%s') g
    INNER JOIN coding_sequences c
    ON c.sequence_id = g.sequence_id
    """%(group_name_a)
    output_contents = output_contents.replace('\n','')
    cursor.execute(output_contents)
    set1 = cursor.fetchall()
    
    
    output_contents = """SELECT c.sequence_id, c.rank, c.orientation, c.left_end, c.right_end, c.sequence_length, c.taxon_id, g.group_id
    FROM
    (SELECT group_id, sequence_id
    FROM groups
    WHERE group_id = '%s') g
    INNER JOIN coding_sequences c
    ON c.sequence_id = g.sequence_id
    """%group_name_b
    output_contents = output_contents.replace('\n','')
    cursor.execute(output_contents)
    set2 = cursor.fetchall()

    
    for row1 in set1:
	sequence_id_1,rank_1,orientation_1,left_end_1,right_end_1, sequence_length_1, taxon_id_1 = row1[0],row1[1],row1[2],int(row1[3]),int(row1[4]), row1[5], row1[6]
	accession_id_1 = sequence_id_1.split('|')[0]
	
	if not (accession_id_1 == 'NC_000913'):	#Ecoli ref gets counted as well. This is to prevent it. 
	    for row2 in set2:
		sequence_id_2,rank_2,orientation_2,left_end_2,right_end_2, sequence_length_2 = row2[0],row2[1],row2[2],int(row2[3]),int(row2[4]), row2[5]
		accession_id_2 = sequence_id_2.split('|')[0]
		if (accession_id_1 == accession_id_2):
		    number_tried += 1	#This keeps a count of all the number of pairwise differences tried. 
		    tried_gene_pairs += sequence_id_1 + ':' + sequence_id_2 + ':' + taxon_id_1 + ','
		    rank_diff = rank_2 - rank_1
		    if (abs(rank_diff) == 1):
			count+=1
			#Distance between the genes
			if rank_diff == 1: distance = left_end_2 - right_end_1
			if rank_diff == -1: distance = left_end_1 - right_end_2
			orthologous_gene_pairs += sequence_id_1 + '(' + str(sequence_length_1) + ')' +':' + sequence_id_2 + '(' + str(sequence_length_2) + ')' + ':' + str(distance) + ':' + taxon_id_1 + ',' 
			
    
    if not orthologous_gene_pairs: orthologous_gene_pairs = 'Null'
    if not tried_gene_pairs: tried_gene_pairs = 'Null'
    
    
    
    return count, orthologous_gene_pairs, number_tried, tried_gene_pairs
    

def operate_two_ecoli_rows(row_a,row_b):
    """Two ecoli rows are passed - one set1 and another set2. These contain the ranks of the operon_name,ecoli_bnumber,rank, group_id. 
    The operon_names are compared. If they are same, then the rank difference between the two are calculated. They must be one else it is passed again. 
    The walk is defined here as walk_rank1_rank2. 
    If the two criteria as satisfied then the two group_ids are taken and passed to mysql execution to return back the count of the number of times, the rank difference was indeed one. 
    """
    
    #Initializing
    frequency_count = 0
    orthologous_gene_pairs = 'Null'
    ecoli_walk = ''
    number_tried = 0
    tried_gene_pairs = 'Null'
        
    operon_name_a, ecol_sequence_id_a, rank_a, group_name_a = row_a[0],row_a[1],int(row_a[2]),row_a[3]
    operon_name_b, ecol_sequence_id_b, rank_b, group_name_b = row_b[0],row_b[1],int(row_b[2]),row_b[3]
    
    if ((rank_b - rank_a) == 1):
	frequency_count, orthologous_gene_pairs, number_tried, tried_gene_pairs = sql_execute_two_sets(group_name_a,group_name_b)
    else:
	frequency_count = 0
    
    ecoli_walk = operon_name_a + '|' + ecol_sequence_id_a + ':' + operon_name_b + '|' + ecol_sequence_id_b
    walk = ecoli_walk + '\t' + 'walk_%d_%d = '%(rank_a,rank_b) + '\t' + str(frequency_count) + '\t' + str(number_tried) + '\t' + orthologous_gene_pairs + '\t' + tried_gene_pairs
    
    print walk
    
    return walk 


def ecoli_data_file_parsing(file_name):
    """ This opens the file ecoli_data.txt and reads out two lines at a time and passes them to the function operate_two_ecoli_rows. The ecoli_data.txt is a file created in sqlite by using the following command - 
    ***
    SELECT ecoli_gene_operonset.operon_name, c.sequence_id, c.rank, c.gene_id, c.orientation, c.left_end, c.right_end, groups.group_id
    FROM 
    (SELECT sequence_id, rank, locus_tag, gene_id, orientation, left_end, right_end
    FROM coding_sequences
    WHERE accession_number = 'NC_000913') c
    INNER JOIN ecoli_gene_operonset
    ON c.locus_tag = ecoli_gene_operonset.b_number
    INNER JOIN groups
    ON c.sequence_id = groups.sequence_id
    ***
    The output was written to the file - ecoli_data.txt
    """
    
    # Initializing
    operon_name = locus_tag = rank = gene_id = orientation = left_end = right_end = group_id = ''
    ecoli_list = []

    #Reading file_name and obtaining a list
    ifile = open_file(file_name)
    
    lines = ifile.readlines()
    ifile.close()
    
    for line in lines:
	row = line.split('|')
	operon_name = row[0]
	sequence_id = row[1] + '|' + row[2]
	rank = row[3]
	group_id = row[8]
	ecoli_list.append([operon_name, sequence_id, rank, group_id])

    return ecoli_list


def walk_two_genes(ecoli_list):
    """ This function takes in the data generated from ecoli_list file and sends two rows to operate on
    """
    
    # Initializing
    i = 1
    all_walks = []
    
    # 
    row_a = ecoli_list[0]
    
    for row_b in ecoli_list:
	if i>1:
	    operon_name_a, ecol_sequence_id_a, rank_a, group_name_a = row_a[0],row_a[1],int(row_a[2]),row_a[3]
	    operon_name_b, ecol_sequence_id_b, rank_b, group_name_b = row_b[0],row_b[1],int(row_b[2]),row_b[3]
	    if operon_name_a == operon_name_b:
		walk = operate_two_ecoli_rows(row_a,row_b)
		all_walks.append(walk)
	    row_a = row_b
	i += 1

    return all_walks
    


if __name__ == '__main__':
    
    # Calling database using MySQLdb module
    conn = sqlite3.connect('orthomcl.db')
    cursor = conn.cursor()
    
    file_name = 'ecoli_data.txt'
    ecoli_list = ecoli_data_file_parsing(file_name)
    all_walks = walk_two_genes(ecoli_list)
        
    ofile = open_file('allwalks_output.txt','w')
    for walk in all_walks:
	ofile.write(walk+'\n')
	
    ofile.close()
    
    from timeit import Timer
    t = Timer(stmt = "walk_two_genes()",setup = "from __main__ import walk_two_genes",)
    print t.timeit(0)
    print strftime("%d %b %Y %H:%M:%S",localtime())
        
    cursor.close()
    conn.commit()
    conn.close()