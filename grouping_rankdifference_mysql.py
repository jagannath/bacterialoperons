#!usr/bin/env python
# Creating mysql table groups - seqid | group_number from groups.txt.
# Planning to add the group_number to the organisms data table, probably by just column mapping

# Connecting to database - orthomcl
import os
import sys
import MySQLdb
from warnings import filterwarnings



def open_file(file_name,open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called
    """
    
    #Opening and reading/writing the passed file_name """
    try:
	file = open(file_name,open_status)
    except IOError:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened. Exiting "%file_name
	sys.exit(0)
	
    return file

def operate_table(table_name,contents,operation='insert'):
    """Creating Table - The details of the table columns are passed on as the string with contents. The contents are passed as within the triple quote format
    """
    
    table_contents = contents
    if operation == 'create': 
	cursor.execute(""" DROP TABLE IF EXISTS %s """%table_name)	#This works well only when the warning sign is turned off
    else: pass 
    
    cursor.execute(table_contents)
    
    return 
    

def parsing_groups(file_name):
    """ This function parses the groups.txt and makes rows containing a tuple of group_name and seqid of the organisms. There are a total of 5519 groups but a lot more seqids. This will then add the row to the table created. 
    """
    #Initializing
    
    
    #Creating table - Two Columns, Groupname (JEM_1...) and SeqID. 
    table_name = 'groups'
    contents = """ CREATE TABLE %s
		    (
		    sequence_id VARCHAR(30),
		    group_id VARCHAR(8)
		    )
		"""%(table_name)
    operate_table(table_name, contents, 'create')

    ifile = open_file(file_name)
    
    lines = ifile.readlines()
    for line in lines:
	split_line = line.split(' ')	#Space is all that separates the different sequenceids
	group_name = (split_line.pop(0))[:-1]	#It has a ':' in the end
	for sequences in split_line:
	    sequences = sequences.replace('\n','')	#Some seq id has this annoying carriage return
	    
	    #Inserting this information into the table
	    contents = """ INSERT INTO groups (sequence_id, group_id) VALUES ('%s','%s')
			"""%(sequences,group_name)
	    operate_table('groups',contents)

    return "Completed"

def operate_two_ecoli_rows(row_a,row_b):
    """Two ecoli rows are passed - one set1 and another set2. These contain the ranks of the operon_name,ecoli_bnumber,rank, group_id. 
    The operon_names are compared. If they are same, then the rank difference between the two are calculated. They must be one else it is passed again. 
    The walk is defined here as walk_rank1_rank2. 
    If the two criteria as satisfied then the two group_ids are taken and passed to mysql execution to return back the count of the number of times, the rank difference was indeed one. 
    """
    
    #Initializing
    frequency_count = 0
    
    operon_name_a, ecol_sequence_id_a, rank_a, group_name_a = row_a[0],row_a[1],row_a[2],row_a[3]
    operon_name_b, ecol_sequence_id_b, rank_b, group_name_b = row_b[0],row_b[1],row_b[2],row_b[3]
    
    if operon_name_a == operon_name_b:
	if ((rank_b - rank_a) == 1):
	    frequency_count = mysql_execute_two_sets(group_name_a,group_name_a)
	else:
	    frequency_count = 0
    else:
	frequency_count = 0

    walk = 'walk_{0}_{1} = '.format(rank_a,rank_b) + str(frequency_count)
    print walk
    return walk 

def mysql_execute_two_sets(group_name_1,group_name_2):
    """ This function executes the rank difference between two sets (rank difference between same organisms are only considered). If the rank difference is equal to absolute 1, then count increments by one. So the function returns back the frequency of the pairs in the two sets. It is pairwise difference between all the possible pairs between the two sets of the same organisms. This is because more than one gene of an organims can be in the same set. Returns back the count
    """
    #Initializing
    count = 0
    
    output_contents = """SELECT coding_sequences.sequence_id, coding_sequences.rank, groups.group_id
    FROM coding_sequences
    INNER JOIN groups
    ON coding_sequences.sequence_id = groups.sequence_id
    WHERE groups.group_id = '%s' """%group_name_1
    cursor.execute(output_contents)
    set1 = cursor.fetchall()
    
    output_contents = """SELECT coding_sequences.sequence_id, coding_sequences.rank, groups.group_id
    FROM coding_sequences
    INNER JOIN groups
    ON coding_sequences.sequence_id = groups.sequence_id
    WHERE groups.group_id = '%s'"""%group_name_2
    cursor.execute(output_contents)
    set2 = cursor.fetchall()
    
    for row1 in set1:
	sequence_id_1,rank_1 = row1[0],row1[1]
	taxon_id_1 = sequence_id_1.split('|')[0]
	for row2 in set2:
	    sequence_id_2,rank_2 = row2[0],row2[1]
	    taxon_id_2 = sequence_id_2.split('|')[0]
	    if (taxon_id_1 == taxon_id_2):
		rank_diff = rank_2 - rank_1
		if (abs(rank_diff) == 1):
		    count+=1
	
    return count


def operate_mysql_ecoli():
    """ This function operates via mysql on the ecoli_gene_operonset and organisms table to fetch back - (a) operon_name (b) seqid(ecoli genes) (c) Rank (d) Group_id
    The filter criteria is kept for the operon having more than one gene. So would get only about 2500ish genes. There are some E.coli genes which for some strange reason isnt available in the groups. So I am getting it back as NULL. I am pretty sure there are orthologs for that..Will keep it aside for the time. But note that it is breaking down the operon count to about 1900 and there are some even within the operon that arent turning up. Maybe with a larger database I may not have this problem! 
    """
    #Initializing
    allwalks = []
    
    contents = """ SELECT ecoli_gene_operonset.operon_name, coding_sequences.sequence_id, coding_sequences.rank, groups.group_id
		    FROM coding_sequences
		    INNER JOIN ecoli_gene_operonset
		    ON coding_sequences.locus_tag = ecoli_gene_operonset.b_number
		    INNER JOIN groups
		    ON coding_sequences.sequence_id = groups.sequence_id
		    WHERE ecoli_gene_operonset.number_genes_operon>1
		    ORDER BY coding_sequences.rank
		"""

    cursor.execute(contents)
    output = cursor.fetchall()
    ecoli_operon_1,ecoli_seq_id_1,ecoli_rank_1, ecoli_groupid_1 = output[0][0],output[0][1],output[0][2],output[0][3]
    row_0 = output[0]
    i = 1
    j = 1
    
    for row_1 in output:
	if i>1:
	    walk = operate_two_ecoli_rows(row_0,row_1)
	    allwalks.append(walk)
	    row_0 = row_1
	i+=1


    return allwalks

if __name__ == '__main__':

    # Calling database using MySQLdb module
    database = 'orthomcl266_1'
    user_name = 'jaggusql'
    password = 'ishani'
    
    conn = MySQLdb.connect(host = "localhost",
			user = user_name,
			passwd = password,
			db = database)
    cursor = conn.cursor()
    
    table_name = 'groups'
    #cursor.execute(""" DROP TABLE IF EXISTS %s """%table_name)	#This works well only when the warning sign is turned off

    #contents = """ CREATE TABLE %s
    #(
    #sequence_id VARCHAR(15),
    #group_id VARCHAR(8),
    #CONSTRAINT pk_sequence_id PRIMARY KEY (sequence_id)
    #)
    #"""%(table_name)
    #operate_table(table_name, contents, 'create')



    file_name = 'groups.txt'
    status = parsing_groups(file_name)
    print status

    allwalks = operate_mysql_ecoli()
    print allwalks

    ofile = open_file('allwalks_output.txt','w')
    for walk in allwalks:
	ofile.write(walk+'\n')

    ofile.close()
    cursor.close()
    conn.commit()
    conn.close()

