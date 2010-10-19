#!usr/bin/env python

# Script file to parse the Operonset data (downloaded - 9th Sept 2010). I removed the comment rows from the file which is in the current working directory. It tells all about citations and what each column means. It is going to be similar to ecogene_data parsing file. Although operon is more than one gene, the file doesnt say so. I will still be making a column with number of genes. Another point to note - this operonset file has many insertions as well which I encountered earlier. 
# mysql table 'gene_operonset' is generated. It contains 4 columns - 
#(1) b_number (VARCHAR 6) (2) gene_name (VARCHAR 5) (3) operon_name (VARCHAR 40)(4) Number of genes in that operon in which it is a part of (SMALLINT 1). 
# Another table with just some operon information - ecoli_operon_data (just in case i need the info about how it was arranged)
# (1) Operon_name (VARCHAR 40) (2) Number_genes (SMALLINT (1))(3) Orientation VARCHAR (7) (4) how_identified (MEDIUMTEXT)

# Connecting to database - orthomcl
import os
import sys
import sqlite3
import MySQLdb
from warnings import filterwarnings


# Calling database using sqlite3 module
# Calling database using MySQLdb module
database = 'orthomcl266_1'
user_name = 'jaggusql'
password = 'ishani'

# Calling database using MySQLdb module
conn = MySQLdb.connect(host = "localhost",
user = user_name,
passwd = password,
db = database)
cursor = conn.cursor()

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
    
def create_table(table_name,contents):
    """Creating Table - The details of the table columns are passed on as the string with contents. The contents are passed as within the triple quote format
    """

    create_table = contents
    cursor.execute(""" DROP TABLE IF EXISTS %s """%table_name)	#This works well only when the warning sign is turned off
    cursor.execute(create_table)
	
    return 

def parse_operon_file(file_name):
    """ This is the big function that does the parsing to result in two mysql tables
    """

    file_handle = open_file(file_name)
    lines = file_handle.readlines()
    file_handle.close()
    
    for line in lines:
	# Initializing
	gene_name = b_number = genes = number_bs = operon_name = number_genes_operon = ''
	
	row = line.split('\t')
	
	# Inserting the different splitted_rows into the table ecoli_operon_data. This will be a longer table with more operons. I had trimmed many in the unique gene_operonset
	
	insert_table_ecoli_operon = """
				    INSERT INTO ecoli_operon_data
				    (operon_name,number_genes_operon,orientation, how_identified)
				    VALUES
				    ('%s','%s','%s','%s')
				    """%(row[0],row[1],row[2],row[4])
	cursor.execute(insert_table_ecoli_operon)
	
	# I have to check if the operon has more than one gene. So another loop within this structure. I also have to create a table within this loop to ensure that all genes, operons etc are consistent.
	number_bs = row[3].count(',')
	genes = row[3].split(',')	# This works also when there is one gene too in the operon. 
	for i in xrange(0,number_bs):
	    gene_name, b_number = genes[i].split('|')
	    if (b_number.startswith('b')):	# I had to keep this constraint as some of the gene_name had first one ok, but second one without a b_number. This way it ensures that the b_numbers are always unique and each have an operon tag to it! 
		operon_name = row[0]
		number_genes_operon = row[1]
	    
		# Inserting this into the table. I also kept the gene_operonset fixed for the table name. Too many things to pass otherwise
		insert_table ="""
		INSERT INTO ecoli_gene_operonset
		(gene_name,b_number,operon_name,number_genes_operon)
		VALUES
		('%s','%s','%s','%s')
		"""%(gene_name,b_number,operon_name,number_genes_operon)
		cursor.execute(insert_table)

    return
    

gene_operonset_contents = """CREATE TABLE ecoli_gene_operonset
			    (
			    gene_name VARCHAR(10),
			    b_number VARCHAR(6),
			    operon_name VARCHAR(100),
			    number_genes_operon SMALLINT(1),
			    CONSTRAINT pk_b_number PRIMARY KEY (b_number)
			    )
			  """

operon_data_contents = """ CREATE TABLE ecoli_operon_data
			(
			operon_name VARCHAR (100),
			number_genes_operon SMALLINT(1),
			orientation VARCHAR (7),
			how_identified MEDIUMTEXT,
			CONSTRAINT pk_operon_name PRIMARY KEY (operon_name)
			)
			"""

create_table('ecoli_operon_data',operon_data_contents)
create_table('ecoli_gene_operonset',gene_operonset_contents)
file_name = 'operonset.txt'

parse_operon_file(file_name)

cursor.close()
conn.commit()
conn.close()



    
    
  
    