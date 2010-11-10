#! /usr/bin/env python

# Script file to parse the Ecogene data (downloaded - 9th Sept 2010). I removed the first row (which contained the column description. The file is parsed and a
# mysql table 'ecogene_data' is generated. It contains 10 columns - 
#(1) EG number - Ecogene Number - eg_number VARCHAR(8)
#(2) Gene name - 3 - 4 letter - gene_name VARCHAR (5)
#(3) ECK number - 7 chars eck_number VARCHAR(8)
#(4) B number - Blattner number - b# - b_number (same as locus_tag) - Is the primary key - VARCHAR(6)
#(5) Length of gene - typically in 100s - length_gene SMALLINT(4) UNSIGNED
#(6) Orientation - Enumerate - forward/reverse orientation VARCHAR (7)
#(7) Left End - typically 7 numbers - left_end INT UNSIGNED
#(8) Right End - typically 7 numbers - right_end INT UNSIGNED
#(9) Mneumonic - Expanded form of gene name, provides some short description of gene - gene_mneumonic VARCHAR(1000)
#(10) Description - A good description of gene characteristics. - gene_description LONGTEXT

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
#conn = MySQLdb.connect(host = "localhost",
#user = user_name,
#passwd = password,
#db = database)
conn = sqlite3.connect('all_orgs.db')

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
    
def create_table(table_name):
    """Creating Table - The details of the table columns are given in the description at the beginning of the file
    """
	
    create_table = """ CREATE TABLE %s 
	(
	eg_number VARCHAR(8),
	gene_name VARCHAR(10),
	eck_number VARCHAR(8),
	b_number VARCHAR(5),
	length_gene SMALLINT UNSIGNED,
	orientation VARCHAR (8),
	left_end INT UNSIGNED,
	right_end INT UNSIGNED,
	gene_short_desc VARCHAR (500),
	gene_description LONGTEXT,
	CONSTRAINT pk_b_number PRIMARY KEY (b_number)
	)
	"""%(table_name)
    cursor.execute(""" DROP TABLE IF EXISTS ecogene_data """)	#This works well only when the warning sign is turned off
    cursor.execute(create_table)
	
    return 

def parse_ecogene_file(filename,table_name):
    """ This function parses the ecogene_data file. It is given the filename of the ecogene_file. The file is a csv file, tab delimited and without the header line. 
    # Note this important - the header row must be removed!! It contains 10 columns and they are in the order in which the script file described it.
    """
    
    # Initializing
    

    file_handle = open_file(filename)
    lines = file_handle.readlines()
    file_handle.close()
    
    for line in lines:
	#Initializing all the 10 variables here.
	b_number = eg_number = eck_number = length_gene = orientation = left_end = right_end = gene_short_desc = gene_description = ''
	
	row = line.split('\t')
	# There is some ambiguity in some b_numbers. There can be many. These mainly are Insertion elements. So I am removing those rows that contain more than one b_number. I see that they dont map onto the orthology dataset too as well very nicely in the operon set. The b_number must be unique so I put the ifstatement to check that out first. 
	if (row[3].count(';')==0):
	    b_number = row[3]

	    eg_number = row[0]
	    gene_name = row[1]
	    eck_number = row[2]
	    length_gene = row[4]
	
	    if (row[5].find('Counterclockwise') > -1):	# Changing orientation to forward and reverse. Doing Counterclockwise first to ensure no pattern matches first
		orientation = row[5].replace('Counterclockwise','reverse')
	    elif (row[5].find('Clockwise') > -1):
		orientation = row[5].replace('Clockwise','forward')
	    else:
		print "Orientation must be Clockwise or Counterclockwise only!"
		import sys
		sys.exit(0)

	    left_end = row[6]
	    right_end = row[7]
	    gene_short_desc = row[8]
	
	    gene_description = row[9]	
	    gene_description = gene_description.replace('\'','--')	#This is because some of the description carries a single quote which is a bug in mysql
	
	    # Populating the table - table_name (in this case is ecogene_data). 
	    insert_table ="""
	    INSERT INTO %s
	    (eg_number,gene_name,eck_number,b_number,length_gene,orientation,left_end,right_end,gene_short_desc,gene_description)
	    VALUES
	    ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')
	    """%(table_name,eg_number,gene_name,eck_number,b_number,length_gene,orientation,left_end,right_end,gene_short_desc,gene_description)
	    cursor.execute(insert_table)
	
    return
    


create_table('ecogene_data')
table_name = 'ecogene_data'
parse_ecogene_file('ecogenedata.txt',table_name)



cursor.close()
conn.commit()
conn.close()



