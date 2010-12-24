#! usr/bin/python

# Script is to parse the file - NC_000964.opr.door. This is actually a operon predicted file downloaded from website DOOR. It has a funny operon id instead of any name, but since it is prediction, i guess you need something like that. 

# Will be parsing the file and enter it in the sqlite database. 

import sqlite3
from warnings import filterwarnings
import sys


def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened : %s "%(name_file,err.strerror)
	sys.exit(0)

    return file

def get_operon_locus(ifile):
    """
    @param ifile: File handle for the file - NC_000964.opr.door
    @function: Parses the file and 
    @return: list of operon and locus_tag
    """
    #Initializing
    operon_locus = []
    
    lines = ifile.readlines()
    for line in lines:
	operon_id, locus_tag = line.split('\t')[0], line.split('\t')[2]
	operon_locus.append([operon_id,locus_tag])

    return operon_locus

def operate_table(table_name,contents,operation='insert'):
    """
    @function: Creating Table - The details of the table columns are passed on as the string with contents. The contents are passed as within the triple quote format
    """

    table_contents = contents
    
    if operation == 'create': 
	cursor.execute(""" DROP TABLE IF EXISTS %s """%table_name)	#This works well only when the warning sign is turned off
    else: 
	pass 

    cursor.execute(table_contents)

    return 


def table_operon_locus(operon_locus):
    """
    @param operon_locus: A list containing gene_name:operon_name
    @function: Creates a table name - bsub_locus_operonset_door. Has two columns (a) locus_tag (b) operon_name
    """
    table_name = 'bsub_locus_operonset_door'
    contents = ''' CREATE TABLE bsub_locus_operonset_door
		(
		locus_tag VARCHAR (10),
		operon_name VARCHAR (10)
		)
		'''
    operate_table(table_name,contents,'create')
    
    for [operon,locus_tag] in operon_locus:
	insert_table = ''' INSERT INTO bsub_locus_operonset_door
			(locus_tag, operon_name)
			VALUES 
			('%s','%s')
			'''%(locus_tag,operon)
	operate_table(table_name,insert_table)
	
    return True   
    
    
def main():
    """
    @function: Parses the file - NC_000964.opr.door_cut(first line is cut off) and creates a table = bsub_operon_locus
    """
       
    # (1) Opens file - NC_000964.opr.door
    directory = '/project/marcotte/jagannath/projectfiles/Downloads/bsubtilis_database/'
    filename = directory + 'NC_000964.opr.door_cut'
    ifile = open_file(filename)
    
    # (2) Get a list of all [operon_id locus_tag]
    operon_locus_list = get_operon_locus(ifile)
    ifile.close()
    
    # (3) Create and update table
    table_operon_locus(operon_locus_list)



if __name__ == '__main__':
    # Calling database using sqlite module
    conn = sqlite3.connect('bsub.db')
    cursor = conn.cursor()
    
    main()
    
    cursor.close()
    conn.commit()
    conn.close()
    
    import time
    print "Script - bsub_operon_parsing_door.py %s \t Completed \t %s"%(main, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))