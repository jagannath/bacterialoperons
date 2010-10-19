#! usr/bin/env python
# Script file to parse the gbk files (a newer modified version of the older file parsinggbk.py. 
# It will make room to incorporate other identifiers and better creation and updating of table. 

import os
import sys
import MySQLdb
import sqlite3
from warnings import filterwarnings
import timeit
import re



class Gbk_file:
    """ This class aims to do the entire parsing of the gbk file. It is separated into functions - (1) Retrieve the organisms information (2) Retrieve each sequence_information (3) Populate a sequence_information coding_sequences table (4) Populate the organisms information table
    """
    
    def __init__(self, file_name):
	self.file_name = file_name
	file_name_split = file_name.split('/')
	self.organism_accession_number = file_name_split[len(file_name_split)-1]
	self.cwd = os.getcwd()

    def open_file(self,name_file, open_status = 'r'):
	""" This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

	#Opening and reading/writing the passed file_name """
	try:
	    file = open(name_file,open_status)
	except IOError,err:	# If the file cannot be opened i am exiting
	    raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	    sys.exit(0)
    
	return file
    
    def make_dir(self,dir_name,path = os.getcwd()):
	"""This function makes a directory in the path passed to the function. If no path is passed then it takes the current working directory as default. If the directory already exists then it doesnt do anything and returns back the directory path. In both cases the new directory path is returned back """
	
	#Initializing
	new_dir_path = ''
	
	new_dir_path = path + '/' + dir_name
	try:
	    os.mkdir(new_dir_path)
	except OSError:
	    #print "Directory already exists. Overwriting on files in directory if any!"
	    pass
	    
	return new_dir_path
	    
    
    def parse_position_gene(self,position_gene):
	""" This function (created adhoc) will parse the position_gene (with four possibilities - (a) simple one - forward (b) simple one - reverse --complement(..) (c) join - forward --join(..) or (d) complement(join(..,..)) - reverse. It will return 4 variables (a) joined (yes or no) (b) orientation (forward or reverse) (c) left_end and (d) right_end. In forward : left_end < right_end and in reverse : left_end > right_end
	"""
	
	# Initializing
	joined = orientation = left_end = right_end = ''
	
	if position_gene.startswith('complement'):
	    orientation = 'reverse'
	    position_gene = position_gene.replace('complement','')
	    position_gene = (position_gene.replace('(','')).replace(')','')
	
	    if position_gene.startswith('join'):
		joined = 'yes'
		position_gene = position_gene.replace('join','')
		number_positions = position_gene.count(',')
		left_end = (position_gene.split('..'))[0]
		right_end = (position_gene.split('..'))[number_positions+1]
	    else:
		joined = 'no'
		left_end = (position_gene.split('..'))[0]
		right_end = (position_gene.split('..'))[1]
	else:
	    orientation = 'forward'
	    if position_gene.startswith('join'):
		joined = 'yes'
		position_gene = position_gene.replace('join','')
		position_gene = (position_gene.replace('(','')).replace(')','')
		number_positions = position_gene.count(',')
		left_end = (position_gene.split('..'))[0]
		right_end = (position_gene.split('..'))[number_positions+1]
	    else:
		joined = 'no'
		left_end = (position_gene.split('..'))[0]
		right_end = (position_gene.split('..'))[1]
	
	return joined, orientation, left_end, right_end
			    

    def parse_cds(self,all_cds):
	""" This function carries out the parsing of each cds to retrieve the relevant sequence information. 
	"""
	
	#Initializing
	i = 1
	all_cds_identifiers = []
	

	for cds in all_cds:
	    
	    cds = cds.replace('  /','//')	# Some / in the middle of the annotation flanked by some character will not be identified
	    cds = cds.replace(' ','')	#Removes all the empty spaces
	    cds = cds.replace('\n','')	#Removes all the lines
	    cds = cds.split('//')

	    for identifier in cds:
		if (identifier.find('db_xref="GeneID')>-1):	#if there is a mismatch -1 is given back not False or 0! This is Gene ID and not GI 
		    gene_id = identifier[16:-1]
		if ((identifier.startswith('CDS')) and re.match('..',identifier)):
		    position_gene = identifier[3:]
		if (identifier.find('locus_tag')>-1):	#This is annoying as some of the gbk files have something called old_locus_tag
		    if ((identifier.find('old_locus_tag'))>-1): pass
		    else: locus_tag = identifier[11:-1]
		if (identifier.find('gene=')>-1):
		    genename = identifier[6:-1]
		else: genename = 'None'
		if (identifier.find('protein_id')>-1):
		    protein_id = identifier[12:-1]
		if (identifier.find('codon_start')>-1):
		    codonstart = identifier
		if (identifier.find('translation=')>-1):
		    sequence = identifier[13:-1]
		    sequence_length = len(sequence)
		
	    rank = "rank=%d"%(i)

	    joined, orientation, left_end, right_end = self.parse_position_gene(position_gene)	#Function which will parse the position_gene string well.
	    i += 1
		
	    accession_number, taxon_id = self.organism_information[0], self.organism_information[3]
	    sequence_id = accession_number + '|' + locus_tag
	    
	    cds_identifiers = [gene_id, locus_tag, protein_id, sequence_id, rank, joined, orientation, left_end, right_end, sequence, sequence_length, accession_number, taxon_id]
	    all_cds_identifiers.append(cds_identifiers)
	
	return all_cds_identifiers
    
        
    def organism_information(self):
	""" This function finds out all the relevant information for the organism. These are - (1) Accession number or locus or Version - naming accession_number matching the filename. I think these are distinct but could be wrong. 
	(2) Size of the genome size_genome (3) organism_definition (Name of the organism) (4) taxonomic_order - what are all its family, class, order etc (5) taxon_id (this identifies the organism) 
	"""
	# Initializing
	line = organism_definition = taxonomic_order = taxon = taxon_id = ''
	organism_identifiers = []
	flag = None	# Default for flagging
	
	ifile = self.open_file(self.file_name)
	lines = ifile.readlines()
	ifile.close()
	
	print "Processing %s ..."%(self.organism_accession_number)
	
	for line in lines:
	    if line.startswith('LOCUS'):
		first_line = line.split(' ')
		for item in first_line: 
		    if len(item): 
			organism_identifiers.append(item)
					    
		accession_number, genome_size = organism_identifiers[1], organism_identifiers[2]
		
	    if line.startswith('DEFINITION'):
		organism_definition = line.replace('DEFINITION  ','')[:-1]
		    
	    if line.startswith('  ORGANISM'):# This triggers flag to capture the next few lines till it encounters a period (.) at the end of the description. 
		flag = True
	    if flag: taxonomic_order += line[:-1]
	    if (line.find('.')>-1): flag = None
		    
	    if line.startswith('                     /db_xref="taxon:'):
		taxon = line.split('"')
		taxon_id = (taxon[1])[6:]
			
	organism_definition = organism_definition.replace('(','_')
	organism_definition = organism_definition.replace(')','_')
	taxonomic_order = taxonomic_order.replace('(','_')
	taxonomic_order = taxonomic_order.replace(')','_')
	
	organism_information = [accession_number, organism_definition, taxonomic_order,taxon_id, genome_size]
	self.organism_information = organism_information
	
	return organism_information
	
    def cds_information(self):
	""" The function uses the filename that was passed to call this class and extracts information about all the CDS. This is the main parsing function and extracts the information of (1) gene id (2) locus_tag (3) protein_id (4) rank (5) orientation (6) left_position (7) right_position (8) Sequence (10) Sequence_length (11) Accession number (12) taxon_id 
	
	"""
	# Initializing
	flag = False
	one_cd_information = ''
	all_cds_information = []
	cds_identifiers = ''
	all_cds_identifiers = []
	
	accession_number = self.organism_information[0]
	taxon_id = self.organism_information[3]
	
	ifile = self.open_file(self.file_name)
	lines = ifile.readlines()
	ifile.close()
	
	# Extract the CDS alone
	for line in lines:
	        
	    if line.startswith('     CDS'): flag = True
	    
	    dont_if_start_list = ['     gene','ORIGIN','     misc_feature','     repeat_region','     mat_peptide','     sig_peptide']
	    for dont_start in dont_if_start_list:
		if line.startswith(dont_start):	
		    flag = False
		    if len(one_cd_information): all_cds_information.append(one_cd_information)
		    one_cd_information = ''
		    
	    if flag: one_cd_information += line
	    
	
	all_cds_information.pop(0)	#The first item is an empty list. Has to be popped out
	
	#Parsing each CDS for protein and sequence information
	
	self.all_cds_identifiers = self.parse_cds(all_cds_information)
	    
	return self.all_cds_identifiers

    def update_organisms_tables(self,table_name1, table_name2):
	""" This function updates the organisms table already created in the main function. It updates two tables - (1) coding_sequences table with the information of all the genes, their positions etc (2) organisms with the information about the organism, its accession number, taxon_id, taxon_details (and if later - the genome sequence)
	"""
	
	#Initializing
	all_cds_identifiers = self.all_cds_identifiers
	one_cds_identifier = []
	gene_id = locus_tag = protein_id = sequence_id = rank = joined = orientation = left_end = right_end = sequence = sequence_length = accession_number = taxon_id = ''
	organism_information = []
	accession_number = organism_definition = taxonomic_order = taxon_id = genome_size = ''
 	
	for one_cds_identifier in all_cds_identifiers:
	    [gene_id, locus_tag, protein_id, sequence_id, rank, joined, orientation, left_end, right_end, sequence, sequence_length, accession_number, taxon_id] = one_cds_identifier
	    rank = int(rank[5:])
	    # Inserting details of table_name1 : coding_sequences
	    insert_table ="""
	    INSERT INTO %s
	    (gene_id, locus_tag, protein_id, sequence_id, rank, joined, orientation, left_end, right_end, sequence, sequence_length, accession_number, taxon_id)
	    VALUES
	    ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')
	    """%(table_name1, gene_id,locus_tag,sequence_id, protein_id, rank, joined, orientation, left_end, right_end,sequence, sequence_length, accession_number, taxon_id )
	    cursor.execute(insert_table)
	    
	# Inserting details of table_name2 : organisms
	
	organism_information = self.organism_information
	[accession_number, organism_definition, taxonomic_order,taxon_id, genome_size] = organism_information
	
	insert_table ="""
	INSERT INTO %s
	(accession_number, organism_definition, taxonomic_order, taxon_id, genome_size)
	VALUES
	('%s','%s','%s','%s','%s')
	"""%(table_name2,accession_number, organism_definition, taxonomic_order, taxon_id, genome_size)
	try:
	    cursor.execute(insert_table)
	except:
	    print insert_table

	return "Completed"

    def write_file(self,file_type = 'regular'):
	""" This function writes a parsed file and can do it in two types - (a) Normal type where all the identifiers are written down as if in fasta format. 
	(b) Orthomcl complaint where there is a sequence_id that is written followed by the sequence
	"""
	# Initializing
	all_cds_identifiers = self.all_cds_identifiers
	cds_identifiers = ''
	[accession_number, organism_definition, taxonomic_order,taxon_id, genome_size] = self.organism_information
	gene_id = locus_tag = protein_id = sequence_id = rank = joined = orientation = left_end = right_end = sequence = sequence_length = ''
	write_status = 'Incomplete'
	
	if file_type == 'orthomcl_complaint':
	    #parameters_required = [3]	# These numbers are based on what parameters are needed for the first line
	    dir_path = self.make_dir('parsed_gbk_files')
	    new_dir_path = self.make_dir(accession_number,dir_path)
	    os.chdir(new_dir_path)
	    complaint_file_name = accession_number + '.fasta'
	    ofile = self.open_file(complaint_file_name,'w')
	    
	    #dir_path = self.make_dir('orthomcl_complaint_fasta')
	    #os.chdir(dir_path)
	    #complaint_file_name = accession_number + '.fasta'
	
	    for cds_identifiers in all_cds_identifiers:
		#parsed_information = '|'.join('%s'%(cds_identifiers[i]) for i in parameters_required) #Use this when more parameters need to be passed
		sequence_id = cds_identifiers[3]
		sequence = cds_identifiers[9]
		parsed_information = '>' + sequence_id + '\n' + sequence + '\n'
		ofile.write(parsed_information)
		
	    ofile.close()
	    destination_file_path = new_dir_path + '/' + complaint_file_name
	    os.chdir(self.cwd)
	    

	if file_type == 'regular':
	    parameters_required = [0,1,2,3,4,6,7,8,10,11]
	    dir_path = self.make_dir('parsed_gbk_files')
	    new_dir_path = self.make_dir(accession_number,dir_path)
	    os.chdir(new_dir_path)
	    parsed_fasta_file_name = accession_number + '.fasta'
	    ofile = self.open_file(parsed_fasta_file_name,'w')
	    
	    
	    for cds_identifiers in all_cds_identifiers:
		sequence = cds_identifiers[9]
		parsed_information = '>' + '|'.join('%s'%(cds_identifiers[i]) for i in parameters_required) + '\n' + sequence + '\n' 
		ofile.write(parsed_information)

	    ofile.close()
	    destination_file_path = new_dir_path + '/' + parsed_fasta_file_name
	    os.chdir(self.cwd)
	
	return destination_file_path
		    
    def number_cds(self):
	""" This function just gives back the number of coding_sequences present in the genome. Attempt to do this by counting the length of the array list - all_cds_identifiers
	"""
	
	length = len(self.all_cds_identifiers)
	
	return length
	

# END OF CLASS #

def open_list(listfile):
    """ This function opens the listfile which is the file containing the path of all the gbk files. It returns a list (i.e. all the file lines)
    """
    # Initializing
    paths = []
    path = ''
    all_paths = []
    
    ifile = open(listfile,'r')
    paths = ifile.readlines()
    ifile.close()
    
    for path in paths:
	all_paths.append(path[:-1])

    
    return all_paths


def operate_table(table_name,contents,operation='insert'):
    """Creating Table - The details of the table columns are passed on as the string with contents. The contents are passed as within the triple quote format
    """

    table_contents = contents
    if operation == 'create': 
	cursor.execute(""" DROP TABLE IF EXISTS %s """%table_name)	#This works well only when the warning sign is turned off
    else: 
	pass 

    cursor.execute(table_contents)

    return 



def create_tables(table1, table2):
    """ This function creates two tables - table 1 : coding_sequences and table 2: organisms 
    """
    
    contents1 =''' CREATE TABLE coding_sequences
	    (gene_id INT(15),
	    locus_tag VARCHAR(15),
	    protein_id VARCHAR(20),
	    sequence_id VARCHAR(15),
	    rank INT,	
	    joined VARCHAR(4),
	    orientation VARCHAR (7),
	    left_end INT(10),
	    right_end INT (10),
	    sequence LONGTEXT,
	    sequence_length MEDIUMINT,
	    accession_number VARCHAR(15),
	    taxon_id VARCHAR (15)
	    )
	    '''

    contents2 =''' CREATE TABLE organisms
	    (accession_number VARCHAR(15),
	    organism_definition VARCHAR (100),
	    taxonomic_order MEDIUMTEXT, 
	    taxon_id INT,
	    genome_size INT
	    )
	    '''
    
    operate_table(table1,contents1,'create')
    operate_table(table2,contents2,'create')

    return
    

def run_gbk():
    #Initializing
    destination_file_paths = []
    
    file_name = '/home/jaggu/databases/All Bacterial Genomes/Streptococcus_equi_zooepidemicus/NC_012470.gbk'
    table_name1 = 'coding_sequences'
    table_name2 = 'organisms'
    all_created_paths = open_list('list_266_1.lst')

    create_tables(table_name1, table_name2)
    
    for path in all_created_paths:
	gbk = Gbk_file(path)
    	information = gbk.organism_information()
	cds_information = gbk.cds_information()
	table_status = gbk.update_organisms_tables(table_name1,table_name2)
	number_cds = gbk.number_cds()
	#write_status = gbk.write_file(cds_information)
        destination_file_path = gbk.write_file('orthomcl_complaint')
        destination_file_paths.append(destination_file_path)
    
    return destination_file_paths
    
def separate_list(name_file,file_paths, number_parts = 8):
    """ This takes in a lst file and divides it into 8 parts and writes a file with names according to the name passed eg, dest_1.lst for filename passed - dest
    The default number of splits is 8, but can be modified. 
    """
    # Initializing
    i = 0
    big_list = []
    big_list += ['']*number_parts
    
    for path in file_paths:
	if i< number_parts:
	    big_list[i] += path + '\n'
	    i+=1
	    if (i == (number_parts-1)): i = 0
    
    for i in xrange(0,number_parts):
	file_name = name_file + '_%s.lst'%str(i)
	if (os.path.exists(file_name)):
	    os.remove(file_name)

    for i in xrange(0,number_parts):
	file_name = name_file + '_%s.lst'%str(i)
	try:
	    file = open(file_name,'a')
	except IOError,err:	# If the file cannot be opened i am exiting
	    raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	    sys.exit(0)
	file.write(big_list[i])
	file.close()
	
    return


if __name__ == '__main__':

    ## Filtering out the warnings from mysql
    #filterwarnings('ignore',category = MySQLdb.Warning)

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
    
    destination_file_paths = run_gbk()
    
    separate_list('created_fasta',destination_file_paths,8)
    
    from timeit import Timer
    t = Timer(stmt = "run_gbk()",setup = "from __main__ import run_gbk",)
    print t.timeit(0)

    cursor.close()
    conn.commit()
    conn.close()
