#! usr/bin/python

# Script parses the xml file of dbtbs.xml. It makes an operon database -  gene_name  : operon_name


from xml.dom import minidom
import sys
import sqlite3


def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened : %s "%(name_file,err.strerror)
	sys.exit(0)

    return file


def open_xml(file_name):
    """
    @param file_name: This is the name of the xml file that will be opened by minidom
    @return: The parsed output of the xmlparser
    """
    
    try:
	xmldoc = minidom.parse(file_name)
    except IOError,err:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened : %s "%(file_name,err.strerror)
	sys.exit(0)
    
    return xmldoc


def get_gene_operons(operon):
    """
    @param operon: This is the element node of one operon. 
    @function: Parse this to generate a gene:operon list
    @return gene_operon: Returns gene_operon as a small list, which is then flattened. 
    @return operon_name: Returns the name of the operon
    """
    
    #Initializing
    gene_operon = []
    
    [op_name] = operon.getElementsByTagName('name')
    operon_name = op_name.firstChild.data
    
    [gen_names] = operon.getElementsByTagName('genes')
    gene_names = gen_names.firstChild.data
    genes = gene_names.split(',')
    for gene in genes:
	gene_operon.append([gene,operon_name])
    
    return gene_operon, operon_name


def parse_operon_list(operon_list):
    """
    @param operon_list: This is the list of all the child nodes with operon 
    @function: Retrieves the gene names and the operon name
    @return: Gene_name: operon_name (membership in the operon)
    """
    
    #Initializing
    all_gene_operons = []
    all_operons = []
    
    # (1) Run a loop through the operon_list 
    for operon in operon_list:
	gene_operon, operon_name = get_gene_operons(operon)
	number_genes_operon = len(gene_operon)
	all_operons.append([operon_name, number_genes_operon])
	
	for gene_operon_pair in gene_operon:	#This is to break the list and append the gene_operon_pair into a longer list
	    all_gene_operons.append(gene_operon_pair)
		    

    return all_gene_operons, all_operons

def parse_bsub_gbk():
    """
    @function: It parses the file - 'NC_000964.gbk' in the database directory. It retrieves all the cds information and returns it. 
    """
    
    from parse_gbk import Gbk_file
    
    path = '/home/jaggu/databases/All Bacterial Genomes/Bacillus_subtilis/NC_000964.gbk'
    bsub_gbk = Gbk_file(path)
    information = bsub_gbk.organism_information()
    cds_information = bsub_gbk.cds_information()
    
    return cds_information

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

def table_gene_operon(gene_operon):
    """
    @param gene_operon: A list containing gene_name:operon_name
    @function: Creates a table name - bsub_gene_operonset. Has two columns (a) gene_name (b) operon_name
    """
    table_name = 'bsub_gene_operonset'
    contents = ''' CREATE TABLE bsub_gene_operonset
		(
		gene_name VARCHAR (10),
		operon_name VARCHAR (100)
		)
		'''
    operate_table(table_name,contents,'create')
    
    for [gene,operon] in gene_operon:
	insert_table = ''' INSERT INTO bsub_gene_operonset
			(gene_name, operon_name)
			VALUES 
			('%s','%s')
			'''%(gene, operon)
	operate_table(table_name,insert_table)
	
    return True


def table_bsub_cds(cds):
    """
    @param cds: This is the list of coding_sequence information. These have to be written to a table
    @function: Writes the information in the coding_sequence_information into table name bsub_coding_sequences
    """
    
    table_name = 'bsub_coding_sequences'
    
    contents =''' CREATE TABLE bsub_coding_sequences
	    (gene_id INT(15),
	    locus_tag VARCHAR(15),
	    gene_name VARCHAR (10),
	    protein_id VARCHAR(15),
	    sequence_id VARCHAR(30),
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
    operate_table(table_name,contents,'create')
    
    for info in cds:
	[gene_id, locus_tag, gene_name, protein_id, sequence_id, rank, joined, orientation, left_end, right_end, sequence, sequence_length, accession_number, taxon_id] = info

	rank = int(rank[5:])
	
	insert_table ="""
	INSERT INTO %s
	(gene_id, locus_tag, gene_name, protein_id, sequence_id, rank, joined, orientation, left_end, right_end, sequence, sequence_length, accession_number, taxon_id)
	VALUES
	('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')
	"""%(table_name, gene_id, locus_tag, gene_name, protein_id, sequence_id, rank, joined, orientation, left_end, right_end, sequence, sequence_length, accession_number, taxon_id)
	operate_table(table_name,insert_table)

    return True
    
def main():
    """
    @function: Parses the dbtbs.xml file
    """
    
    #(1) Open file - dbtbs.xml: This is the xml file of bacillus subtilis containing all operons and the genes involved in them
    directory = '/project/marcotte/jagannath/projectfiles/Downloads/bsubtilis_database/'
    filename = directory + 'dbtbs.xml'
    xmldoc = open_xml(filename)
    
    #(2) Get list of all child nodes with operon as the tag
    operon_list = xmldoc.getElementsByTagName('operon')
    
    #(3) Parse the operon list and return back a list containing gene name, operon name; Another list is operon_name: number of genes 
    gene_operon_list, operon_number_genes_list = parse_operon_list(operon_list)
    
    #(4) Parse file - NC_000964.gbk; This is the genbank file in database directory under Bacillus subtilis. It contains all the identifiers including gene name
    cds_information = parse_bsub_gbk()
    
    #(5) Open a new sqlite3 db - bsub.db; 
    
    #(6) Create and update tables - 
    # (a) bsub_gene_operonset = gene_name : operon pair 
    table_gene_operon(gene_operon_list)
    # (b) bsub_coding_sequences
    table_bsub_cds(cds_information)
    
    return


if __name__ == '__main__':
    # Calling database using sqlite module
    conn = sqlite3.connect('bsub.db')
    cursor = conn.cursor()
    
    main()
    
    cursor.close()
    conn.commit()
    conn.close()
    
    import time
    print "Script - bsub_operon_parsing.py %s \t Completed \t %s"%(main, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))