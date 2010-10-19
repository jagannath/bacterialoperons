# Python File for Parsing gbk files

# Script takes in the database files in the gbk format. It uses the path suggested by the listfile. It parses the .gbk file into appropriate format and returns the parsed files in a separate directory. 
# The script also connects to the localhost and uses the already existing database organisms. It populates the table with the all the relevant details on the geneid, protein id, position info, rank, 


from pprint import pprint
import os
import sys
import MySQLdb
from warnings import filterwarnings

# Filtering out the warnings from mysql
filterwarnings('ignore',category = MySQLdb.Warning)

i = 0;
allcds_info = []
global acc_no
database = 'orthomcl266'
user_name = 'jaggusql'
password = 'ishani'

# Calling database using MySQLdb module
conn = MySQLdb.connect(host = "localhost",
			user = user_name,
			passwd = password,
			db = database)
cursor = conn.cursor()

# Creating a MySQL Table in the database - organisms. The table name will be the acc_no (organisms) and will have the same amount of information as that in the CDS
# Creating a table in the database organisms. The name of the table is the acc_no. It will carry the columns - 
# (a) Gene_id = varchar(15)  (b) position_gene = varchar(50) (c) genename = varchar(10) (d) protein_id = varchar(15)
# (e) codonstart = varchar(15) (f) rank = varchar(10) (g) sequence = varchar(1000) (i) Locus_tag - This is the unique identifier. (j) taxon_id and (k) sequence_id
cursor.execute(""" DROP TABLE IF EXISTS organisms """)
table_create = """ CREATE TABLE organisms
(gene_id INT(15),
locus_tag VARCHAR(15),
sequence_id VARCHAR(15),
protein_id VARCHAR(15),
codonstart VARCHAR(15),
rank INT,
sequence LONGTEXT,
acc_no VARCHAR(15),
taxon_id VARCHAR (15),
joined ENUM('yes','no'),
orientation ENUM ('forward','reverse'),
left_end INT(10),
right_end INT (10)
)
"""
cursor.execute(table_create)


def make_dir(dir_name,path = os.getcwd()):
    """This function makes a directory in the path passed to the function. If no path is passed then it takes the current working directory as default. If the directory already exists then it doesnt do anything and returns back the directory path. In both cases the new directory path is returned back """
    
    #Initializing
    new_dir_path = ''
    
    new_dir_path = path + '/' + dir_name
    try:
	os.mkdir(new_dir_path)
    except OSError:
	print "Directory already exists. Overwriting on files in directory if any!"
	pass
    
    return new_dir_path
    


def open_gbkfile(filename,taxon_number):
    """This function opens a gbk file and passes on the file handle - 'file'"""
    file = open(filename)
    print filename
    return [organism_ids(file),cd_ids(file,taxon_number)]

def cd_ids(src,taxon_number):
    """This looks for the geneinformation under CDS - gi(a number)|position|/locustag|/genename|/proteinid|/product or gene annotation|/translation=sequence"""
    # Declaration Statements
    count = 0
    allcds_info = []

    line = src.readline()
  
    while not line.startswith('//'):
	if line.startswith('     CDS'):
	    onecds_info = '>'
	    count+=1
	    flag = True
	    while flag:
		onecds_info+=line
		line = src.readline()
		if line.startswith('     gene') or line.startswith('ORIGIN') or line.startswith('     misc_feature') or line.startswith('     repeat_region') or line.startswith('     sig_peptide') or line.startswith('     mat_peptide'):
		    allcds_info.append([onecds_info])
		    flag = False
	else:
	    line = src.readline()
    
    return parse_cds_info(allcds_info,taxon_number)

def organism_ids(src):
    """This is to obtain information about the organism - LOCUS and DEFINITION. Here I am not making use of any while loop as I know that the first two lines carry all the needed information """
    flag = True
    line = src.readline()
    organism = 'Taxonomy : '
    if line.startswith('LOCUS'):
	locus_info = line[12:].split(12*' ')
	global acc_no
	acc_no = locus_info[0]
	line = src.readline()
    if line.startswith('DEFINITION'):
	strain = line[12:]
    while flag:
	line = src.readline()
	if line.startswith('  ORGANISM'):	#This is to get the taxonomy data of the organism. The section can be of varied lines but next section is REFERENCE
	    in_organism = True
	    flag = False
	    while in_organism:
		organism+= line+';'
		line = src.readline()
		if line.startswith('REFERENCE'):
		    in_organism=False
	organism = organism.replace('\n','')
	organism = organism.replace(' ','')
	organism = organism.replace(';;',';')
	organism = organism.replace('ORGANISM','')
    organism_information = acc_no + '|' + strain[:-1] + '|' + organism + '\n'
    return organism_information


def parse_position_gene(position_gene):
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


def parse_cds_info(allcds,taxon_number):
    """This is to parse the CDS information and put in fasta format"""
    all_cds_parsed = []
    number_cds = len(allcds)
    table_name = acc_no

    for i in range(0,number_cds):
	[cds] = allcds[i]                 # Now every cds is a string, amenable to manipulation 
	cds = cds.replace('  /','//')	# Some / in the middle of the annotation flanked by some character will not be identified
	cds = cds.replace(' ','')	#Removes all the empty spaces
	cds = cds.replace('\n','')	#Removes all the lines
	cds = cds.split('//')
	# Now I am going to extract just the needed strings, manipulate them and finally stitch them together 
	for j in range(0,len(cds)):
	    # Initialising the identifiers as None so that everything will have a value
#	    gene_id = position_gene = genename = protein_id = codonstart = sequence = 'Nothing'
	    if (cds[j].find('db_xref="GeneID')>-1):	#if there is a mismatch -1 is given back not False or 0! This is Gene ID and not GI 
		gene_id = cds[j][16:-1]
	    if ((cds[j].find('CDS') and cds[j].find('..'))>-1):
		position_gene = cds[j][4:]
	    if (cds[j].find('locus_tag')>-1):	#This is annoying as some of the gbk files have something called old_locus_tag
		if ((cds[j].find('old_locus_tag'))>-1): pass
		else: 
		    locus_tag = cds[j][11:-1]
            if (cds[j].find('gene=')>-1):
		genename = cds[j][6:-1]
	    else: genename = 'None'
	    if (cds[j].find('protein_id')>-1):
		protein_id = cds[j][12:-1]
	    if (cds[j].find('codon_start')>-1):
		codonstart = cds[j]
	    if (cds[j].find('translation=')>-1):
		sequence = cds[j][13:-1]
	rank = "rank={0}".format(i+1)
	
	joined, orientation, left_end, right_end = parse_position_gene(position_gene)	#Function which will parse the position_gene string well.
# The Annotation for every sequence is - 
#	>gi|gene_id|ref|proteinid|(rank)--pos=(Number)--acc_no\nSequence
	parsed = '>gid|' + gene_id + '|ref|' + protein_id + '|' +  rank + '--' + 'pos=' + position_gene + '--' + acc_no + '|locus|' + locus_tag + '|position=' + 'joined--' + joined + '--' + orientation + '--' + left_end + '--' + right_end +  '\n' + sequence
	all_cds_parsed.append(parsed)
	# Inserting the information of each row in the database 
	# There is another information needed - SEQUENCE_ID_X which is a four letter Taxon_id|gene_id. I can get the taxon_id from sifting through a list that I made between the accno and the taxon id. Will verify the accno matches that in the list and then that is the taxon_id and hence the sequence id. I am not keeping any exceptions here. 
	rank = int(rank[5:])

	#for row in taxon_id_list:
	    #if (row[0] == acc_no): taxon_id = row[1]
	taxon_id = acc_no
	sequence_id = taxon_id + '|' + locus_tag

	insert_table ="""
	INSERT INTO organisms
	(gene_id, locus_tag,sequence_id, protein_id, codonstart, rank, sequence,acc_no,taxon_id,joined, orientation, left_end,right_end)
	VALUES
	('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')
	"""%(gene_id,locus_tag,sequence_id, protein_id, codonstart, rank, sequence,acc_no, taxon_id, joined, orientation, left_end, right_end)
	cursor.execute(insert_table)
	
    return all_cds_parsed

def newly_written_file(organism_information, all_cds_information):
    """This function is to make a new file with the accession_number.f containing the Organism information in the first line and Fasta formatted all CDS """
    new_file_name = organism_information.split('|')[0] + '.faa'
    directory_path ='ParsedgbkFiles'
    dir_path = make_dir(directory_path)
    new_dir_path = make_dir(new_file_name[:-4],dir_path)
    print new_dir_path
    new_file_name = new_dir_path + '/' + new_file_name
    file = open(new_file_name, 'w')
#    file.write(organism_information)
    for i in range(0,(len(all_cds_information)-1)):
	file.write(all_cds_information[i]+'\n')
    file.close()
    return new_file_name
    
#making_newfiles(file_name)
def open_listfile(listfile):
    """This function simply opens up the all_filenames.lst file where the path for all the gbk files were available"""
    #Initializing
    i = 0
    
    all_newfiles_list = []
    file = open(listfile)
    lines = file.readlines()
    for path in lines:
	newfilepath = importing_parsing_gbk(path,str(i))
	i+=1
    # newfile with the entire path was created. It would be useful to make a file having all the list written
	all_newfiles_list.append(newfilepath)
    file.close()

    # File created - newfileslist.lst for having the path of all the newly parsed files. This will be extremely helpful for further accessing these files
    file = open('newfileslist_faa.lst','w')
    for newfilename in all_newfiles_list:
	file.write(newfilename+'\n')
    file.close()
    return "Complete"

def importing_parsing_gbk(file_name,taxon_number):
    """This imports the path sent from function - open_listfile and sends the filename to be opened by open_gbkfile. Returns nothing back"""
    file_name = file_name[:-1]		#This is to remove the \n at the end of evry line
    [organism_information, all_cds_information] = open_gbkfile(file_name,taxon_number)
    number_genes = len(all_cds_information)
    newfile = newly_written_file(organism_information, all_cds_information)
    print "File Created - ", newfile

    if (os.path.isfile('destination_list.lst') and taxon_number == 0):
	os.remove('destination_list.lst')

    afile = open('destination_list.lst','a')
    afile.write(taxon_number+'|'+ newfile + '\n')
    afile.close()
    return newfile


if os.path.isfile('destination_list.lst'):
    os.remove('destination_list.lst')


listfile = "list_266_1.lst"
check = open_listfile(listfile)
print "Completed"

cursor.close()
conn.commit()
conn.close()

    
