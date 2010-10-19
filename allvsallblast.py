#!/usr/bin/env python

# Script to do a blast all vs all against the same organisms database
# 1. Take the parsed file (acno.faa) and blastall against a formatted database
# 2. The results is in CSV format and is parsed for the query name, hit name and evalue
# 3. The query name, hit name are converted to a long name with syntax - accno-gino-rank
# 4. The evalue is converted to log(10) evalue. 
# 5. A file is generated - as a base file (just as backup) where the results of blast are parsed and stored
# 6. Making the matrix file .mtx which is just a big matrix (rows - query genes; columns - database genes). A large number of entries are zeros and the other few 
#    that answered the blast hits have the -log10(eval).
# 7. Another file containing the list of all the CDS is also made. The format is the same as that in the CSV format file (blast_results)
# 8. A new directory 'Blast Results' is made in each query directory. It contains all the blast results file where it was the query ((a) the blast_results file 
#    (.bstprs) (b) the matrix file (just 0s and -log10eval) (c) File containing all the CDS parsed in appropriate format 


# Modules Imported
import sys
import os
from numpy import *
import csv
import math

def make_dir(dir_name,path = os.getcwd()):
    """This function makes a directory in the path passed to the function. If no path is passed then it takes the current working directory as default. If the directory already exists then it doesnt do anything and returns back the directory path. In both cases the new directory path is returned back """
    
    #Initializing
    new_dir_path = ''
    
    new_dir_path = path + '/' + dir_name
    try:
	os.mkdir(new_dir_path)
    except OSError:
	print "Directory %s already exists. Overwriting on files in directory if any!"%dir_name
	pass
    
    return new_dir_path

def open_file(file_name,open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file
    """
    
    #Opening and reading/writing the passed file_name """
    try:
	file = open(file_name,open_status)
    except IOError:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened. Exiting "%file_name
	sys.exit(0)

    return file
    
    
def number_CDS(file_name):
    """ This function counts the number of CDS in the bacterial genome. Note this is only from the parsed .faa file and not for gbk. It will count the number of '>' which indicates the CDS. Thus it should be able to give the number of genes. It returns back the number of genes. """
    
    #Initializing 
    number_genes = 0
    contents = ''
    
    file_handle = open_file(file_name)
    contents = file_handle.read()
    number_genes = contents.count('>')
    
    file_handle.close()
  
    return number_genes
    
def open_list_file(listfile):
    """ This function opens the listfile (which contains all the path for all the newly parsed files)
    It takes in the filename and passes back a list containing all the paths present in the lst file """
    
    # Initializing
    paths = []
    number_paths = 0
    
    file_handle = open_file(listfile)	
    lines = file_handle.readlines()
    for line in lines:
	paths.append(line[:-1])	# :-1 is to remove the \n newline character
	number_paths += 1
    
    file_handle.close()
    return number_paths, paths
    
def parsed_faa_file(file_name,path = os.getcwd()):
    """ This function opens the file passed to it (which is a .faa file) containing the sequences of all the CDS of the genome. The function will parse the file and create a list of all the CDS in the format generated by the csv_row_parsing function, i.e accno|gi|pos=*|rank. The Rank assigned is the critical parameter that will be used for setting the row and column numbers. It also takes in the path where you want the file to be stored. The default is the current working directory. It then generates a file .fprs containing this list of CDS."""
   
    #Initializing
    parsed_cds_list = []
    gi = ''
    pos = ''
    rank = ''
    accno = ''
    
    file_handle = open_file(file_name)
    lines = file_handle.readlines()
    
    for line in lines:	# The parsing is almost identical to that used in csv_row_parsing as well
	if line.startswith('>'):
	    gi = line.split('|')[1]	# Note this is a string although it is a number
	    pos = (line.split('|')[4]).split('--')[1]
	    rank = ((line.split('|')[4]).split('--')[0]).split('=')[1]
	    accno = (line.split('|')[4]).split('--')[2][:-1]	#[:-1] is to remove the endline character associated with the accno
	    egg = accno + '|' + gi + '|' + pos + '|' + rank
	    parsed_cds_list.append(egg)
    file_handle.close()
    
    # Writing the file in the path passed.The file name will be .fprs (for .faa file parsed)
    result_file_name = path + '/' + accno + '.fprs' 
    
    ofile = open_file(result_file_name,'w')
    for i in range(0,len(parsed_cds_list)):
	ofile.write(parsed_cds_list[i]+'\n')

    ofile.close()
    print "File Parsed (.faa) : %s"%result_file_name
    return 

def csv_row_parsing(row):
    """ This function parses the row passed from the CSV file that was generated by the blast results. It returns - (a) -log10(evalue) (b) query: accno|gi|pos=*|rank (c) database hit: accno|gi|pos=*|rank . What is returned is a (a) float (b) string and (c) string """
    
    # Initializing
    evalue_string = ''
    evalue_num = 0.0
    log_evalue = 0
    parsed_result = []	# This is the list that will be returned
    gi = ''
    pos = ''
    rank = ''
    accno = ''
    
    # Working with the evalue and converting it into -log10 (after all the string - float conversions). The log_evalue is actually -log10(eval) and returned as a float
    evalue_string = row[10]
    evalue_num = float(evalue_string.replace('e','E'))
    try:
	log_evalue = - (math.log10(evalue_num))
    except ValueError:
	log_evalue = 0
	pass
          
    # Working with the query and database attributes
    query_attr = row[0]
    database_attr  = row[1]
    
    for attr_set in [query_attr,database_attr]:
	gi = attr_set.split('|')[1]	# Note this is a string although it is a number
	pos = (attr_set.split('|')[4]).split('--')[1]
	rank = ((attr_set.split('|')[4]).split('--')[0]).split('=')[1]
	accno = (attr_set.split('|')[4]).split('--')[2]
	foo = accno + '|' + gi + '|' + pos + '|' + rank
	parsed_result.append(foo)

    parsed_result.append(log_evalue)
    return parsed_result
    
def parse_blast_results(result_file_name,num_query,num_database):
    """ This function does the parsing of the blast result which is in the CSV format. The parameters passed are the result_file_name(csv format with .bst), the number of CDS in query and in the database. 
    Plan to make a file in this syntax - 
    1. Line 1 and 2 tell the number of CDS in the query and the database respectively
    2. The resultant file (named result_file_name.bstprs) will have query:accno|gi|pos|rank	databse:accno|gi|pos|rank	log10(evalue)"""

    # Initializing
    parsed_row = ''
    row = ''
    
    print result_file_name
    # This is opening the blastresult_file file for reading
    ifile = open_file(result_file_name)
    
    # Naming a parsed_blast_file will be the result_file_name.bstprs. So will have to do a string manipulation to remove the file extension and substitute with the .bstprs extension. Using '.' to split
    parsed_blast_file = result_file_name.split('.')[0] + '.bstprs'
    ofile = open_file(parsed_blast_file,'w')
    
    # In the CSV file - querydetails, databasedetails and evalue are the 1st, 2nd and 11th column respectively
    # 1st column = gi|number|ref|number|rank=number--pos=number--organismaccno
    # 2nd column = same as above
    # 11th column = It is in format '*e-10'. So will replace the small e with the big E and then change the string to float. Take log10 of that convert it back to string
    csvread = csv.reader(ifile,delimiter = '\t')
    csvwrite = csv.writer(ofile,delimiter = '#')
    csvwrite.writerow(['Number of genes - Query: {0}'.format(num_query)])
    csvwrite.writerow(['Number of genes - Database: {0}'.format(num_database)])
    csvwrite.writerow(['Query','Database','-log10(Evalue)'])
    
    for row in csvread:
	parsed_row = csv_row_parsing(row)
	csvwrite.writerow(parsed_row)
	
    ifile.close()
    ofile.close()
    print "Blast Results Parsed : %s"%parsed_blast_file
    # The parsed blast result file is completed
    
    # Making the big matrix file which will mainly be zeros but for the eval of the blast hits, a very sparse matrix. I will be passing the lengths of the query and database as well
    # the outputfile handle (why again open the file) - ofile
    make_matrix_file(parsed_blast_file, num_query,num_database)
    

    return
    
def make_matrix_file(file_name,num_query,num_database):
    """ This is the function for making the elaborate matrix. I am going to use numpy to first create the zero matrix and then use the rank information of the hits in the parsed_csv_file to place the -log10eval. 
    """
    
    # Initializing
    mat = zeros((num_query,num_database))

    ifile = open_file(file_name)
    lines = ifile.readlines()

    try:
	for line in lines:
	    if (line.find('|')>-1):
		query_rank = int((line.split('#')[0]).split('|')[3])
		db_rank = int((line.split('#')[1]).split('|')[3])
		leval = float(line.split('#')[2])
		mat[query_rank,db_rank] = leval
    except IndexError:
	pass
    
    ifile.close()

    # Making the matrix file file_name.mtx
    
    # Changing the file name
    matrix_file_name = file_name[:-7] 	# Had to remove .bstprs from the filename
    
    save(matrix_file_name,mat)	# Saving as a binary file. Will need the load option to open it. 
    print "Matrix File Made : %s .npy" %matrix_file_name
    return
    
def blast_these_two(query,database,num_query,num_database):
    """ This function does the blastall of the query against the database. The two parameters passed are strings carrying the path to the respective directories of the organism. The directories already have the formatted database and the .faa file. This .faa file is the fasta files containing all the genes in the genome and is blastp against the database. The function takes in the query and database paths and returns the path of the results file which is in CSV format """
    
    # Initializing
    path_length = 0
    query_accno = ''
    database_accno = ''
    
    # Sometimes if I was to change paths of the files this set of commands will ensure i get the correct query and database_accno
    temp = query.split('/')
    path_length = len(temp)
    query_accno = temp[path_length-1][:-4]
    
    temp = database.split('/')
    path_length = len(temp)
    database_accno = temp[path_length-1][:-4]
    
    #defining a result_file_name
    # Going to make a directory called blastresults_csv in the query directory. And saving the result file in that
    query_dir = os.path.dirname(query)
    dir_name = 'blastresults_'+query_accno
    new_dir_path = make_dir(dir_name,query_dir)
    
    result_file_name = new_dir_path + '/' + query_accno + 'vs' + database_accno + '.bst'
    # Opening the os commands to run blastall. 
    evalue = '0.00001' # evalue (set at) = E-5
    command_line = 'blastall -p blastp -d ' + database + ' -i ' + query + ' -e ' + evalue + ' -m 8 ' + ' -o ' + result_file_name
    os.system(command_line)
    
    # Parsing the result_file.bst which is in csv format.
    parse_blast_results(result_file_name,num_query,num_database)

    # Calling function to create the file with the list of all the CDS parsed in appropriate format. You pass in the file name and the directory in which the file parsed must be stored
    parsed_faa_file(database,new_dir_path)
    

    return "Completed"


listfile = 'newfileslist_faa.lst'
number_paths, paths = open_list_file(listfile)

for i in range(0,(number_paths)):
    for j in range(0,(number_paths)):
	query = paths[i]
	database = paths[j]
	query_number_CDS = number_CDS(query)
	database_number_CDS = number_CDS(database)
	print "Number of genes - Query : %d"%query_number_CDS
	print "Number of genes - Database : %d"%database_number_CDS
	confirm = blast_these_two(query,database,query_number_CDS,database_number_CDS)
	print confirm

print "Entire operation completed"
