# This script is to make a formatdb on all the Fasta files of the 10 interesting organisms

import os
import re

def make_dir(dir_name,path = os.getcwd()):
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
	


def open_list_file(listfile):
    """ This function opens the listfile (which contains all the path for all the newly parsed files)
    It takes in the filename and passes back a list containing all the paths present in the lst file """
    
    # Initializing
    paths = []
    number_paths = 0

    file_handle = open(listfile,'r')	
    lines = file_handle.readlines()
    for line in lines:
	paths.append(line[:-1])	# :-1 is to remove the \n newline character
	number_paths += 1
	
    file_handle.close()
    
    return number_paths, paths

def shell_format_database(file_paths,number):
    """ This function makes a shell script for performing a formatdb on all the files in the directory in which the parsed gbk file is stored
    """
    
    # Initializing
    all_command_lines = []
    
    # Creating shell file and writing the firstline
    file_name = 'format_database_%s.sh'%(str(number))
    ifile = open(file_name,'w')
    ifile.write('#!/bin/bash \n')
    
    for path in file_paths:
	dir_path, fasta_file = os.path.split(path)
	log_file = dir_path + '/' + fasta_file[:-5] + 'log'
	input_file = path
	command_line = 'formatdb -i ' + input_file + ' -oF ' + '-pT -l ' + log_file + '\n'
	all_command_lines.append(command_line)
    
    ifile.writelines(all_command_lines)
    ifile.close()
    
    return

def blast_all(file_paths,number=0):
    
    #Initializing
    all_command_lines = []
    
    # Creating shell file and writing the firstline
    file_name = 'blastall_%s.sh'%(str(number))
    ifile = open(file_name,'w')
    ifile.write('#!/bin/bash \n')
    
    for query in file_paths:
	for database in file_paths:
	    query_dir, fasta_file_query = os.path.split(query)
	    database_dir, fasta_file_database = os.path.split(database)
	    
	    #defining a result_file_name
	    # Going to make a directory called blastresults_csv in the a new directory called blastresults. And saving all result file in that
	    cwd = os.getcwd()
	    new_dir_path = make_dir('blast_results')
	    print new_dir_path
	    query_name = fasta_file_query[:-6]
	    database_name = fasta_file_database[:-6]
	    
	    
	    result_file_name = new_dir_path + '/' + query_name + 'vs' + database_name + '.blast'
	    #print result_file_name
	    evalue = '0.00001' # evalue (set at) = E-5
	    command_line = 'blastall -p blastp -d ' + database + ' -i ' + query + ' -e ' + evalue + ' -m 8 ' + ' -o ' + result_file_name + '\n'
	    all_command_lines.append(command_line)

    ifile.writelines(all_command_lines)
    ifile.close()

    return


if __name__ == '__main__':
    
    cwd = os.getcwd()
    for number in xrange(0,4):
	list_file = 'created_fasta_%s.lst'%(str(number))
        number_paths, paths = open_list_file(list_file)
	shell_format_database(paths,number)
	blast_all(paths,number)
	
    print "completed"
	
