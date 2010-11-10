# This script is to make a formatdb and blast command for shell on all the Fasta files from the organisms with the complete genomes (not plasmids)

import os
import os.path
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
	


#def open_list_file(listfile):
    #""" This function opens the listfile (which contains all the path for all the newly parsed files)
    #It takes in the filename and passes back a list containing all the paths present in the lst file """
    
    ## Initializing
    #paths = []
    #number_paths = 0

    #file_handle = open(listfile,'r')	
    #lines = file_handle.readlines()
    #for line in lines:
	#paths.append(line[:-1])	# :-1 is to remove the \n newline character
	#number_paths += 1
	
    #file_handle.close()
    
    #return number_paths, paths

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

def blast_all_script(file_paths,number=0):
    
    #Initializing
    all_command_lines = []
    
    # Creating shell file and writing the firstline
    new_blast_scripts_dir = make_dir('blast_sh')
    file_name = new_blast_scripts_dir + '/' + 'blastall_%s.sh'%(str(number))
    ifile = open(file_name,'w')
    #ifile.write('#!/bin/bash \n')
    print "Writing shell scripts ..."
    for query in file_paths:
	for database in file_paths:
	    query_dir, fasta_file_query = os.path.split(query)
	    database_dir, fasta_file_database = os.path.split(database)
	    
	    #defining a result_file_name
	    # Going to make a directory called blastresults. And saving all result file in that
	    cwd = os.getcwd()
	    new_dir_path = make_dir('blast_results')
	    query_name = fasta_file_query[:-6]
	    database_name = fasta_file_database[:-6]
	    
	    
	    result_file_name = new_dir_path + '/' + query_name + 'vs' + database_name + '.blast'
	    #print result_file_name
	    evalue = '0.001' # evalue (set at) = E-5
	    command_line = 'blastall -p blastp -d ' + database + ' -i ' + query + ' -e ' + evalue + ' -m 8 ' + ' -o ' + result_file_name + '\n' + 'echo ' + result_file_name + ' > out.1 \n'
	    
	    all_command_lines.append(command_line)

    ifile.writelines(all_command_lines)
    ifile.close()

    return

def get_list_fasta():
    """ This performs an os.walk on the directory parsed_gbk_files_genomes and picks up all the .fasta files out
    """
    #Initializing
    fasta_file_list = []
    base_dir = os.getcwd() + '/parsed_gbk_files_genomes'	#Even when in untar the tar file this is the base directory
    ecoli_path = ''
    
    for root, dirs, files in os.walk(base_dir):
	for item in files:
	    if item.endswith('.fasta'):
		fasta_file = root + '/' + item
		fasta_file_list.append(fasta_file)
		if item == 'NC_000913.fasta': ecoli_path = fasta_file
		
    return len(fasta_file_list), fasta_file_list, ecoli_path

if __name__ == '__main__':
    
    cwd = os.getcwd()
    #list_file = 'created_fasta_%s.lst'%(str(number))
    number_paths, paths,ecoli_ref_path = get_list_fasta()
    #shell_format_database(paths,number)
    number_organisms = 200	
    shorter_list = paths[0:100]	#The first 200 organisms
    #Will just generate a single blastall_#.sh file which will contain all the combinations of the shorter_list
    # Check if it contains the ecoli_ref. If not appends it to the list
    if ecoli_ref_path in shorter_list:
	pass
    else:
	shorter_list.append(ecoli_ref_path)
    
    
    number = 1
    blast_all_script(shorter_list,number)
    
    print "completed"
	

