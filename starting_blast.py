#! usr/bin/bash

# Adhoc programme to retrieve all the blast results that were generated when blastall stopped due to server going down. 
# So i am planning to read out all the written files in the blast_results directory. I will then compare the script lines in the bash file and get a list of all lines executed in the bash file. 

import os
import sys


def get_blast_result_files():
    """ This function gets the list of all the files in the blast_results directory and returns back the list along with the size of the file. 
    """
    
    # Initializing
    file_details = []
    path = os.getcwd() + '/blast_results'
    count = 0
    
    directory_list = os.listdir(path)

    for file_name in directory_list:
	query, database = file_name.split('vs')[0],file_name.split('vs')[1][:-6]
	file_path = path + '/' + file_name
	size = os.path.getsize(file_path)
	
	if size < 1000: 
	    count+=1
	else:
	    file_details.append(file_path)

	
    print count
    return file_details

def get_non_executed_bash(files,bash_file,result_file):
    """
    This function aims to check the blastall_# file to see which of the lines were executed.
    Lines that are not executed - i.e. arent in the get_files list is returned back
    command line: 
    blastall -p blastp -d /project/marcotte/jagannath/projectfiles/bacterialoperons/parsed_gbk_files_genomes/NC_009443/NC_009443.fasta -i /project/marcotte/jag
    annath/projectfiles/bacterialoperons/parsed_gbk_files_genomes/NC_009443/NC_009443.fasta -e 0.001 -m 8  -o /project/marcotte/jagannath/projectfiles/bacteria
    loperons/blast_results/NC_009443vsNC_009443.blast
    """
    # Initializing
    non_executed_lines = []
    
    ifile = open(bash_file,'r')
    lines = ifile.readlines()
    ifile.close()
    
    for line in lines:
	if line.startswith('blastall'):
	    split_line = line.split(' ')
	    path = split_line[13]
	    for file_name in files:
		if path == file_name:
		    pass
		else:
		    if line in non_executed_lines: pass
		    else:	
			print path
			non_executed_lines.append(line)
    
    ofile = open(result_file,'w')
    ofile.write(non_executed_lines)
    ofile.close()
    
    return non_executed_lines

def make_bash_files(lines):
    """ This function makes 10 bash files - blastall_bash_#.sh and divides the gigantic list into 10 parts
    """
    #initializing
    i = 0
    name_list = []
    
    for line in lines:
	name_list = 'list_%d'%(i)
	name_list.append(line)

	i+=1
    
    print name_list
    return
    
if __name__ == '__main__':
    all_lines = []
    i = 0
    argument = int(sys.argv[1])
    print "Argument = %d" %(argument)
    
    files = get_blast_result_files()
    bash_files = ['blastall_0.sh','blastall_1.sh','blastall_2.sh','blastall_3.sh','blastall_4.sh']
    
    bash_file = bash_files[argument]
    result_file = 'blastall_bash_%d.sh'%(argument)
    non_executed_lines = get_non_executed_bash(files,bash_file,result_file)
    
    import time
    print "Argument = %d"%(argument)
    print time.strftime("%d %b %Y %H:%M:%S",time.localtime())
    
    #make_bash_files(all_lines)
    