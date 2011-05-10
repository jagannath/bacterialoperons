#! usr/bin/python

# This script obtains the list of all operons in Ecoli genome and writes a command to 6 .sh files to be run.

import os

def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """
    
    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	raise AssertionError("File %s cannot be opened : %s "%(name_file,err.strerror))
	sys.exit(0)
    
    return file
    

if __name__ == '__main__':
    
    ifile = open_file('ecoli_homolog_data.txt')
    lines = ifile.readlines()
    ifile.close()
    
    all_operon_list =([line.split('\t')[0] for line in lines]) 
    operon_list = set(all_operon_list)
    print len(operon_list)
    count = 0
    i = 0
    j=0
    
    # For sh_genome_fgoc. Iteration through all organisms calling gen_fgoc_calc.py
    #path = os.getcwd() + '/parsed_gbk_files_genomes'
    #all_acc_nbrs = os.listdir(path)
   

    for operon in all_operon_list:
	#if all_operon_list.count(operon)>1 :
	j+=1
	i+=1
	dir_path = os.getcwd() + '/sh_files/'
	file_name = dir_path + 'ecoli_homolog_shuffle' + str(i) + '.sh'
	ofile = open_file(file_name,'a')
	ofile.write("python operon_shuffling.py %s shuffle \n"%(operon))
	if i == 6:
	    i=0
	
    print j
