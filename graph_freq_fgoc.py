#! usr/bin/python

# Script makes a graph for the frequency distribution of the fgoc scores calculated

import numpy as np
import pylab as p
import matplotlib
import sys
import os

def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened : %s "%(name_file,err.strerror)
	sys.exit(0)

    return file


def get_fgoc_list(ifile):
    """
    @param ifile: File_handle for the file containing the details of the walks. 
    @function: Retrieves all the fgoc scores from the file. 
    @return: list of all the fgoc scores
    """
    
    #Initializing
    lines = ifile.readlines()
    ifile.close()
    
    fgoc_list = [float((item.split('\t')[3])) for item in lines if not item.startswith('///')]
    
    return fgoc_list




def main():
    
    #(1) Open file - 'walk_in_bsub_operons.txt'
    file_name = 'walk_in_bsub_operons.txt'
    ifile = open_file(file_name)
    
    #(2) Obtain list of all fgoc scores
    fgoc_list = get_fgoc_list(ifile)
    print fgoc_list.count(0)
    print len(fgoc_list)
    
    #(3) Draw histogram
    
    bins = 500
    fig = p.figure(figsize=(15,15))
    ax1 = fig.add_subplot(1,1,1)
    n,bins,patches = p.hist(fgoc_list,bins,color='g')
    p.xlabel('fraction of conserved gene_pair to the total number of combined gene occurrences in B.subtilis')
    p.ylabel('frequency')
    p.title('Histogram of conservation of gene pairs of within B.subtilis operons ')
    fig_dir = os.getcwd() + '/figures/'
    fig_name = fig_dir + 'bsub_fgoc_within_operon.png'
    
    p.savefig(fig_name,format = "png", orientation='landscape')
    






if __name__ == '__main__':
    main()
    
    