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
    # it is split('\t')[1] for between operons and [3] for within operons
    
    pairs_conserved = [(item.split('\t')[1], item.split('\t')[3]) for item in lines if not item.startswith('///') and float(item.split('\t')[3]) > 0.95]
    #print pairs_conserved
    # Writing to file
    #ofile = open_file('within_operons_conserved_pairs(homology).txt','w')
    #for item in pairs_conserved:
	#ofile.write(str(item)+'\n')
    #ofile.close()
    #print len(pairs_conserved)
    
    return fgoc_list


def main():
    
    #(1) Open file - 'walk_in_bsub_operons.txt'
    file_name = 'walk_in_ecoli_operons_ver6_homolog.txt'
    ifile = open_file(file_name)
    
    #(2) Obtain list of all fgoc scores
    fgoc_list = get_fgoc_list(ifile)
    cutofflist = [item for item in fgoc_list if item == 1]
    print len(cutofflist)
    print fgoc_list.count(0)
    print len(fgoc_list)
    
    
    
    
    
    #(3) Draw histogram
    
    #bins = 500
    #fig = p.figure(figsize=(15,15))
    #ax1 = fig.add_subplot(1,1,1)
    #n,bins,patches = p.hist(fgoc_list,bins,color='g')
    #p.xlabel('fraction of conserved gene_pair to the total number of combined gene occurrences in E.coli (homology)')
    #p.ylabel('frequency')
    #p.title('Histogram of conservation of gene pairs within E.coli operons ')
    #ax2 = p.axes([.5, .5, .3, .3])
    #p.title('Histogram excluding zero')
    ##ax2 = fig.add_subplot(0.65,0.6,0.2)
    #n,bins,patches = p.hist(cutofflist,bins,color='b')
    #p.setp(ax2)
    #fig_dir = os.getcwd() + '/figures/'
    #fig_name = fig_dir + 'ecoli_fgoc_within_operon(homology)_ver6.png'
    #p.savefig(fig_name,format = "png", orientation='landscape')
    ##p.show()






if __name__ == '__main__':
    main()
    
    