#! usr/bin/python

# Script to reduce the redundancy in the 16S rRNA

import os
import sys




def main(argument):
    """
    """
    
    all_orgs = []
    ifile = open('/project/marcotte/jagannath/projectfiles/Downloads/rRNA/prok_msa.fasta','r')
    lines = ifile.readlines()
    ifile.close()
    flag = True

    ofile = open('prok_reduced_16SrRNA.fasta','w')
    
    for line in lines:
	header = False
	if line.startswith('>'):
	    acc_nbr_with_dot = line.split(' ')[1]
	    acc_nbr = acc_nbr_with_dot.split('.')[0]
	    
	    if acc_nbr in all_orgs:
		flag = False
		pass
	    else:
		all_orgs.append(acc_nbr)
		ofile.write('>' + acc_nbr + '\n')
		flag = True
		header = True
	if flag:
	    if header:
		pass
	    else:
		ofile.write(line)
	
    return True

if __name__ == '__main__':
    
    main(sys.argv[1:])
    
    import time
    print "Script - operon_shuffling.py %s \t Completed \t %s"%(main, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    