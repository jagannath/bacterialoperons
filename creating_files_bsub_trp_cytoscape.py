#! usr/bin/python

# Script for making file for cytoscape

import os
import sys












def main(argument):
    """
    """
    
    ifile = open(os.getcwd() + '/tryptophan_study/all_possible_gene_pairs.trp')
    lines = ifile.readlines()
    for line in lines:
	group_pair, fgoc = line.split('\t')[1],line.split('\t')[2]
	group_pair_eval = eval(group_pair)
	grp_a, grp_b = group_pair_eval[0],group_pair_eval[1]
	ofile1 = open('trp_bsub.sif','a')
	ofile1.write(grp_a + '\t' + '(ab)' + '\t' + grp_b + '\n')
	ofile2 = open('trp_bsub.attredges','a')
	ofile2.write(grp_a + ' ' + '(ab)' + ' ' + grp_b + ' ' + ' = ' + ' ' + fgoc + '\n')
	print fgoc

    ofile1.close()
    ofile2.close()



if __name__ == '__main__':
    # Calling database using sqlite module
    
    main(sys.argv[1:])
    
    import time
    print "Script - operon_shuffling.py %s \t Completed \t %s"%(main, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    