#! usr/bin/python

from __future__ import division
import os
import sys
try:
   import cPickle as pickle
except:
   import pickle
import shelve

# This script parses the mclInput file to be inputed for mcl, making numpy matrix or cytoscape graph


def main(argument):
    [argv] = argument
    
    if argv == '1':
	# Parses the mclInput file and writes two files - edge attribute and sif format.
	file_name = 'mclInput_fgoc_ko_allvsall.txt'
	ifile = open(file_name)
	lines = (ifile.read()).splitlines()
	ifile.close()
	count = 0
	# Write the two files
	ofile1 = open('genevsgene.sif','w')
	ofile2 = open('genevsgene.eda','w')
	ofile3 = open('network_ecoli_allvsall_ko.cytoscape','w')
	for line in lines:
	    count+=1
	    print "%d %s "%(count,line)
	    split_len = len(line[:-1].split('\t'))
	    if split_len == 3:
		b1,b2,fgoc = line[:-1].split('\t')
		try:
		    assert len(b1) == 5
		    assert len(b2) == 5
		    fgoc = float(fgoc)
		    if not fgoc == 0:
			ofile1.write(b1 + ' bn ' + b2 + '\n')
			ofile2.write(b1 + ' (bn) ' + b2 + ' = ' + str(fgoc) +'\n')
			ofile3.write(b1 + ' ' + b2 + ' ' + str(fgoc) + '\n')
		except AssertionError:
		    pass

	ofile1.close()
	ofile2.close()
	
    if argv == '2':
	# Attribute file
	ifile = open('ecoli_cluster.txt')
	lines = (ifile.read()).splitlines()
	ifile.close()
	ofile = open('genevsgene.noa','w')
	#ofile.write('Gene Function'+'\n')
	for line in lines:
	    b_nbr, rank, func = line.split('\t')[1], line.split('\t')[3], line.split('\t')[8]
	    ofile.write(b_nbr + '\t' + rank + '|' + func + '\n')
	
	ofile.close()
    
    if argv == '3':
	# Making cluster attributes from the mcl output file so that I have cluster Number : Gene
	# (1) Opening the mcl output file
	file_name = 'mclOutput_ko_allvsall.mclout'
	ifile = open(file_name)
	lines = (ifile.read()).splitlines()
	ifile.close()
	count = 0
	cluster_name = ''
	ofile = open('kegg_all_ecoligenes.clust','w')
	# (2) Writing out as Cluster id bnumber
	for line in lines:
	    count+=1
	    gene_list = line.split('\t')
	    for gene in gene_list:
		info_towrite = gene + '\t' + cluster_name + str(count) + '\n'
		print info_towrite
		ofile.write(info_towrite)

    if argv == '4':
	# Want to rank the genes according to the maximum number of connectedness. 
	file_name = os.getcwd() + '/directed_fgoc_ecoli_ko/'+ 'directed_fgoc_ecoli_ko.mclinput'
	ifile = open(file_name)
	lines = ifile.read().splitlines()
	ifile.close()
	cutoff = 0.2
	file_name = os.getcwd() + '/cytoscape_files/' + 'directed_fgoc_ecoli_ko'+str(cutoff) + '.cytoscape'
	ofile = open(file_name,'w')
	for line in lines:
	    if line.startswith('b') and len(line.split('\t'))==8:
		#if float(line.split(' ')[2]) > cutoff:
		if float(line.split('\t')[3]) > cutoff:
		    ofile.write(line + '\n')
		    print line
	ofile.close()

    if argv == '5':
	# Using network_ecoli_allvsall_ko.cytoscape to pull up the ecoli pairs that must be considered and using that to calculate the newfgoc score, will call it the dirfgoc score (directed fgoc).
	file_name = 'network_ecoli_allvsall_ko.cytoscape'
	ifile = open(file_name)
	lines = ifile.read().splitlines()
	ifile.close()
	# Parsing input file to write to output in sh_files directory
	path_dir = os.getcwd() + '/sh_files/'
	index = 0
	for line in lines:
	    pair = line.split(' ')
	    index += 1
	    if index == 16: index = 1
	    
	    file_name = path_dir + 'calc_dir_fgoc_ecoli_ko_' + str(index) + '.sh'
	    info_line = 'python get_fgoc_pair.py ' + str(pair[0]) + ' ' + str(pair[1]) + '\n'
	    ofile = open(file_name,'a')
	    print info_line
	    ofile.write(info_line)
	    ofile.close()

    if argv == '6':
	file_name = os.getcwd() + '/directed_fgoc_ecoli_ko/'+ 'ecoli_kegg_dir_convdiv.fgoc'
	ifile = open(file_name)
	lines = ifile.read().splitlines()
	ifile.close()
	cutoff = 0.00
	file_name = os.getcwd() + '/cytoscape_files/' + 'ecoli_kegg_dir'+ '_withreverse_' + 'b1296_to_b1305_' + 'all_attributes' + '.cytoscape'
	ofile = open(file_name,'w')
	#geneset = ['b1268','b1267','b1266','b1265','b1264','b1263','b1262','b1261','b1260']
	geneset = ['b1296','b1297','b1298','b1299','b1301','b1302','b1303','b1304','b1305']
	for line in lines:
	    gene_a, gene_b, dir_fgoc,conv_fgoc, div_fgoc = line.split('\t')[0], line.split('\t')[1], line.split('\t')[3], line.split('\t')[4], line.split('\t')[5]
	    if gene_a in geneset:
		dir_fgoc = float(dir_fgoc)
		if dir_fgoc > cutoff:
		    info_towrite = gene_a + ' ' + 'dr' + ' ' + gene_b + ' ' + str(dir_fgoc) + '\n'
		    ofile.write(info_towrite)
		    info_towrite = gene_a + ' ' + 'cv' + ' ' + gene_b + ' ' + str(conv_fgoc) + '\n'
		    ofile.write(info_towrite)
		    info_towrite = gene_a + ' ' + 'dv' + ' ' + gene_b + ' ' + str(div_fgoc) + '\n'
		    ofile.write(info_towrite)
	ofile.close()


if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - make_genevsgene_mat.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    
