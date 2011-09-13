# Simple script to write all Ecoli Operon genes and Non Operon Genes

path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/txt_file_backup'
fname = 'operonset.txt'

ifile = open(path_dir + '/' + fname)
lines = ifile.readlines()
ifile.close()

fname = 'NotOperonGenes.ecoli'
ofile = open(fname,'w')
header = '# Operon Name' + '\t' + 'Nbr Genes' + '\t' + 'Name' + '\t' + 'LocusTag' + '\n'
ofile.write(header)
for line in lines:
    operonName, nbrGenes, geneList = [line.split('\t')[i] for i in [0,1,3]]
    if int(nbrGenes) == 1:
	geneName, geneLocusTag = geneList.split('|')
	if not geneLocusTag == ',':
	    ofile.write(operonName + '\t' + nbrGenes + '\t' + geneName + '\t' + geneLocusTag + '\n')
ofile.close()
    