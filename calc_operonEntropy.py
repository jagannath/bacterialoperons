#! usr/bin/python

# This attempts to calculate the entropy associated with each (equally sized) operons and make them into those respective bins and probably do some kind of histograms for each. I will do only for Ecoli 

from __future__ import division
import os
import sys
try:
    import cPickle as pickle
except:
    import pickle
import shelve
import re
import random
import pdb
from operator import itemgetter
import numpy as np
import pylab as p
import matplotlib
import itertools
from matplotlib.patches import Polygon


class Operon():
    """ This class operon is used to compute the entropy associated with the gene pairs """
    def __init__(self,operonName):
	self.operonName = operonName
	self.geneList, self.orientation = operonName_geneList_dict[operonName]
	if self.orientation == 'reverse': self.geneList.reverse()	# If reverse then change the ordering scheme
    def actualOrder(self): return (' --> '.join(item for item in self.geneList) + ' | ' + self.orientation)
    def getPairs(self):
	geneList = self.geneList
	gene1 = geneList[0]
	genePairs = []
	i = 0
	for gene2 in geneList:
	    if i > 0:
		genePairs.append([gene1,gene2])
		gene1 = gene2
	    i+=1
	return genePairs
    def nbrGenes(self):	return len(self.geneList)
    def broken(self,genePair):
	""" Checks if the pair is broken by checking if each pair is different in rank by one"""
	brokenStatus = False
	geneA, geneB = locus_cds_dict[(genePair[0])], locus_cds_dict[(genePair[1])]
	if self.orientation == 'forward': rankDiff = int(geneB[1]) - int(geneA[1])
	if self.orientation == 'reverse': rankDiff = int(geneA[1]) - int(geneB[1])
	if not rankDiff == 1: brokenStatus = True
	return brokenStatus
    def permutate(self):
	return list(itertools.permutations(self.geneList,2))
    def hasAllKO(self):
	allKeggGroups = [self.getKeggInfo(gene) for gene in self.geneList]
	if 'UNKNOWN' in allKeggGroups: return False
	else: return True
    def getEntropy(self):
	""" For each pair convert to the string ([koA,koB]) and obtain the FGOC score. Obtain list of all Entropy and compute sum(plogp) """
	import math
	allFGOC = []
	entropy = 0
	permList = self.permutate()
	for (geneA,geneB) in permList: 
	    koA,koB = self.getKeggInfo(geneA), self.getKeggInfo(geneB)
	    try: fgoc = fgocPair_dict[str([koA,koB])]	# Since I have ensured that it is only Non broken and Kegg Containing operons
	    except KeyError: fgoc = 1
	    allFGOC.append(fgoc)
	for fgoc in allFGOC:
	    if fgoc == 0: fgoc = 1	# Just to remedy the math value error. The answer is still the same
	    entropy += fgoc * math.log(fgoc)
	return (-1 * entropy)
    def getKeggInfo(self, locusTag):
	""" @function: Check if there is a kegg group associated with the locusTag. If not return Unknown Kegg, else the Kegg group """
	try:
	    keggGroup = locustag_group_dict[locusTag] 
	except KeyError:
	    keggGroup = 'UNKNOWN'
	return keggGroup


class Graphics():
    
    def __init__(self,figHandle):
	self.figHandle = figHandle
    def drawHistogram(self, numberList, **kwargs):
	""" Given a list and titles draw a histogram using numpy """
	figHandle = self.figHandle
	color = 'g'
	xtitle = ytitle = title =  ''
	axis = []
	option = 'show'
	pos = 111
	font = 14
	for key, val in kwargs.items():
	    if key == 'color': color = val
	    if key == 'xtitle': xtitle = val
	    if key == 'ytitle': ytitle = val
	    if key == 'title': title = val
	    if key == 'option': option = val
	    if key == 'pos': pos = val
	    if key == 'axis': axis = axis
	ax = figHandle.add_subplot(pos)
	n, bins, patches= ax.hist(numberList,color=color)
	ax.set_title(title)
	ax.set_ylabel(ytitle)
	ax.set_xlabel(xtitle)
	if axis: ax.axis(axis)
	save_path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/figures'
	figure_name = title + '.png'
	if option == 'save': p.savefig(save_path_dir + '/' + figure_name)
	return [n, bins, figHandle]
    
    def drawPlot(self,xaxisRange, nbrList, **kwargs):
	figHandle=self.figHandle
	color = 'b'
	xtitle = ytitle = title =  ''
	axis = []
	option = 'show'
	pos = 111
	font = 14
	for key, val in kwargs.items():
	    if key == 'color': color = val
	    if key == 'xtitle': xtitle = val
	    if key == 'ytitle': ytitle = val
	    if key == 'title': title = val
	    if key == 'option': option = val
	    if key == 'pos': pos = val
	    if key == 'axis': axis = axis
	ax = figHandle.add_subplot(pos)
	ax.plot(xaxisRange, nbrList,color=color)
	ax.set_title(title)
	ax.set_xlabel(xtitle)
	ax.set_ylabel(ytitle)
	ax.axhline(y=1, linewidth=2,color='r')
	if axis: ax.axis(axis)
	return figHandle


def pickle_file(file_name, contents, path_dir=os.getcwd()):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    #pkl_dir = os.getcwd() + '/pkl_files/'
    pkl_dir = path_dir + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ofile = open(pkl_file,'wb')
    pickle.dump(contents, ofile)
    ofile.close()
    
    return True

def unpickle_file(file_name, path_dir=os.getcwd()):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    pkl_dir = path_dir + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ifile = open(pkl_file)
    contents = pickle.load(ifile)
    ifile.close()
    
    return contents



def parseArgument(args):
    """ parses argument to obtain a key val dictionary """
    argDict = {'org':''}
    for arg in args:
	key,val = arg.rstrip().split('=')
	argDict[key] = val
    return argDict

def checkInputStatus(argDict, argument):
    """ Checks the input arguments. For organism checks if there is a file detailing the org.operon in the operon directory"""
    operon_dir = '/project/marcotte/jagannath/projectfiles/Downloads/operons'
    status = True
    if not os.path.isfile(operon_dir + '/' + argDict['org'] + '.operon'):
	print "Organism Operon Annotated file not found"
	status = False
    if not argDict['org']: 
	print "Organism not entered in argument"
	status = False
    if not status : sys.exit(1)
    else: return True


def getOperon_dict(organism):
    """ This creates and operonName : [[geneList], orientation] """
    path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/txt_file_backup'
    fname = 'ecoli_homolog_data.txt'
    ifile = open(path_dir + '/' + fname)
    lines = ifile.readlines()
    ifile.close()
    operonName_geneInfo = {}
    operonExistList = []
    for line in lines:
	prevList = ''
	operonName, locusTag, orientation = (line.split('\t')[i] for i in [0,1,4])
	if operonName in operonExistList:
	    prevList = operonName_geneInfo[operonName]
	else:
	    operonName_geneInfo[operonName] = [[locusTag,orientation]]
	    operonExistList.append(operonName)
	if prevList:
	    prevList.append([locusTag,orientation])
	    operonName_geneInfo[operonName] = prevList
    
    operonName_geneList_dict = {}
    for operonName, geneList in operonName_geneInfo.items():
	allGenes = [item[0] for item in geneList]
	allOrientation = [item[1] for item in geneList]
	assert allOrientation.count(allOrientation[0]) == len(allOrientation)	# Just to ensure no extraneous opposite direction genes are present
	operonName_geneList_dict[operonName] = [allGenes,allOrientation[0]]
	
    return operonName_geneList_dict

def loadRelevantFiles(argDict):
    """ Load all relevant files and dictionaries """
    loadedFileList = []
    # locusTag: Kegg Group
    path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/graph_genomes/all_orgs_dictionaries'
    file_name = 'locustag_group_dict.pkl'
    locustag_group_dict = unpickle_file(file_name,path_dir)
    loadedFileList.append(file_name)
    file_name = 'ALL.pairs.FGOC_dict.pkl'
    fgocPair_dict = unpickle_file(file_name,  path_dir)
    loadedFileList.append(file_name)
    file_name = 'ALL.pairs.dirFGOC_dict.pkl'
    dirFGOCPair_dict = unpickle_file(file_name,  path_dir)
    loadedFileList.append(file_name)
    # Unshelve file containing locus_tag: detail
    shelve_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/shelve_files'
    file_name = 'locus:cds_information_dictionary'
    locus_cds_dict = shelve.open(shelve_dir + '/' + file_name,flag='r')
    loadedFileList.append(file_name)
    # Print Loaded Files
    for file_name in loadedFileList: print " File loaded : %s "%(file_name)
    return [locustag_group_dict, fgocPair_dict, dirFGOCPair_dict, locus_cds_dict]

def drawBoxPlot(data):
    #Plotting BoxPlots (for each of the List)
    p.rcParams.update({'font.size':15})
    fig = p.figure(figsize=(10,10))

    ax1 = fig.add_subplot(111)
    bp = p.boxplot(data,notch=1, vert=1, whis=1.5)
    p.setp(bp['boxes'],color='k')
    p.setp(bp['whiskers'],color='black')
    p.setp(bp['fliers'],color='green',marker='D',markersize=3.0)
    xtickNames = p.setp(ax1,xticklabels=['Three Genes \n \n N=83', 'Four Genes \n \n N=57','Five Genes \n \n N=32','Six Genes \n \n N=22','Seven Genes \n \n N=11','> Seven Genes \n \n N=16'])
    p.setp(xtickNames)
    p.ylabel('Entropy Of Operon per Gene Pair')
    p.title('Entropy Of E.coli operons')
    save_path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/figures'
    figure_name = argDict['org']+'.boxPlot.EntropyPerGene.png'
    p.savefig(save_path_dir + '/' + figure_name)
    p.show()
    
    return


def main():
    operonLength_entropy_dict = {'3':[],'4':[],'5':[],'6':[],'6plus':[]}
    list3 = []
    list4 = []
    list5 = []
    list6 = []
    list7 = []
    list7plus = []
    for name in operonName_geneList_dict.keys():
	operon = Operon(name)
	if operon.nbrGenes() > 2:
	    allGenePairs = operon.getPairs()
	    allBrokenStatus = []
	    nbrGenes = operon.nbrGenes()
	    for genePair in allGenePairs: allBrokenStatus.append(operon.broken(genePair))
	    if not (True in allBrokenStatus) and (operon.hasAllKO()): 
		entropy = operon.getEntropy()
		normEntropy = entropy/nbrGenes
		if nbrGenes == 3: list3.append([entropy, name, nbrGenes, normEntropy])
		elif nbrGenes == 4:list4.append([entropy, name, nbrGenes, normEntropy])
		elif nbrGenes == 5:list5.append([entropy, name, nbrGenes, normEntropy])
		elif nbrGenes == 6: list6.append([entropy, name, nbrGenes, normEntropy])
		elif nbrGenes == 7: list7.append([entropy, name, nbrGenes, normEntropy])
		else: list7plus.append([entropy, name, nbrGenes])

    operonLength_entropy_dict = {'3':list3, '4':list4,'5':list5,'6':list6,'7':list7, '7plus':list7plus}
    entropy3 = [item[0]/item[2] for item in operonLength_entropy_dict['3']]
    entropy4 = [item[0]/item[2] for item in operonLength_entropy_dict['4']]
    entropy5 = [item[0]/item[2] for item in operonLength_entropy_dict['5']]
    entropy6 = [item[0]/item[2] for item in operonLength_entropy_dict['6']]
    entropy7 = [item[0]/item[2] for item in operonLength_entropy_dict['7']]
    entropy7plus = [item[0]/item[2] for item in operonLength_entropy_dict['7plus']]

    data = [entropy3, entropy4, entropy5, entropy6, entropy7, entropy7plus]
    drawBoxPlot(data)
    
    sortEntropy = sorted(operonLength_entropy_dict['7plus'], key=itemgetter(0))
    bottomOperon = [item[1] for item in sortEntropy[0:5]]
    print bottomOperon
    topOperon = [item[1] for item in sortEntropy[len(sortEntropy)-5:]]
    print topOperon

if __name__ == '__main__':
    
    argument = sys.argv[1:]
    #pdb.set_trace()
    print "Processing %s ..."%(argument)
    # Setting some global dictionaries
    argDict = parseArgument(argument)
    checkInputStatus(argDict,argument)
    operonName_geneList_dict = getOperon_dict(argDict['org'])
    loadedList = loadRelevantFiles(argDict)
    [locustag_group_dict, fgocPair_dict, dirFGOCPair_dict, locus_cds_dict] = loadedList
    main()
    
    
