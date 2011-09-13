#! usr/bin/python

# This is a much faster and more cleaner way to observe the distribution of FGOC scores between and within operons. Using downloaded operon prediction text file from MicrobesOnline to determine which genes to consider within operon and not. A consequetive pairwise genes within operon is all in an operon and if there is a false then that gene pair is thrown outside operon. I will need to do a little editing on the operon prediction tool as if on reverse strand I need to exchange the gene pair order. I will then make a histograms of FGOC score between and within operon gene pairs. I will also normalise the distribution and calculate an odds ratio (FGOC (within): FGOC (between)) and I will then get a graph


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

# *** Opening Classes *** #
class Organism():
    """ This class uses the organism as query and can calculate the FGOC score for gene pairs within operon, between operon and even a defined gene pair queried. """
    def __init__(self, organism):
	self.organism = organism
	self.withinOperonPairs, self.betweenOperonPairs = self.getAllGenePairs(operonLines)
    def getWithinOperonPairs(self):
	return self.withinOperonPairs
    def getBetweenOperonPairs(self):
	return self.betweenOperonPairs
    def getdirFGOC(self,genePair):
	""" Uses the genePair and retrieves the Kegg Orthogroups and the corresponding dirFGOC score from ALL.MAG.dictionary """
	keggPair = [self.getKeggInfo(genePair[0]), self.getKeggInfo(genePair[1])]
	if keggPair[0] == 'UNKNOWN' or keggPair[1] == 'UNKNOWN':return None
	else: return dirFGOCPair_dict[str(keggPair)]
    def getFGOC(self, genePair):
	keggPair = [self.getKeggInfo(genePair[0]), self.getKeggInfo(genePair[1])]
	if keggPair[0] == 'UNKNOWN' or keggPair[1] == 'UNKNOWN':return None
	else: 
	    try: fgoc = fgocPair_dict[str(keggPair)]
	    except KeyError: fgoc = 'Pair not there'
	    return fgoc

    def getKeggInfo(self, locusTag):
	""" @function: Check if there is a kegg group associated with the locusTag. If not return Unknown Kegg, else the Kegg group """
	try:
	    keggGroup = locustag_group_dict[locusTag] 
	except KeyError:
	    keggGroup = 'UNKNOWN'
	return keggGroup
	
    def getAllGenePairs(self,lines):
	""" Parses the operon lines and gets a list of pairs - [within operon gene pairs], [between operon gene pairs]. The criteria of reversing the gene pairs based on orientation is considered too """
	withinOperonPairs = []
	betweenOperonPairs = []
	for line in lines:
	    geneA, geneB, operonStatus = (line.split('\t')[i] for i in [2,3,6])
	    try:
		orientation_a, orientation_b = locus_cds_dict[geneA][2], locus_cds_dict[geneB][2]
		orientation_status = self.get_orientation_status(orientation_a,orientation_b)
		if operonStatus == 'TRUE':# It is classified as operon. Use orientation information to reverse the pairs
		    if orientation_status == 'reverse':withinOperonPairs.append([geneB, geneA]) # Reversing direction of the pairs
		    elif orientation_status == 'forward': withinOperonPairs.append([geneA,geneB])
		    else: pass # This should never happen. but just in case!
		else:
		    if orientation_status == 'forward': betweenOperonPairs.append([geneA, geneB])
		    if orientation_status == 'reverse': betweenOperonPairs.append([geneB, geneA])
		    if orientation_status == 'opposite':
			if orientation_a == 'forward': betweenOperonPairs.append([geneA, geneB])
			if orientation_a == 'reverse': betweenOperonPairs.append([geneB, geneA])
	    except KeyError:
		pass
	return withinOperonPairs, betweenOperonPairs
    def drawHistogram(self, numberList, figHandle, **kwargs):
	""" Given a list and titles draw a histogram using numpy """
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
	givenbins = [x*0.02 for x in range(0,52)]
	n, bins, patches= ax.hist(numberList,bins=givenbins,color=color)
	ax.set_title(title)
	ax.set_ylabel(ytitle)
	ax.set_xlabel(xtitle)
	ax.set_xlim(0,1)
	ax.set_ylim(0,75)
	if axis: ax.axis(axis)
	save_path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/figures'
	figure_name = title + '.png'
	if option == 'save': p.savefig(save_path_dir + '/' + figure_name)
	return [n, bins, figHandle]
    
    def drawPlot(self,xaxisRange, nbrList, figHandle, **kwargs):
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
	
    def get_orientation_status(self,orientation_a, orientation_b):
	""" function: Compares the orientation of the two genes. It can be forward, reverse or opposite only """
	orientation_status = ''
	if orientation_a == 'forward' and orientation_b == 'forward':
	    orientation_status = 'forward'
	elif orientation_a == 'reverse' and orientation_b == 'reverse':
	    orientation_status = 'reverse'
	else:
	    orientation_status = 'opposite'
	return orientation_status

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

def loadRelevantFiles(argDict):
    """ All relevant files loaded and passed on as either file path """
    loadedFileList = []
    # locustag : keggGroup dictionary
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
    # list of Lines of predicted operon file
    path_dir = '/project/marcotte/jagannath/projectfiles/Downloads/operons'
    file_name = argDict['org'] + '.operon'	# Organism name is in lower letters
    ifile = open(path_dir+'/'+ file_name)
    operonLines = ifile.readlines()[1:]# First line is the header
    ifile.close()
    loadedFileList.append(file_name)
    # Print the files loaded
    for file_name in loadedFileList: print " File loaded : %s "%(file_name)
    return [locustag_group_dict, locus_cds_dict, fgocPair_dict, dirFGOCPair_dict, operonLines]

def parseArgument(args):
    """ parses argument to obtain a key val dictionary """
    argDict = {'org':''}
    for arg in args:
	key,val = arg.rstrip().split('=')
	argDict[key] = val
    return argDict

def checkInputStatus(argDict, argument):
    """ Checks the input arguments. For organism checks if there is a file for org.operon in the operon directory"""
    operon_dir = '/project/marcotte/jagannath/projectfiles/Downloads/operons'
    status = True
    if not os.path.isfile(operon_dir + '/' + argDict['org'] + '.operon'):
	print "Organism Operon prediction file not found"
	status = False
    if not argDict['org']: 
	print "Organism not entered in argument"
	status = False
    if not status : sys.exit(1)
    else: return True


def removeInf(nbrList):
    """ There may be numbers in the list with Inf. Put that equal to the maximum val obtained in the list """
    newList = []
    for nbr in nbrList:
	if nbr == float('infinity'):
	    nbrList.sort(reverse=True)
	    withoutInf = [item for item in nbrList if not item == float('infinity')]
	    newList.append(max(withoutInf))
	else: newList.append(nbr)
    return newList

def main():
    """ argDict and loadedList are in memory and can be accessed anywhere """
    query_org_id = argDict['org']
    
    # Opening class Organism
    org = Organism(query_org_id)
    betweenOperonPairs = org.getBetweenOperonPairs()
    withinOperonPairs = org.getWithinOperonPairs()
    allBetweenFGOC = [org.getFGOC(pair) for pair in betweenOperonPairs]
    allWithinFGOC = [org.getFGOC(pair) for pair in withinOperonPairs]
    betweenFGOC = [fgoc for fgoc in allBetweenFGOC if not (fgoc == None or fgoc == 'Pair not there')]
    withinFGOC = [fgoc for fgoc in allWithinFGOC if not (fgoc == None or fgoc == 'Pair not there')]

    # Drawing figure
    p.rcParams.update({'font.size':14})
    fig = p.figure(figsize=(10,10))
    withinFreq, withinBins, fig = org.drawHistogram(withinFGOC,fig,pos=221, option='binFreq', xtitle='FGOC Score', ytitle='Frequency',title='Within Operon Gene pairs')
    betweenFreq, betweenBins, fig = org.drawHistogram(betweenFGOC,fig,pos=222, option='binFreq',xtitle='FGOC Score',title='Between Operon Gene pairs')
    # Normalising and calculating the odds ratio
    withinFreqNorm = [item/sum(withinFreq) for item in withinFreq]
    betweenFreqNorm = [item/sum(betweenFreq) for item in betweenFreq]
    oddsRatioList = [withinFreqNorm[i]/betweenFreqNorm[i] for i in xrange(0,len(betweenFreq))]
    oddsRatioList = removeInf(oddsRatioList)
    
    xaxisRange = betweenBins[:-1]
    fig = org.drawPlot(xaxisRange,oddsRatioList,fig,axis=[0,1,-2,10],pos=212,xtitle='FGOC Score',ytitle='Odds Ratio (Within/Between) operons')
    
    save_path_dir = '/project/marcotte/jagannath/projectfiles/bacterialoperons/figures'
    figure_name = argDict['org']+'.within&betweenOperonPairs.fgoc.png'
    p.savefig(save_path_dir + '/' + figure_name)
    p.show()

if __name__ == '__main__':
    
    argument = sys.argv[1:]
    #pdb.set_trace()
    print "Processing %s ..."%(argument)
    # Setting some global dictionaries
    argDict = parseArgument(argument)
    checkInputStatus(argDict,argument)
    loadedList = loadRelevantFiles(argDict)
    [locustag_group_dict, locus_cds_dict, fgocPair_dict, dirFGOCPair_dict, operonLines] = loadedList
    
    
    main()
    
    
    import time
    print "Script - pathwayToCytoscape.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))

