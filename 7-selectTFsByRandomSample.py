#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 22:42:52 2017

@author: pitagoras
"""
import pandas as pd
import numpy as np
import os.path
from scipy import stats
from multiprocessing import Pool
import matplotlib.pyplot as plt

scriptpath = os.path.dirname(__file__)

genesPerTissueFolder = "../input/genesPerTissue/"
filteredBedIntersectPath = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
tfByGenesAndTissuesPath = "../results/tfByGenesAndTissues.tsv"
mostRepresentativeTFsPerTissuePath = "../results/mostRepresentativeTFsPerTissue.tsv"
leastRepresentativeTFsPerTissuePath = "../results/leastRepresentativeTFsPerTissue.tsv"
tissueNames = ["adipose_tissue", "adrenal_gland", "brain", "breast", "colon",
           "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary",
           "prostate", "skeletal_muscle", "testis", "thyriod"]

percent = 10.0/2.0
nSamples = 2000

#key: tissueName, value: set of genes representative to tissueName
genesPerTissue = dict()
#set of all geneNames
geneNames = set()
#set of tfNames
tfs = set()
#key: tfName, value: set of genes with at least one binding site to tfName
genesWithBSTo = dict()
#key: tissueName, value: list of random subsets of 'geneNames'
geneSamples = dict()
#key: (tfName, geneName), value: number of binding sites
tfbsCount = dict()

def fileRelativeToParentDir(file):
    return os.path.join(os.path.dirname(scriptpath), file)

def readGenesFromTissue(tissueName):
    if tissueName not in genesPerTissue:
        genesPerTissue[tissueName] = set()
        
    fileName = genesPerTissueFolder + "/" + tissueName + ".txt"
    f = open(fileName, 'r')
    
    for line in f:
        gene = line.replace('\n', '')
        genesPerTissue[tissueName].add(gene)
        
        geneNames.add(gene)
    f.close()

def readAllTissueGenes():
    for tissueName in tissueNames:
            readGenesFromTissue(tissueName)

def listGenesWithBS(row):
    if row['tfName'] in genesWithBSTo:
        genesWithBSTo[row['tfName']].add(row['geneName'])

def countTFBS(row):
    tfbsCount[(row['tfName'], row['geneName'])] = row['count']
                
def readTfAndGeneNames():
    filteredBedIntersect = pd.read_csv(filteredBedIntersectPath, sep="\t")
    print("\tListing genes")
    for geneName in filteredBedIntersect['geneName'].unique().tolist():
        geneNames.add(geneName)
    print("\tListing TFs")    
    for tfName in filteredBedIntersect['tfName'].unique().tolist():
        tfs.add(tfName)
        
    for tfName in tfs:
        genesWithBSTo[tfName] = set()
    print("\tListing genes with binding site to TFs")
    filteredBedIntersect.apply(lambda row: listGenesWithBS(row), axis=1)
    print("\tCounting BS for each TFxGene")
    filteredBedIntersect.apply(lambda row: countTFBS(row), axis=1)

def createGeneSamplesForTissue(tissueName,samplesAmount=nSamples):
    nGenesForTissue = len(genesPerTissue[tissueName])
    geneSamples[tissueName] = list()
    todo = samplesAmount
    while(todo > 0):
        geneSet = set()
        genesChosen = np.random.choice(list(geneNames), nGenesForTissue)
        for geneName in genesChosen:
            geneSet.add(geneName)
        geneSamples[tissueName].append(geneSet)
        todo -= 1
        
def createAllGeneSamples():
    for tissueName in tissueNames:
        createGeneSamplesForTissue(tissueName)
        print("\tSampled " + tissueName + " tissue!")
            
def countBSToTFinGeneSet(geneSet, tfName):
    total = 0
    for geneName in geneSet:
        if (tfName, geneName) in tfbsCount:
            total += tfbsCount[(tfName, geneName)]
            #print("tfbsCount[("+tfName+","+geneName+")]="+str(tfbsCount[(tfName, geneName)]))
    return total

def getBsCountInSamples(tfName, tissueName):
    genesWithBS = genesWithBSTo[tfName]
    bsCountInSamples = []
    for sample in geneSamples[tissueName]:
        geneSet = sample.intersection(genesWithBS)
        bsCount = countBSToTFinGeneSet(geneSet, tfName)
        bsCountInSamples.append(bsCount)
    return bsCountInSamples
            
def qualifyTFAtTissue(tfName, tissueName):
    genesWithBS = genesWithBSTo[tfName]
    genesInTissueWithBs = genesPerTissue[tissueName].intersection(genesWithBS)
    bsCountInTissue = countBSToTFinGeneSet(genesInTissueWithBs, tfName)
    bsCountInSamples = getBsCountInSamples(tfName, tissueName)
    bsCountInSamples.append(bsCountInTissue)
    
    npArray = np.asarray(bsCountInSamples)
    std = npArray.std()
    median = np.median(npArray)
    perc = stats.percentileofscore(bsCountInSamples, bsCountInTissue, 'weak')
    return (perc, median, std, bsCountInTissue, len(genesInTissueWithBs))

def makeRows(tfName):
    rowsCreated = []
    for ts in tissueNames:
        newRow = dict()
        newRow['tfName'] = tfName
        newRow['tissue'] = ts
        newRow['tissue_genes'] = len(genesPerTissue[ts])
        perc, median, std, bsCountInTissue, genesInTissueWithBs = qualifyTFAtTissue(tfName, ts)
        newRow['tissue_genes_with_bs'] = genesInTissueWithBs
        newRow['tissue_tfbs'] = bsCountInTissue
        newRow['mt_median'] = median
        newRow['mt_std'] = std
        newRow['mt_percentile'] = perc
        rowsCreated.append(newRow)
    return rowsCreated

def createMainDF():
    pool = Pool(processes = 5)
    rawRows = pool.map(makeRows, tfs)
    pool.close()
    pool.join()
    rows = []
    for rowSet in rawRows:
        rows += rowSet
    return pd.DataFrame(rows, columns=['tfName', 'tissue', 'tissue_genes', 
                                       'tissue_genes_with_bs', 'tissue_tfbs',
                                       'mt_median', 'mt_std', 'mt_percentile'])

def createDFs():
    tfAndTissueDF = createMainDF()
    
    tfAndTissueDF['is_outside_std'] = tfAndTissueDF.apply(lambda row: 
        (row['tissue_tfbs'] < row['mt_median'] - row['mt_std'])
        or (row['tissue_tfbs'] > row['mt_median'] + row['mt_std']), axis = 1)
    tfAndTissueDF = tfAndTissueDF[tfAndTissueDF.is_outside_std == True]
        
    topDF = tfAndTissueDF[tfAndTissueDF.mt_percentile >= (100.0 - percent)]
    topDF = topDF.sort_values(['mt_percentile'], ascending=False)
    topDF.to_csv(mostRepresentativeTFsPerTissuePath, sep='\t', index=False)
    
    bottomDF = tfAndTissueDF[tfAndTissueDF.mt_percentile <= percent]
    bottomDF = bottomDF.sort_values(['mt_percentile'], ascending=True)
    bottomDF.to_csv(leastRepresentativeTFsPerTissuePath, sep='\t', index=False)

def histForTFandTissue(tfName, tissue, nbins=120, samples=nSamples):
    createGeneSamplesForTissue(tissue,samplesAmount=samples)
    vals = getBsCountInSamples(tfName, tissue)
    n, bins, patches = plt.hist(vals,bins=nbins)
    plt.show()

def loadAll():
    print("Loading input")    
    readAllTissueGenes()
    readTfAndGeneNames()
    print("Creating samples")
    createAllGeneSamples()

loadAll()
createDFs()