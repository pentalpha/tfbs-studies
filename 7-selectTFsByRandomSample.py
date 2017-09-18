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

scriptpath = os.path.dirname(__file__)

genesPerTissueFolder = "../input/genesPerTissue/"
filteredBedIntersectPath = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
tfByGenesAndTissuesPath = "../results/tfByGenesAndTissues.tsv"
mostRepresentativeTFsPerTissuePath = "../results/mostRepresentativeTFsPerTissue.tsv"
leastRepresentativeTFsPerTissuePath = "../results/leastRepresentativeTFsPerTissue.tsv"
tissueNames = ["adipose_tissue", "adrenal_gland", "brain", "breast", "colon",
           "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary",
           "prostate", "skeletal_muscle", "testis", "thyriod"]

percent = 5.0
nSamples = 2000

genesPerTissue = dict()
tissuePerGene = dict()
geneNames = set()
tfsInTissue = dict()
tfs = set()
genesWithBSTo = dict()
geneSamples = dict()
tfbs = dict()

filteredBedIntersect = pd.read_csv(filteredBedIntersectPath, sep="\t")

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
            
            if gene not in tissuePerGene:
                tissuePerGene[gene] = list()
            tissuePerGene[gene].append(tissueName)
            
            geneNames.add(gene)
        f.close()

def readAllTissueGenes():
        for tissueName in tissueNames:
                readGenesFromTissue(tissueName)
                
def addTFBScount(row, ts):
        tfbs[(row['tfName'], ts)] = row[ts]
                
def readTfAndGeneNames():
    for geneName in filteredBedIntersect['geneFullName'].unique().tolist():
        geneNames.add(geneName)
        
    tfByGenesAndTissues = pd.read_csv(tfByGenesAndTissuesPath, sep='\t')
    tfByGenesAndTissues = tfByGenesAndTissues[tfByGenesAndTissues.allTissues > 8]
    
    for ts in tissueNames:
        tfByGenesAndTissues.apply(lambda row: addTFBScount(row, ts), axis=1)
        
    for tfName in tfByGenesAndTissues['tfName'].unique().tolist():
        tfs.add(tfName)
        
    for tfName in tfs:
        genesWithBSTo[tfName] = set()

def createGeneSamplesForTissue(tissueName):
    nGenesForTissue = len(genesPerTissue[tissueName])
    geneSamples[tissueName] = list()
    todo = nSamples
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
        print("Sampled " + tissueName + " tissue!")
                        
def qualifyTFAtTissue(tfName, tissueName):
    genesWithBS = genesWithBSTo[tfName]
    genesInTissue = genesPerTissue[tissueName].intersection(genesWithBS)
    sampleIntersections = []
    for sample in geneSamples[tissueName]:
        nGenes = len(sample.intersection(genesWithBS))
        sampleIntersections.append(nGenes)
    res = stats.percentileofscore(sampleIntersections, len(genesInTissue), 'weak')
    return res

def qualifyTFForEachTissue(tfName):
    row = dict()
    row['tfName'] = tfName
    for tissueName in tissueNames:
        row[tissueName] = qualifyTFAtTissue(tfName, tissueName)
    return row

def getDFforAllTFs():
    print("<Creating rows>")
    pool = Pool(processes = 5)
    rows = pool.map(qualifyTFForEachTissue, tfs)
    pool.close()
    pool.join()
    print("</Created rows>")
    df = pd.DataFrame(rows, columns=(['tfName'] + tissueNames))
    return df

def separateDFInTuples(generalDF):
    rows = []
    def makeRows(originalRow):
        for ts in tissueNames:
            newRow = dict()
            tfName = originalRow['tfName']
            newRow['tfName'] = tfName
            newRow['tissue'] = ts
            if (originalRow['tfName'], ts) in tfbs:
                newRow['tfbs_count'] = tfbs[(tfName, ts)]
            else:
                newRow['tfbs_count'] = 0
            newRow['tfbs_count_all'] = len(genesWithBSTo[tfName])
            newRow['percentile'] = originalRow[ts]
            rows.append(newRow)
    generalDF.apply(lambda row: makeRows(row), axis=1)
    return pd.DataFrame(rows, columns=['tfName', 'tissue', 'tfbs_count',
                                       'tfbs_count_all', 'percentile'])

print("Loading input")    
readAllTissueGenes()
readTfAndGeneNames()
print("Creating samples")
createAllGeneSamples()
print("Reading genes per TF, this may take a while")

def countGenesWithBS(row):
    if row['tfName'] in genesWithBSTo:
        genesWithBSTo[row['tfName']].add(row['geneFullName'])    
filteredBedIntersect.apply(lambda row: countGenesWithBS(row), axis=1)

generalDF = getDFforAllTFs()
tfAndTissueDF = separateDFInTuples(generalDF)
topDF = tfAndTissueDF[tfAndTissueDF.tfbs_count > 0]
topDF = topDF[topDF.percentile >= (100.0 - percent)]
topDF = topDF.sort_values(['percentile'], ascending=False)
topDF.to_csv(mostRepresentativeTFsPerTissuePath, sep='\t', index=False)

bottomDF = tfAndTissueDF[tfAndTissueDF.percentile <= percent]
bottomDF = bottomDF.sort_values(['percentile'], ascending=True)
bottomDF.to_csv(leastRepresentativeTFsPerTissuePath, sep='\t', index=False)
