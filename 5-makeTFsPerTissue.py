# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from multiprocessing import Pool
import os.path

scriptpath = os.path.dirname(__file__)

genesPerTissueFolder = "../input/genesPerTissue/"
filteredBedIntersect = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
tissueNames = ["adipose_tissue", "adrenal_gland", "brain", "breast", "colon",
           "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary",
           "prostate", "skeletal_muscle", "testis", "thyriod"]

print("</Reading input data>")
filteredBedIntersect = pd.read_csv(filteredBedIntersect, sep="\t")

genesPerTissue = dict()
tissuePerGene = dict()
geneNames = set()
tfsInTissue = dict()
tfs = set()

def fileRelativeToParentDir(file):
    return os.path.join(os.path.dirname(scriptpath), file)

def readGenesFromTissue(tissueName):
        if tissueName not in genesPerTissue:
            genesPerTissue[tissueName] = list()
            
        fileName = genesPerTissueFolder + "/" + tissueName + ".txt"
        f = open(fileName, 'r')
        
        for line in f:
            gene = line.replace('\n', '')
            genesPerTissue[tissueName].append(gene)
            
            if gene not in tissuePerGene:
                tissuePerGene[gene] = list()
            tissuePerGene[gene].append(tissueName)
            
            geneNames.add(gene)
        f.close()

def readAllTissueGenes():
        for tissueName in tissueNames:
                readGenesFromTissue(tissueName)
                
def addTfbsToTfCounts(row):
    geneName = row['geneFullName']
    tfName = row['tfName']
    tfs.add(tfName)
    for tissue in tissueNames:
        if geneName in genesPerTissue[tissue]:
            if tfName not in tfsInTissue[tissue]:
                tfsInTissue[tissue][tfName] = 1
            else:
                tfsInTissue[tissue][tfName] += 1

def writeTFlistForTissue(tissue):
    print("<Started listing for " + tissue + ">")
    tissueTFs = "results/TFsPerTissue/" + tissue + "_tfs.txt"
    with open(fileRelativeToParentDir(tissueTFs), 'w') as output:
        output.write("tfName\tfrequency\n")
        for tfName in tfsInTissue[tissue]:
            output.write(tfName + "\t" + str(tfsInTissue[tissue][tfName]) + "\n")
        output.close()
    print("</TF list saved at " + tissueTFs + ">")

print("Reading tissue genes lists")
readAllTissueGenes()
for tissueName in tissueNames:
        tfsInTissue[tissueName] = dict()

print("Looking for relations between genes and TFs in tissues\nThis can take a while.")
filteredBedIntersect.apply(lambda row: addTfbsToTfCounts(row), axis=1)

for tissueName in tissueNames:
        writeTFlistForTissue(tissueName)