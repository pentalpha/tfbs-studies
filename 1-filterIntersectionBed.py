#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 15:32:58 2017

@author: pitagoras
"""
import pandas as pd
import numpy as np

geneNames = set()
bindingSiteFrequency = dict()

humanGenesFPKMInTissues = "../input/table_Human_body_map_ze_FPKM.txt"
bedIntersectPath = "../results/bedIntersectWaWbTFBSinGenes.bed"
filteredBedIntersectPath = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"

def readGeneNames():
    print("Reading gene names: ")
    df = pd.read_csv(humanGenesFPKMInTissues, sep="\t")
    df.apply(lambda row: geneNames.add(row["geneName"]), axis=1)
    print(str(len(geneNames)) + " names")

def decideBetwenGeneNames(row):
    names = row['geneFullName'].split(sep=',')
    for name in names:
        if name in geneNames:
            return name
    return ''

def addGeneAandBNamesToTSVFile(filePath):
    bedIntersectDf = pd.read_csv(filePath, sep='\t')
    bedIntersectDf = bedIntersectDf.iloc[np.random.choice(bedIntersectDf.index,
                                              int(len(bedIntersectDf)*1))]
    print("Choosing new gene names")
    bedIntersectDf['geneFullName'] = bedIntersectDf.apply(lambda row: decideBetwenGeneNames(row), axis=1)
    print("Droping gene names that were not found in tissues")
    newBedIntersectDf = bedIntersectDf[bedIntersectDf.geneFullName != '']
    print(str(len(bedIntersectDf.index)-len(newBedIntersectDf.index)) + " rows dropped because of unknown gene name")
    bedIntersectDf = newBedIntersectDf
    print("Droping useless columns")
    bedIntersectDf.drop('columnX', axis=1, inplace=True)
    bedIntersectDf.drop('columnY', axis=1, inplace=True)
    bedIntersectDf.drop('columnZ', axis=1, inplace=True)
    bedIntersectDf.drop('tfbsChr', axis=1, inplace=True)
    bedIntersectDf.drop('tfbsPosB', axis=1, inplace=True)
    bedIntersectDf.drop('tfbsPosA', axis=1, inplace=True)
    bedIntersectDf.drop('geneChr', axis=1, inplace=True)
    bedIntersectDf.drop('genePosA', axis=1, inplace=True)
    bedIntersectDf.drop('genePosB', axis=1, inplace=True)
    return bedIntersectDf

def startBsFrequency(row):
    if (row['tfName'],row['geneFullName']) not in bindingSiteFrequency:
        bindingSiteFrequency[(row['tfName'],row['geneFullName'])] = 0

def increaseBsFrequency(row):
    bindingSiteFrequency[(row['tfName'],row['geneFullName'])] += 1    

readGeneNames()
print("Treating data:")
treatedDf = addGeneAandBNamesToTSVFile(bedIntersectPath)
print("Starting to count binding site frequencies")
treatedDf.apply(lambda row: startBsFrequency(row), axis=1)
print("Counting binding site frequencies")
treatedDf.apply(lambda row: increaseBsFrequency(row), axis=1)
rows = []
for entry in bindingSiteFrequency:
    newRow = dict()
    newRow['tfName'] = entry[0]
    newRow['geneName'] = entry[1]
    newRow['count'] = bindingSiteFrequency[entry]
    rows.append(newRow)

treatedDf = pd.DataFrame(rows, columns=['tfName','geneName','count'])
print("Treated data")
print("Writing data")
treatedDf.to_csv(filteredBedIntersectPath, sep='\t', index=False)