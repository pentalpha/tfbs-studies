# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

tfFPKMInTissuesPath = "../results/tfFPKMinTissues.tsv"
zScoreOfTFxTissuePath = "../results/zScoreOfTFxTissue.tsv"
tissueNames = ["adipose_tissue", "adrenal_gland", "brain", "breast", "colon",
               "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary",
               "prostate", "skeletal_muscle", "testis", "thyriod"]

def getZScoreRows(row, newDFRows):
    mean = row["mean"]
    std = row["stDeviation"]
    for tissueName in tissueNames:
        newRow = dict()
        newRow["tfName"] = row["tfName"]
        newRow["tissue"] = tissueName
        newRow["bindingSites"] = row["bindingSites"]
        newRow["fpkmSum"] = row[tissueName]
        newRow["mean"] = mean
        newRow["stDeviation"] = std
        if(std != 0.0):
            newRow["zScore"] = (row[tissueName] - mean) / std
        else:
            newRow["zScore"] = 0.0
        newDFRows.append(newRow)

def calcExpressionForEachTissue(tFactorDF):
    newDFRows = []
    tFactorDF.apply(lambda row: getZScoreRows(row, newDFRows), axis=1)
    zScoreDf = pd.DataFrame(newDFRows, columns=['tfName', 'tissue', 'bindingSites',
                                                    'fpkmSum', 'mean', 'stDeviation', 
                                                    'zScore'])
    #zScoreDf.sort_values()
    zScoreDf.sort_values(['zScore'], inplace=True, ascending=False)
    return zScoreDf

fpkmDf = pd.read_csv(tfFPKMInTissuesPath, sep='\t')
zScoreDf = calcExpressionForEachTissue(fpkmDf)
zScoreDf.to_csv(zScoreOfTFxTissuePath, sep='\t', index=False)