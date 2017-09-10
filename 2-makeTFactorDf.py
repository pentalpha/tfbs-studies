# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from multiprocessing import Pool


filteredBedIntersectPath = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
humanGenesFPKMInTissues = "../input/table_Human_body_map_ze_FPKM.txt"
tfFPKMInTissuesPath = "../results/tfFPKMinTissues.tsv"
tissueNames = ["adipose_tissue", "adrenal_gland", "brain", "breast", "colon",
               "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary",
               "prostate", "skeletal_muscle", "testis", "thyriod"]

#random sample size: 15%
percent = 0.15

print("<Reading input data>")
filteredBedIntersectDf = pd.read_csv(filteredBedIntersectPath, sep="\t")
filteredBedIntersectDf = filteredBedIntersectDf.iloc[
        np.random.choice(filteredBedIntersectDf.index, int(len(filteredBedIntersectDf)*percent))]
geneTFRelations = set()
tfsOfGene = dict()

geneFpkmDf = pd.read_csv(humanGenesFPKMInTissues, sep="\t")
geneFpkmDf = geneFpkmDf.iloc[np.random.choice(geneFpkmDf.index,
                                              int(len(geneFpkmDf)*percent))]
print("</Loaded input data>")

def getTFRow(tFactorName):
    print("\tStarting for " + tFactorName)
    tfDF = filteredBedIntersectDf[filteredBedIntersectDf.tfName == tFactorName]
    newRow = dict()
    newRow["tfName"] = tFactorName
    nGenes = 0
    for tissue in tissueNames:
        newRow[tissue] = 0.0
    for index, row in tfDF.iterrows():
        gName = row['geneFullName']
        #this for loo is suposed to have only 1 iteration, unless there are 2
        #or more genes with he same name
        for geneIndex, geneRow in geneFpkmDf[geneFpkmDf.geneName == gName].iterrows():
            nGenes += 1
            for tissue in tissueNames:
                newRow[tissue] += geneRow[tissue]
    newRow["genesWithBS"] = nGenes
    print("\tDone for " + tFactorName)
    return newRow

print("<Creating rows>")
pool = Pool(processes = 5)
rows = pool.map(getTFRow, filteredBedIntersectDf.tfName.unique())
pool.close()
pool.join()
print("</Created rows>")

print("<Writing dataframe>")
columnNames = ["tfName"]
for tissue in tissueNames:
    columnNames.append(tissue)
columnNames.append("genesWithBS")
df = pd.DataFrame(rows, columns=columnNames)
print("Computing mean values column")
df['mean'] = df[np.array(list(tissueNames))].mean(axis=1)
df['stDeviation'] = df[np.array(list(tissueNames))].std(axis=1)
df.to_csv(tfFPKMInTissuesPath, sep="\t", index=False)
print("</Dataframe saved at " + tfFPKMInTissuesPath + ">")