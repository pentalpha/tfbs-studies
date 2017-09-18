# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

def main():
    filteredBedIntersect = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
    tfByGenesAndTissues = "../results/tfByGenesAndTissues.tsv"
    tissueNames = ["adipose_tissue", "adrenal_gland", "brain", "breast", "colon",
                   "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary",
                   "prostate", "skeletal_muscle", "testis", "thyriod"]

    print("<Reading input data>")
    filteredBedIntersect = pd.read_csv(filteredBedIntersect, sep="\t")
    
    #start dict with 0 values for the frequencies
    tfbsByTF = dict()
    for tfName in filteredBedIntersect['tfName'].unique():
        tfbsByTF[tfName] = dict()
        for tissueName in tissueNames:
            tfbsByTF[tfName][tissueName] = 0
            
    #read tf frequency on tissues from lists        
    for tissue in tissueNames:
        tissueTFsPath = "../results/TFsPerTissue/" + tissue + "_tfs.txt"
        tissueTFs = pd.read_csv(tissueTFsPath, sep='\t')
        for index, row in tissueTFs.iterrows():
            tfbsByTF[row['tfName']][tissue] += row['frequency']
    
    #create rows for new df        
    tfRows = []
    for tfName in filteredBedIntersect['tfName'].unique():
        newRow = dict()
        newRow['tfName'] = tfName
        for tissueName in tissueNames:
            newRow[tissueName] = tfbsByTF[tfName][tissueName]
        tfRows.append(newRow)
        
    tfDf = pd.DataFrame(tfRows, columns=(['tfName'] + tissueNames))
    tfDf.to_csv(tfByGenesAndTissues, sep="\t", index=False)
    print("</Dataframe saved at " + tfByGenesAndTissues + ">")

if __name__ == "__main__":
    main()
