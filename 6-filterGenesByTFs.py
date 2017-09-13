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
    filteredBedIntersect = filteredBedIntersect['tfName'].value_counts().reset_index()
    filteredBedIntersect.columns = ['tfName', 'all_genes']

    for tissue in tissueNames:
        tissueTFs = "../results/TFsPerTissue/" + tissue + "_tfs.txt"
        tissueTFs = pd.read_csv(tissueTFs, header=None)
        tissueTFs = tissueTFs[0].value_counts().reset_index()
        tissueTFs.columns = ['tfName', 'frequency']
        frequencyTFs = []
        for tf in filteredBedIntersect['tfName']:
            tf = tissueTFs.loc[tissueTFs['tfName'] == tf].reset_index(drop=True)
            if tf.empty:
                frequencyTFs.append(0)
            else:
                frequencyTFs.append(tf['frequency'][0])
            print (frequencyTFs)
        filteredBedIntersect[tissue] = frequencyTFs

    filteredBedIntersect = pd.DataFrame(filteredBedIntersect).to_csv(tfByGenesAndTissues, sep="\t", index=False)
    print("</Dataframe saved at " + tfByGenesAndTissues + ">")

if __name__ == "__main__":
    main()
