# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import sys

def main():
    filteredBedIntersect = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
    tissueNames = ["adipose_tissue", "adrenal_gland", "brain", "breast", "colon",
               "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary",
               "prostate", "skeletal_muscle", "testis", "thyriod"]

    print("<Reading input data>")
    filteredBedIntersect = pd.read_csv(filteredBedIntersect, sep="\t")
    for tissue in tissueNames:
        tissueGenes = "../input/genesPerTissue/" + tissue + ".txt"
        tissueTFs = "../results/TFsPerTissue/" + tissue + "_tfs.txt"
        tfs = pd.DataFrame()
        tissueGenes = pd.read_csv(tissueGenes, header=None)
        tissueGenes = tissueGenes[0].values
        for gene in tissueGenes:
            tfsGene = filteredBedIntersect.loc[filteredBedIntersect['geneFullName'] == gene]
            tfsGene = tfsGene.drop('geneFullName', 1)
            tfsGene = tfsGene.reset_index(drop=True)
            tfs = tfs.append(tfsGene)

        tfs = pd.DataFrame(tfs).to_csv(tissueTFs, sep="\t", header=None, index=False)
        print("</Dataframe saved at " + tissueTFs + ">")

if __name__ == "__main__":
    main()
