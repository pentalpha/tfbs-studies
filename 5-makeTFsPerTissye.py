# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import sys

def main():
    filteredBedIntersect = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
    tissueGenes = "../input/genesPerTissue/" + sys.argv[1]
    tissueTFs = "../results/TFsPerTissue/" + sys.argv[1].split('.', 1)[0] + "_tfs.txt"

    print("<Reading input data>")
    filteredBedIntersect = pd.read_csv(filteredBedIntersect, sep="\t")
    TFsDF = pd.DataFrame()
    tissueGenes = pd.read_csv(tissueGenes, header=None)
    tissueGenes = tissueGenes[0].values
    for gene in tissueGenes:
        tfsGene = filteredBedIntersect.loc[filteredBedIntersect['geneFullName'] == gene]
        tfsGene = tfsGene.drop('geneFullName', 1)
        tfsGene = tfsGene.reset_index(drop=True)
        TFsDF = TFsDF.append(tfsGene)

    TFsDF = pd.DataFrame(TFsDF).to_csv(tissueTFs, sep="\t", header=None, index=False)
    print("</Dataframe saved at " + tissueTFs + ">")

if __name__ == "__main__":
    main()
