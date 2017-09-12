# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

def main():
    filteredBedIntersect = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
    tfByGenes = "../results/tfByGenes.tsv"

    print("<Reading input data>")
    filteredBedIntersect = pd.read_csv(filteredBedIntersect, sep="\t")
    print(filteredBedIntersect.columns)
    filteredBedIntersect = filteredBedIntersect['tfName'].value_counts().reset_index()
    filteredBedIntersect.columns = ['tfName', 'genes']

    df = pd.DataFrame(filteredBedIntersect).to_csv(tfByGenes, sep="\t", index=False)
    print("</Dataframe saved at " + tfByGenes + ">")

if __name__ == "__main__":
    main()
