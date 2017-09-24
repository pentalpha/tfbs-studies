# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import scipy.stats as stats

def main():
    genesPerTissueFolder = "../input/genesPerTissue/"
    tfByGenesAndTissuesPath = "../results/tfByGenesAndTissues.tsv"
    chiSquaredTFsPerTissuePath = "../results/chiSquaredTFsPerTissue/"
    filteredBedIntersect = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
    tissueNames = ["adipose_tissue", "adrenal_gland", "brain", "breast", "colon",
                   "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary",
                   "prostate", "skeletal_muscle", "testis", "thyriod"]

    print("<Reading input data>")

    tfByGenesAndTissuesPath = pd.read_csv(tfByGenesAndTissuesPath, sep="\t")
    filteredBedIntersect = pd.read_csv(filteredBedIntersect, sep="\t")
    filteredBedIntersect = filteredBedIntersect['geneFullName'].value_counts().reset_index()
    genesAllTissue = filteredBedIntersect.shape[0]
    #Ignore when the code try to "divide by zero" or "divide by NaN".
    np.seterr(divide='ignore', invalid='ignore')

    for tissue in tissueNames:
          chiSquaredTFs = pd.DataFrame(columns=['tfName','Chi-Square Test','P-Value', 'State'])
          genesTissue = sum(1 for line in open("../input/genesPerTissue/" + tissue + ".txt"))
          state = ""
          for date, tf in tfByGenesAndTissuesPath.T.iteritems():
              observed = tf[tissue]
              expected = tf["allGenes"] * genesTissue/genesAllTissue
              chi_squared_test = stats.chisquare(f_obs = observed, f_exp = expected)[0]
              p_value = 1 - stats.chi2.cdf(x=chi_squared_test, df=1)
              if observed > expected:
                  state = "Depleted"
              else:
                  state = "Increased"
              if p_value <= 0.05:
                  chiSquaredTFs = chiSquaredTFs.append([{'tfName':tf["tfName"], 'Chi-Square Test':chi_squared_test, 'P-Value':p_value, 'State':state}])
          chiSquaredTFs = chiSquaredTFs.sort_values(by='P-Value')
          chiSquaredTFs = pd.DataFrame(chiSquaredTFs).to_csv(chiSquaredTFsPerTissuePath + tissue + ".tsv", sep="\t", index=False)
          print("</Dataframe saved at " + chiSquaredTFsPerTissuePath + tissue + ".tsv" + ">")

if __name__ == "__main__":
    main()
