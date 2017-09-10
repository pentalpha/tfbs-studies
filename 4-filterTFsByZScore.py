# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

zScoreOfTFxTissuePath = "../results/zScoreOfTFxTissue.tsv"
zScoreOfTFxTissueTopPath = "../results/zScoreOfTFxTissue-Top.tsv"
zScoreOfTFxTissueBottomPath = "../results/zScoreOfTFxTissue-Bottom.tsv"
tissueNames = ["adipose_tissue", "adrenal_gland", "brain", "breast", "colon",
               "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary",
               "prostate", "skeletal_muscle", "testis", "thyriod"]

zScoreDf = pd.read_csv(zScoreOfTFxTissuePath, sep='\t')

#top and bottom percentage to be taken:
percent = 0.05

zScoreTopQuantile = zScoreDf['zScore'].quantile(1.0-percent)
zScoreBottomQuantile = zScoreDf['zScore'].quantile(percent)
zScoreTop = zScoreDf[zScoreDf.zScore >= zScoreTopQuantile]
zScoreTop.to_csv(zScoreOfTFxTissueTopPath, sep="\t", index=False)
zScoreBottom = zScoreDf[zScoreDf.zScore <= zScoreBottomQuantile]
zScoreBottom.to_csv(zScoreOfTFxTissueBottomPath, sep="\t", index=False)
