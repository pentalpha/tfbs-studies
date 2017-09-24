# tfbs-studies
Bioinformatics study on transcription factors binding sites (TFBS) and their association with specific tissues.

This repository needs another directory, at the same parent directory, called "input" (or "../input", relatively"). It is available at http://177.20.147.141/~pitagoras/TF-findings/input/.

The results are stored in the "../results/" directory. The scripts must be executed in order(1-..., 2-..., 3-... etc).

## Tasks:
- [Task 1. Finding TFBSs for each gene](https://github.com/pentalpha/tfbs-studies#task-1-finding-tfbss-for-each-gene);
- [Task 2: Finding Transcription Factors related with tissues using FPKM values](https://github.com/pentalpha/tfbs-studies#task-2-finding-transcription-factors-related-with-tissues-using-fpkm-values);
- [Task 3A: Find TFs related to tissues using χ2 test](https://github.com/pentalpha/tfbs-studies#task-3a-find-tfs-related-to-tissues-using-χ2-test);
- [Task 3B: Finding TFs related to tissues using Monte Carlo approach](https://github.com/pentalpha/tfbs-studies#task-3b-finding-tfs-related-to-tissues-using-monte-carlo-approach-compare-with-random-samples);

## Task 1. Finding TFBSs for each gene

First we need to know that **transcription** is the process where a gene's DNA sequence is copied (transcribed) into 
an RNA molecule. So when enzyme RNA polymerase, which makes a new RNA molecule from a DNA template, his must attach 
to the DNA of the gene. It attaches at a spot called the **promoter** and his can attach to the promoter region only
with the help of proteins called **transcription factors** inside their **binding sites**. So transcription factors
are proteins that regulate the transcription of genes (activating or repressing), that is their copying into RNA, on
the way to making a protein. So the binding sites are regions that are linked to transcriptional regulation of a gene. 

For this task, our input data is:
- Transcription factors file: TFBS-ENCODE.bed
- Genes reference file: coord_upstream5kb_refGene.bed

This command has shown to be efficient at intersecting the two inputs ([0-createIntersectionBed.sh](0-createIntersectionBed.sh)):
```sh
bedtools intersect -wa -wb -a TFBS-ENCODE.bed -b coord_upstream5kb_refGene.bed > bedIntersectWaWbTFBSinGenes.bed
```

Ouput size: +- 196MB

1.969.951 lines, almost 2 million transcription factor binding sites related to some gene on the reference input.

Available at: http://177.20.147.141/~pitagoras/TF-findings/results/bedIntersectWaWbTFBSinGenes.bed

## Task 2: Finding Transcription Factors related with tissues using FPKM values
Input data:
- Human TFBSs intersected with genes in .bed ([link](http://177.20.147.141/~pitagoras/TF-findings/results/bedIntersectWaWbTFBSinGenes.bed))
- Expression (FPKM) of human genes in each tissue, tabular table in .txt ([link](http://177.20.147.141/~pitagoras/TF-findings/input/table_Human_body_map_ze_FPKM.txt))

Desired Output:
- TFs most or less represented in a single tissue, compared with the others

The idea here is finding out Transcription Factors which, compared with their representation in all tissues, are very or less related to specific tissues. Since this is a more complex task, we will separate it in several subtasks:

### Task 2.1: Filtering intersection .bed archive from previous task
The intersection .bed archive has a lot of information that are not useful for this task. We just want to know if the transcription factors had intersections with genes (and what genes). So we will delete most columns, leaving just:

tfName | geneName | count
--- | --- | ---
TAF1 | MIR1302-11 | 12
TAF1 | WASH7P | 5
HEY1 | MIR1302-2 | 1
HEY1 | MIR1302-10 | 2
... | ...

All this says is that MIR1302-11 and WASH7P have binding sites to TAF1, MIR1302-2 and MIR1302-10 have to HEY1 and so on...
The names at "gene" are an abreviation of the gene name column at the original .bed file. Gene names not found at the reference FPKM file are discarted.

For each individual gene, there can be more than one binding site to the same transcription factor. The 'count' column is this number of binding sites.

The script used to do this was [1-filterIntersectionBed.py](1-filterIntersectionBed.py).

Output:
- [bedIntersectWaWbTFBSinGenesFiltered.tsv](http://177.20.147.141/~pitagoras/TF-findings/results/bedIntersectWaWbTFBSinGenesFiltered.tsv)

### Task 2.2: Transcription Factor representativity in tissues
> Here’s how you do it for RPKM: First, count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor. Then, divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM). Finally, divide the RPM values by the length of the gene, in kilobases. This gives you RPKM. (RNA-Seq blog)

> FPKM is very similar to RPKM. RPKM was made for single-end RNA-seq, where every read corresponded to a single fragment that was sequenced. FPKM was made for paired-end RNA-seq. With paired-end RNA-seq, two reads can correspond to a single fragment, or, if one read in the pair did not map, one read can correspond to a single fragment. The only difference between RPKM and FPKM is that FPKM takes into account that two reads can map to one fragment (and so it doesn’t count this fragment twice). (also RNA-Seq blog)

Now we want to calculate, for each transcription factor, the sum of the FPKMs (of each tissue) of the genes which have binding sites to that transcription factor. The resulting table would have the following columns:

tfName | adipose_tissue | adrenal_gland | ... | genesWithBS | bindingSites | mean | stDeviation
--- | --- | --- | --- | --- | --- | ---
Transcription Factor name | Sum of the FPKMs in adipose tissue | the same for adrenal_gland | more sums for more tissues | Count of genes with binding sites to TF | Total of binding sites to the TF in all tissues | Mean FPKM sum | Standard deviation of FPKM sum

This is similar to the reference expression file with the FPKM values, but this table is about the transcription factors, not the genes.

The function to create the columns is:
```py
def getTFRow(tFactorName):
    print("\tStarting for " + tFactorName)
    tfDF = filteredBedIntersectDf[filteredBedIntersectDf.tfName == tFactorName]
    newRow = dict()
    newRow["tfName"] = tFactorName
    nGenes = 0
    bindingSites = 0
    for tissue in tissueNames:
        newRow[tissue] = 0.0
    for index, row in tfDF.iterrows():
        gName = row['geneName']
        nBindingSites = row['count']
        #this for loo is suposed to have only 1 iteration, unless there are 2
        #or more genes with he same name
        for geneIndex, geneRow in geneFpkmDf[geneFpkmDf.geneName == gName].iterrows():
            nGenes += 1
            bindingSites += nBindingSites
            for tissue in tissueNames:
                newRow[tissue] += geneRow[tissue]*nBindingSites
    newRow["genesWithBS"] = nGenes
    newRow["bindingSites"] = bindingSites
    print("\tDone for " + tFactorName)
    return newRow
```
In [2-makeTFactorDf.py](2-makeTFactorDf.py)

These columns are computed using Python's "multiprocessing.Pool", to maximize performance in multi-core systems.

Output:
- [tfFPKMinTissues.tsv](http://177.20.147.141/~pitagoras/TF-findings/results/tfFPKMinTissues.tsv)

### Task 2.3: Filtering pairs of Tissue with TF using Z-Score and Quantile
> When a frequency distribution is normally distributed, we can find out the probability of a score occurring by standardising the scores, known as standard scores (or z scores). The standard normal distribution simply converts the group of data in our frequency distribution such that the mean is 0 and the standard deviation is 1. (Laerd Statistics)

For a transcription factor *tf* and a tissue *ts*, the z-score is:
```py
z-score = (fpkmSum(tf, ts) - meanFpkmSum(tf)) / standardDeviationForFpkmSums(tf)
```
(Note that, for each tissue, the transcription factor will have a different Z-Score.)

We want to create a table with the following data:

tfName | tissue | genesWithBS | fpkmSum | mean | stDeviation | zScore
--- | --- | --- | --- | --- | --- | ---
Transcription factor name | Tissue name | Number of genes with binding sites to TF in this tissue | Sum of the tissue FPKM values of all the genes with a binding site to the TF | Mean FPKM sum of TF for all tissues | Standard deviation in FPKM sums | The Z-Score

The rows for such a table are created using the following function:
```py
def getZScoreRows(row, newDFRows):
    mean = row["mean"]
    std = row["stDeviation"]
    for tissueName in tissueNames:
        newRow = dict()
        newRow["tfName"] = row["tfName"]
        newRow["tissue"] = tissueName
        newRow["genesWithBS"] = row["genesWithBS"]
        newRow["fpkmSum"] = row[tissueName]
        newRow["mean"] = mean
        newRow["stDeviation"] = std
        if(std != 0.0):
            newRow["zScore"] = (row[tissueName] - mean) / std
        else:
            newRow["zScore"] = 0.0
        newDFRows.append(newRow)
```
In [3-zScoreOf-TFxTissue.py](3-zScoreOf-TFxTissue.py).

In these rows, the Z-Score quantifies how much the representativity of a TF in a tissue deviates from what would be normal for that TF. So the TFs we are looking for must have very big or very low Z-Scores. Now, all we need to do is filter the rows according to the Z-Scores.

For a given quantitative column in a Pandas DataFrame, the function "df['columnName'].quantile(0.XX)" returns a value in the column, so that XX% of the values are smaller or equal to it. Using this value, we can filter the top 5% and bottom 5% of the Z-Score rows:
```py
#top and bottom percentage to be taken:
percent = 0.05

zScoreTopQuantile = zScoreDf['zScore'].quantile(1.0-percent)
zScoreBottomQuantile = zScoreDf['zScore'].quantile(percent)
zScoreTop = zScoreDf[zScoreDf.zScore >= zScoreTopQuantile]
zScoreTop.to_csv(zScoreOfTFxTissueTopPath, sep="\t", index=False)
zScoreBottom = zScoreDf[zScoreDf.zScore <= zScoreBottomQuantile]
zScoreBottom.to_csv(zScoreOfTFxTissueBottomPath, sep="\t", index=False)
```
In [4-filterTFsByZScore.py](4-filterTFsByZScore.py).

Outputs:
- [zScoreOfTFxTissue-Top.tsv](http://177.20.147.141/~pitagoras/TF-findings/results/zScoreOfTFxTissue-Top.tsv)
- [zScoreOfTFxTissue-Bottom.tsv](http://177.20.147.141/~pitagoras/TF-findings/results/zScoreOfTFxTissue-Bottom.tsv)
- [zScoreOfTFxTissue.tsv](http://177.20.147.141/~pitagoras/TF-findings/results/zScoreOfTFxTissue.tsv)

## Task 3A: Find TFs related to tissues using χ2 test

### Task 3A.1 Make TFs files per Tissue
Using the file named bedIntersectWaWbTFBSinGenesFiltered.tsv, we can know the TFs related with every gene and and using the files .txt found [here](http://177.20.147.141/~pitagoras/TF-findings/input/genesPerTissue/), we can know the genes founded in every tissue. We just want to know what are the transcription factors pertaining in every tissue based on their genes .txt file. So for every gene pertained in a tissue, we will be going to add all the transcription factors related with that gene, leaving just:

| tfName   |
|----------|
| ...      |
| NF-YA    |
| Pol2-4H8 |
| CTCF     |
| CTCF     |
| ...      |

The script used to do this was [5-makeTFsPerTissue.py](5-makeTFsPerTissue.py).

Output:
- [adipose_tissue_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/adipose_tissue_tfs.txt)
- [adrenal_gland_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/adrenal_gland_tfs.txt)
- [brain_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/brain_tfs.txt)
- [breast_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/breast_tfs.txt)
- [colon_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/colon_tfs.txt)
- [heart_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/heart_tfs.txt)
- [kidney_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/kidney_tfs.txt)
- [leukocyte_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/leukocyte_tfs.txt)
- [liver_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/liver_tfs.txt)
- [lung_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/li_tfs.txt)
- [lymph_node_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/lymph_node_tfs.txt)
- [ovary_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/ovary_tfs.txt)
- [prostate_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/prostate_tfs.txt)
- [skeletal_muscle_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/skeletal_muscle_tfs.txt)
- [testis_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/testis_tfs.txt)
- [thyriod_tfs.txt](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/thyriod_tfs.txt)

### Task 3A.2 Filtering quantity of genes for every TF in every tissue
Using the file named bedIntersectWaWbTFBSinGenesFiltered.tsv, we can know the TFs related with every gene and using the files .txt found [here](http://177.20.147.141/~bif/luiseduardo/tfbs-studies/results/TFsPerTissue/), we can know all the TFs founded and their frequence in every tissue. We just want to know quantity of genes for every TF in every tissue. So for every TF found in a tissue, we will be going to count his frequency, leaving just:

| tfName   | all_genes | adipose_tissue | adrenal_gland | ... | testis | thyriod |
|----------|-----------|----------------|---------------|-----|--------|---------|
| Pol2     | 76789     | 12             | 26            | ... | 739    | 13      |
| Pol2-4H8 | 63732     | 4              | 30            | ... | 628    | 11      |
| TAF1     | 58550     | 3              | 26            | ... | 563    | 9       |
| ...      | ...       | ...            | ...           | ... | ...    | ...     |

The script used to do this was [6-filterGenesByTFs.py](6-filterGenesByTFs.py).

Output:
- [tfByGenesAndTissues.tsv](http://10.7.5.38/~bif/luiseduardo/tfbs-studies/results/tfByGenesAndTissues.tsv)

## Task 3B: Finding TFs related to tissues using Monte Carlo approach (compare with random samples)

Input:
- [Lists of genes related to each tissue](http://work.bioinformatics-brazil.org/~pitagoras/TF-findings/input/genesPerTissue/);
- Count of BS between to a TF in a gene: [bedIntersectWaWbTFBSinGenesFiltered.tsv](http://work.bioinformatics-brazil.org/~pitagoras/TF-findings/results/bedIntersectWaWbTFBSinGenesFiltered.tsv);

The input data is used to create several sets and dictionaries:
```py
#key: tissueName, value: set of genes representative to tissueName
genesPerTissue = dict()
#set of all geneNames
geneNames = set()
#set of tfNames
tfs = set()
#key: tfName, value: set of genes with at least one binding site to tfName
genesWithBSTo = dict()
#key: tissueName, value: list of random subsets of 'geneNames'
geneSamples = dict()
#key: (tfName, geneName), value: number of binding sites
tfbsCount = dict()
```

Our objective here is to, for each pair of tissue (TS) and transcription factor (TF), compare the number of binding sites to TF on the set of genes related TS, with the number of binding sites to TF on random sets of genes. To do that, we first create 2000 random sets of genes with the same size as the set of genes related to TS. Then we calculate the number of binding sites (to TF) on each of these random sets and save the results in an array. Our final number is the percentile of random binding site counts below or equal to the observed number of binding sites on the set of genes related to TS.

All that is calculated using the following method, given a pair of TSxTF:
```py
def getBsCountInSamples(tfName, tissueName):
    genesWithBS = genesWithBSTo[tfName]
    bsCountInSamples = []
    for sample in geneSamples[tissueName]:
        geneSet = sample.intersection(genesWithBS)
        bsCount = countBSToTFinGeneSet(geneSet, tfName)
        bsCountInSamples.append(bsCount)
    return bsCountInSamples
            
def qualifyTFAtTissue(tfName, tissueName):
    genesWithBS = genesWithBSTo[tfName]
    genesInTissueWithBs = genesPerTissue[tissueName].intersection(genesWithBS)
    bsCountInTissue = countBSToTFinGeneSet(genesInTissueWithBs, tfName)
    bsCountInSamples = getBsCountInSamples(tfName, tissueName)
    bsCountInSamples.append(bsCountInTissue)
    
    npArray = np.asarray(bsCountInSamples)
    std = npArray.std()
    median = np.median(npArray)
    perc = stats.percentileofscore(bsCountInSamples, bsCountInTissue, 'weak')
    return (perc, median, std, bsCountInTissue, len(genesInTissueWithBs))
```

Next, all we need to do is filter which pairs of TSxTF are relevant. The criterion is:
- Percentile being outside of the standard deviation (in the random samples);
- To be among the most representative to the tissue: percentile >= 95.0;
- To be among the least representative to the tissue: percentile <= 5.0;

Output:
- [mostRepresentativeTFsPerTissue.tsv](http://work.bioinformatics-brazil.org/~pitagoras/TF-findings/results/mostRepresentativeTFsPerTissue.tsv);
- [leastRepresentativeTFsPerTissue](http://work.bioinformatics-brazil.org/~pitagoras/TF-findings/results/leastRepresentativeTFsPerTissue.tsv);

All the code is in [7-selectTFsByRandomSample.py](7-selectTFsByRandomSample.py).

## References

Khan Academy. [Gene regulation](https://www.khanacademy.org/science/biology/gene-regulation).

RNA-Seq Blog. [RPKM, FPKM and TPM clearly explained](http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/).

Lard Statistics. [Standard Score](https://statistics.laerd.com/statistical-guides/standard-score.php).
