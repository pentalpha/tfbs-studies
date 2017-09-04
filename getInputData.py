from pandas import pandas as pd
import numpy as np

bedIntersectPath = "../results/bedIntersectWaWbTFBSinGenes.bed"
filteredBedIntersectPath = "../results/bedIntersectWaWbTFBSinGenesFiltered.tsv"
tfbsCountsInTissuesPath = "../results/tfbsCountsInTissues.tsv"
tfAndTissuesPath = "../results/tfAndTissues.tsv"
tfAndTissuesTopPath = "../results/tfAndTissues-top.tsv"

genesPerTissueFolder = "../input/genesPerTissue"
tissueNames = ["adrenal_gland", "brain", "breast", "colon", "heart", "kidney", "leukocyte", "liver", "lung", "lymph_node", "ovary", "prostate", "skeletal_muscle", "testis", "thyriod"]
genesPerTissue = dict()
tissuePerGene = dict()
geneNames = set()
tfsInTissue = dict()
tfs = set()

def readGenesFromTissue(tissueName):
        if tissueName not in genesPerTissue:
            genesPerTissue[tissueName] = list()
            
        fileName = genesPerTissueFolder + "/" + tissueName + ".txt"
        f = open(fileName, 'r')
        
        for line in f:
            gene = line.replace('\n', '')
            genesPerTissue[tissueName].append(gene)
            
            if gene not in tissuePerGene:
                tissuePerGene[gene] = list()
            tissuePerGene[gene].append(tissueName)
            
            geneNames.add(gene)
        f.close()

def readAllTissueGenes():
        for tissueName in tissueNames:
                readGenesFromTissue(tissueName)
                
def decideBetwenGeneNames(row):
    names = row['geneFullName'].split(sep=',')
    for name in names:
        if name in geneNames:
            return name
    return ''

def addGeneAandBNamesToTSVFile(filePath):
    bedIntersectDf = pd.read_csv(filePath, sep='\t')
    
    print("Choosing new gene names")
    bedIntersectDf['geneFullName'] = bedIntersectDf.apply(lambda row: decideBetwenGeneNames(row), axis=1)
    print("Droping gene names that were not found in tissues")
    bedIntersectDf = bedIntersectDf[bedIntersectDf.geneFullName != '']
    print("Droping useless columns")
    bedIntersectDf.drop('columnX', axis=1, inplace=True)
    bedIntersectDf.drop('columnY', axis=1, inplace=True)
    bedIntersectDf.drop('columnZ', axis=1, inplace=True)
    bedIntersectDf.drop('tfbsChr', axis=1, inplace=True)
    bedIntersectDf.drop('tfbsPosB', axis=1, inplace=True)
    bedIntersectDf.drop('tfbsPosA', axis=1, inplace=True)
    bedIntersectDf.drop('geneChr', axis=1, inplace=True)
    bedIntersectDf.drop('genePosA', axis=1, inplace=True)
    bedIntersectDf.drop('genePosB', axis=1, inplace=True)
    return bedIntersectDf

def createFilteredDf():
    print("Treating data:")
    treatedDf = addGeneAandBNamesToTSVFile(bedIntersectPath)
    print("Treated data")
    print("Writing data")
    treatedDf.to_csv(filteredBedIntersectPath, sep='\t', index=False)

def addTfbsToTfCounts(row):
    geneName = row['geneFullName']
    tfName = row['tfName']
    tfs.add(tfName)
    for tissue in tissueNames:
        if geneName in genesPerTissue[tissue]:
            if tfName not in tfsInTissue[tissue]:
                tfsInTissue[tissue][tfName] = 1
            else:
                tfsInTissue[tissue][tfName] += 1

def mostUsedInTissue(row):
    mostUsed = tissueNames[0]
    for name in tissueNames:
        if row[name] > row[mostUsed]:
            mostUsed = name
    return mostUsed

def calcDiffFactor(row):
    mostUsed = row['mostUsedInTissue']
    mean = row['mean']
    return (row[mostUsed]/mean)

def countTfbsInTissues(df):
    for tissueName in tissueNames:
        tfsInTissue[tissueName] = dict()
    
    df.apply(lambda row: addTfbsToTfCounts(row), axis=1)
    
    print("Creating dataframe")
    cols = ['tfName']
    for tissueName in tissueNames:
        cols.append(tissueName)
    cols.append('totalBS')
    rows = []
    for tfName in tfs:
        row = dict()
        row['tfName'] = tfName
        totalBS = 0
        for tissueName in tissueNames:
            amount = 0
            if tfName in tfsInTissue[tissueName]:
                amount = tfsInTissue[tissueName][tfName]
            row[tissueName] = amount
            totalBS += amount
        row['totalBS'] = totalBS
        rows.append(row)
    tfDf = pd.DataFrame(rows, columns=cols)
    tfDf['median'] = tfDf[np.array(list(tissueNames))].median(axis=1)
    tfDf['mean'] = tfDf[np.array(list(tissueNames))].mean(axis=1)
    #tfDf['mostUsedInTissue'] = tfDf.apply(lambda row: mostUsedInTissue(row), axis=1)
    #tfDf['diffFactor'] = tfDf.apply(lambda row: calcDiffFactor(row), axis=1)
    #tfDf.sort(['diffFactor'], inplace=True, ascending=False)
    return tfDf

def getUniquenessPerTF(row, newDFRows):
    mean = row["mean"]
    for tissueName in tissueNames:
        tfbsCount = row[tissueName]
        newRow = dict()
        diffFactor = (tfbsCount / mean)
        newRow["tfName"] = row["tfName"]
        newRow["tissue"] = tissueName
        newRow["diffFactor"] = diffFactor
        newDFRows.append(newRow)
        

def calcExpressionForEachTissue(tFactorDF):
    newDFRows = []
    tFactorDF.apply(lambda row: getUniquenessPerTF(row, newDFRows), axis=1)
    diffFactorDF = pd.DataFrame(newDFRows, columns=['tfName', 'tissue', 'diffFactor'])
    diffFactorDF.sort(['diffFactor'], inplace=True, ascending=False)
    return diffFactorDF
        
print("Reading tissue genes lists")
readAllTissueGenes()
rawDF = pd.read_csv(filteredBedIntersectPath, sep='\t')
tFactorDF = countTfbsInTissues(rawDF)
print(tFactorDF.head())
expressionOnTissuesDF = calcExpressionForEachTissue(tFactorDF)

print("Writing result files")
tFactorDF.to_csv(tfbsCountsInTissuesPath, sep='\t', index=False)
#expressionOnTissuesDF.to_csv(tfAndTissuesPath, sep='\t', index=False)

diffFactorQuantile = expressionOnTissuesDF['diffFactor'].quantile(0.95)
tfAndTissuesTop = expressionOnTissuesDF[expressionOnTissuesDF.diffFactor >= diffFactorQuantile]
#tfAndTissuesTop.to_csv(tfAndTissuesTopPath, sep='\t', index=False)
#createFilteredDf()

#df_sliced.to_csv(filteredBedIntersectPath, sep='\t')
#print(str(genesIdentified/genesRead))
