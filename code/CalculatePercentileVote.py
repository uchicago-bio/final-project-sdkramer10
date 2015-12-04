from __future__ import division
import pandas as pd
import numpy as np
import sys
from TumorTruths import TumorTruths
from sklearn.metrics import f1_score

def calculateF1Score(svnFileName, truthFileName, percentile):
    svnCalls = pd.read_csv(svnFileName, delimiter='\t', dtype={'CHROM':pd.np.str})
    truths = TumorTruths(truthFileName)
    resultsArray = truths.GetTruths(svnCalls)
    svnCalls = svnCalls[svnCalls.columns[3:-1]]

    numPositivePredictions = svnCalls.sum(axis=1)

    fractionalPercentile = percentile / 100
    transformToPredictions = np.vectorize(lambda x: 1 if x >= (fractionalPercentile * len(svnCalls.columns)) else 0) 
    result = transformToPredictions(numPositivePredictions)
    return f1_score(result, resultsArray)    

svnFiles = ['SNVCalls_IS1.txt', 'SNVCalls_IS2.txt', 'SNVCalls_IS3.txt', 'SNVCalls_IS4.txt']
truthFiles = ['synthetic.challenge.set1.tumor.all.truth.vcf', 'synthetic.challenge.set2.tumor.all.truth.vcf', 'synthetic.challenge.set3.tumor.20pctmasked.truth.vcf', 'synthetic.challenge.set4.tumour.25pctmasked.truth.vcf']

for k in range(len(svnFiles)):
    print calculateF1Score(svnFiles[k], truthFiles[k], float(sys.argv[1]))
