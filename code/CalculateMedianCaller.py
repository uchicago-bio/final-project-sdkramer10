from __future__ import division
import pandas as pd
import numpy as np
from sklearn.metrics import f1_score
from TumorTruths import TumorTruths

svnFiles = ['SNVCalls_IS1.txt', 'SNVCalls_IS2.txt', 'SNVCalls_IS3.txt', 'SNVCalls_IS4.txt']
truthFiles = ['synthetic.challenge.set1.tumor.all.truth.vcf', 'synthetic.challenge.set2.tumor.all.truth.vcf', 'synthetic.challenge.set3.tumor.20pctmasked.truth.vcf', 'synthetic.challenge.set4.tumour.25pctmasked.truth.vcf']

for k in range(len(svnFiles)):

    svnCalls = pd.read_csv(svnFiles[k], delimiter='\t', dtype={'CHROM':pd.np.str})
    truths = TumorTruths(truthFiles[k])
    resultsArray = truths.GetTruths(svnCalls)

    scores = list()
    for column in svnCalls[svnCalls.columns[3:-1]]:
        scores.append(f1_score(svnCalls[column], resultsArray))
    
    print np.median(np.array(scores))
