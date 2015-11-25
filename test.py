import pandas as pd
import vcf
import numpy as np
import sklearn.metrics

def GetResultsArray(svnCalls, results):
    
    mutations = set()
    for record in results:
        mutations.add(record.POS)
        
    truthColumn = np.zeros(len(svnCalls))

    i = 0
    for pos in svnCalls['END']:
        if pos in mutations:
            truthColumn[i] = 1
        i += 1
        
    return truthColumn

svnCalls = pd.read_csv('SNVCalls_IS1.txt', delimiter='\t', dtype={'CHROM':pd.np.str})
results = vcf.Reader(open('synthetic.challenge.set1.tumor.all.truth.vcf', 'rb'))
resultsArray = GetResultsArray(svnCalls, results)

scoresByCaller = dict()
for caller in svnCalls.columns[3:-1]:
    score = sklearn.metrics.f1_score(svnCalls[caller], resultsArray)
    scoresByCaller[caller] = score

sorted_scores = sorted(scoresByCaller, key=scoresByCaller.get, reverse=True)
bestScore = scoresByCaller[sorted_scores[0]]
worstScore = scoresByCaller[sorted_scores[5]]
print sorted_scores[0:5]
weights = list()
for caller in sorted_scores[0:5]:
    score = scoresByCaller[caller]
    weight = (1/ (bestScore - worstScore)) * score  - worstScore
    print weight
