from __future__ import division
import pandas as pd
import numpy as np
from TumorTruths import TumorTruths
from sklearn.metrics import f1_score
from DistanceClustering import DistanceClustering
from sklearn.metrics.pairwise import manhattan_distances

class SomaticMutationClassifier:
    
    def __init__(self, num_clusters, percentile):
        self.num_clusters = num_clusters
        self.transformToPredictions = np.vectorize(lambda x: 1 if x >= percentile else 0) 
    
    def fit(self, df):
        kMediansClustering = DistanceClustering(self.num_clusters, manhattan_distances, np.median)
        self.clusters, self.centers = kMediansClustering.fit(df.T.values)
        
    def predict(self, df):
        predictions = np.zeros(len(df))
        callersInPipeline = list()
        
        for i in range(0, self.num_clusters):
            weight = len([1 for num in self.clusters  if num == i])
            distances = manhattan_distances(df.T, self.centers[i])
            bestIndex = np.argmin(distances)
            bestCallerName = df.columns[bestIndex]
            representativeCaller = df[bestCallerName]
            callersInPipeline.append(bestCallerName)
            predictions += weight * representativeCaller

        predictions /= len(self.clusters)
        predictions = self.transformToPredictions(predictions) 
        
        return predictions, callersInPipeline

svnFiles = ['SNVCalls_IS1.txt', 'SNVCalls_IS2.txt', 'SNVCalls_IS3.txt', 'SNVCalls_IS4.txt']
truthFiles = ['synthetic.challenge.set1.tumor.all.truth.vcf', 'synthetic.challenge.set2.tumor.all.truth.vcf', 'synthetic.challenge.set3.tumor.20pctmasked.truth.vcf', 'synthetic.challenge.set4.tumour.25pctmasked.truth.vcf']

for i in range(len(svnFiles)):
    svnCalls = pd.read_csv(svnFiles[i], delimiter='\t', dtype={'CHROM':pd.np.str})
    truths = TumorTruths(truthFiles[i])

    resultsArray = truths.GetTruths(svnCalls)
    svnCalls = svnCalls[svnCalls.columns[3:-1]]
    
    somaticMutationClassifier = SomaticMutationClassifier(num_clusters=5, percentile=0.40)
    somaticMutationClassifier.fit(svnCalls)
    predictions, callersInPipeline = somaticMutationClassifier.predict(svnCalls)
    print f1_score(predictions, resultsArray)
