import pandas as pd
import vcf
import numpy as np

#svnCalls= pd.read_csv('SNVCalls_IS1.txt', delimiter='\t', dtype={'CHROM':pd.np.str})
#print svnCalls.columns
results = vcf.Reader(open('synthetic.challenge.set1.tumor.all.truth.vcf', 'rb'))
mutations = set()
for record in results:
    mutations.add(record.POS)
   # break
#genomicFeatures = pd.read_csv('GenomicFeatures_IS1.txt', delimiter='\t')
