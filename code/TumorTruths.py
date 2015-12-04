import vcf
import numpy as np

class TumorTruths:

    def __init__(self, vcfFile):
        self.results = vcf.Reader(open(vcfFile, 'rb'))

    def GetTruths(self, svnCalls):
        mutations = set()
        for record in self.results:
	    mutations.add(record.POS)
		       
        truthColumn = np.zeros(len(svnCalls))

        i = 0
	for pos in svnCalls['END']:
	    if pos in mutations:
	        truthColumn[i] = 1
	    i += 1
								        
        return truthColumn
