from IRIS.Info_RNA import *

import sys

import RNA


class MFE:
    
    def __init__(self, rna, ostream = sys.stdout):
        self.rna = rna
        self.ostream = ostream
        self.ensemble = Structure_Ensemble(self.rna)
        
        
    def predict(self):
        self.ensemble.load([RNA.fold(self.rna.seq)[0]], np.array([1]))
        
    
class Nonredundant_Sampling:
    
    def __init__(self, rna, ostream = sys.stdout):
        self.rna = rna
        self.ostream = ostream
        self.ensemble = Structure_Ensemble(self.rna)
    
    
    def generate(self, K):
        md = RNA.md()
        md.uniq_ML = 1
        fc = RNA.fold_compound(self.rna.seq, md)
        (ss, mfe) = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        self.ensemble.load(list(fc.pbacktrack(K, RNA.PBACKTRACK_NON_REDUNDANT)), np.ones(K))
        