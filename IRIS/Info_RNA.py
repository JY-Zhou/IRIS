from IRIS.Utils import *
from IRIS.Info_PARIS import *
from IRIS.Info_icSHAPE import *
from IRIS.Info_Coevolution import *

import RNA


class Info_RNA:
    
    def __init__(self, ident, seq):
        self.ident = str(ident)
        self.seq = str(seq)
        self.seq = self.seq.replace('T', 'U')
        self.N = len(self.seq)
        self.PARIS = None
        self.icSHAPE = None
        self.coevolution = None
        

class Structure_Ensemble:
    
    def __init__(self, rna):
        self.rna = rna
        
        
    def clear(self):
        self.K = 0
        self.structures = []
        self.alpha = np.zeros(0)
        
    
    def load(self, structure_list, alpha_list):
        self.clear()
        self.K = len(structure_list)
        self.structures = structure_list
        self.alpha = alpha_list
        
        
    def add(self, structure, alpha):
        self.K += 1
        self.structures.append(structure)
        self.alpha = np.append(self.alpha, alpha)
    
    
    def normalize_alpha(self):
        return self.alpha / np.sum(self.alpha)
    
    
    def assign(self, alpha):
        if alpha.size != self.K:
            raise ValueError('Error: the size of alpha not match!')
            
        self.alpha = alpha
        
        
    def get_tensor(self):
        X = np.zeros((self.rna.N, self.rna.N, self.K))
        for k in range(self.K):
            ptable = RNA.ptable(self.structures[k])
            for i in range(1, len(ptable)):
                if ptable[i] > 0:
                    X[i-1, ptable[i]-1, k] = 1
        return X
    
    
    def get_bpp_per_pair(self):
        alpha = self.normalize_alpha()
        bpp_matrix = np.zeros((self.rna.N, self.rna.N))
        for k in range(self.K):
            ptable = RNA.ptable(self.structures[k])
            for i in range(1, len(ptable)):
                if ptable[i] > 0:
                    bpp_matrix[i-1, ptable[i]-1] += alpha[k]
        return bpp_matrix
    
    
    def get_bpp_per_base(self):
        alpha = self.normalize_alpha()
        bpp_vector = np.zeros(self.rna.N)
        for k in range(self.K):
            ptable = RNA.ptable(self.structures[k])
            for i in range(1, len(ptable)):
                if ptable[i] > 0:
                    bpp_vector[i-1] += alpha[k]
        return bpp_vector
    
    
    def get_mean_free_enery(self):
        alpha = self.normalize_alpha()
        mean_free_energy = 0
        for k in range(self.K):
            mean_free_energy += alpha[k] * RNA.eval_structure_simple(self.rna.seq, self.structures[k])
        return mean_free_energy
    
    
    def get_mean_PARIS_support(self, PARIS_support):
        alpha = self.normalize_alpha()
        X = self.get_tensor()
        mean_PARIS_support = 0
        for k in range(self.K):
            mean_PARIS_support += alpha[k] * np.mean(PARIS_support[X[:,:,k] > 0])
        return mean_PARIS_support
    
    
    def show(self):
        # Sort representative structures by proportions
        alpha = self.normalize_alpha()
        sorted_ensemble = sorted([(alpha[k], self.structures[k]) for k in range(self.K)], key = lambda x: - x[0])
        
        for k in range(self.K):
            print('+++ Structure %d: %f\n%s' % (k+1, *sorted_ensemble[k]))
            
    
    def plot(self, PARIS_support, output_dir):
        # Sort representative structures by proportions
        alpha = self.normalize_alpha()
        sorted_ensemble = sorted([(alpha[k], self.structures[k]) for k in range(self.K)], key = lambda x: - x[0])
        
        # Plot structures as rainbow plots
        figure_paths = []
        plotter = Structure_Plotter(self.rna.seq, PARIS_support)
        for k in range(self.K):
            output_file = output_dir + '/Structure%d_%.3f.png' % (k+1, sorted_ensemble[k][0])
            print('--- Plotting structure %d (%.3f) ... ' % (k+1, sorted_ensemble[k][0]))
            plotter.plot(sorted_ensemble[k][1], output_file)
            figure_paths.append(output_file)
        return figure_paths
    
    
    def save(self, output_file):
        # Sort representative structures by proportions
        alpha = self.normalize_alpha()
        sorted_ensemble = sorted([(alpha[k], self.structures[k]) for k in range(self.K)], key = lambda x: - x[0])
        with open(output_file, 'w') as f:
            for k in range(self.K):
                print('+++ Structure %d:\t%f\n%s' % (k+1, *sorted_ensemble[k]), file = f)