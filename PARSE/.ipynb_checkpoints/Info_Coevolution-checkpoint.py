import numpy as np

from PARSE.Utils import *


class Info_Coevolution:
    
    def __init__(self, rna, stockholm_dir, covariation_dir):
        self.rna = rna
        rfam_str = self.rna.ident.split('|')[1].split('_')
        self.rfam_id = rfam_str[0]
        self.rfam_acc = '/'.join(rfam_str[1:])
        
        # Load multiple sequence alignment
        msa_seq = ''
        with open(stockholm_dir + '/' + self.rfam_id + '.sto', 'r') as f:
            for line in f:
                line = line.strip()
                if self.rfam_acc in line:
                    msa_seq = line.strip().split(' ')[-1]
        
        # Load normalized mutual information matrix
        msa_cov_matrix = []
        with open(covariation_dir + '/' + self.rfam_id + '.cov', 'r') as f:
            for line in f:
                msa_cov_matrix.append([float(x) for x in line.strip().split('\t')])
        msa_cov_matrix = np.array(msa_cov_matrix)
        
        # Remove insertions in the multiple sequence alignment
        matched_matrix = []
        for i in range(len(msa_seq)):
            if msa_seq[i] != '-': # '-' means an insertion
                matched_matrix.append([])
                for j in range(len(msa_seq)):
                    if msa_seq[j] != '-':
                        if msa_cov_matrix[i, j] == 0:
                            matched_matrix[-1].append(Utils.EPS)
                        else:
                            matched_matrix[-1].append(msa_cov_matrix[i, j])
        self.covariation_matrix = np.array(matched_matrix)
    
    
    def evaluate_KL_divergence(self, ensemble):
        bpp_matrix = ensemble.get_bpp_per_pair()
        KL = 0.0
        for i in range(ensemble.rna.N):
            for j in range(ensemble.rna.N):
                if j > i:
                    a = bpp_matrix[i, j]
                    b = self.covariation_matrix[i, j]
                    if a > 0:
                        KL += a * np.log(a) - a * np.log(b) - a + b
        return KL
    
    
    def plot_covariation_matrix(self, xticks):
        fig, axs = plt.subplots(1, 1, figsize = (8, 4), constrained_layout=True)
        ax = axs
        ax.set_title('Evolutionary conservation', fontsize = 16)
        Utils.plot_triangle_matrix(ax, self.covariation_matrix, Utils.white2red_cmap)
        cbar1 = fig.colorbar(cm.ScalarMappable(mpl.colors.Normalize(0, np.max(self.covariation_matrix)), Utils.white2red_cmap), aspect = 20, ax = axs, location = 'left')
        cbar1.set_label(label = 'Normalized mutual information', fontsize = 16)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, y = 0, fontsize = 14)
        plt.show()