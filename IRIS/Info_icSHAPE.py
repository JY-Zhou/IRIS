import numpy as np
from scipy.stats import beta
import matplotlib.pyplot as plt

from IRIS.Utils import *

import RNA


class Info_icSHAPE:
    
    def __init__(self, icSHAPE_score_entry):
        substr = icSHAPE_score_entry.split('\t')
        
        self.rpkm = float(substr[2])
        self.scores = np.array([np.nan for _ in range(len(substr)-3)], dtype = np.float32)
        for i in range(len(self.scores)):
            score = substr[i+3]
            if score != 'NULL':
                self.scores[i] = float(score)
    
    
    
class Train_icSHAPE_Distribution:
    
    def __init__(self, rna_for_training, bprna_structure_dir):
        print('=== Load training data for icSHAPE score distribution ===')
        print('--- Collect %d RNA with curated secondary structure in %s' % (len(rna_for_training), bprna_structure_dir))
        
        self.paired_scores = []
        self.unpaired_scores = []
        for rna in rna_for_training:
            bprna_file_name = rna.ident.split('|')[0]
            # Collect the curated RNA secondary structure from bpRNA interface.
            with open(bprna_structure_dir + '/' + bprna_file_name + '.dbn', 'r') as f:
                # The annotated RNA secondary structure is located at the 5th line in DBN files.
                next(f); next(f); next(f); next(f);
                
                # Covert the dot-bracket notation into pair table representation
                ptable = RNA.ptable(f.readline().strip())
                for i in range(1, len(ptable)):
                    icSHAPE_score = rna.icSHAPE.scores[i-1]
                    
                    # Give an infinitesimal offset to the boundary value
                    if not np.isnan(icSHAPE_score):
                        if icSHAPE_score < Utils.EPS:
                            icSHAPE_score = Utils.EPS
                        if icSHAPE_score > 1.0 - Utils.EPS:
                            icSHAPE_score = 1.0 - Utils.EPS
                            
                        # Collect icSHAPE scores for paired/unpaired bases
                        if ptable[i] > 0:
                            self.paired_scores.append(icSHAPE_score)
                        else:
                            self.unpaired_scores.append(icSHAPE_score)
                        
        print('--- %d paired bases\n--- %d unpaired bases' % (len(self.paired_scores), len(self.unpaired_scores)))
        print('=== Finished ===\n\n')
        
        
    def learn_parameters(self):
        print('=== Train icSHAPE score distribution ===')
        self.parameters = {}
        # Fit scores for paired/unpaired bases
        self.parameters['paired_beta'] = beta.fit(self.paired_scores, floc = 0, fscale = 1)
        self.parameters['unpaired_beta'] = beta.fit(self.unpaired_scores, floc = 0, fscale = 1)
        for pars in self.parameters:
            print('--- %s\t' % pars, self.parameters[pars])
        print('=== Finished ===\n\n')
        
        
    def plot_distributions(self):
        fig, ax1 = plt.subplots(figsize = (7,4), constrained_layout = True)
        # Plot the histogram of icSHAPE scores
        ax1.hist(self.paired_scores, color = Utils.VRED, bins = 50, alpha = 0.5, density = True, label = 'Paired')
        ax1.hist(self.unpaired_scores, color = Utils.VBLUE, bins = 50, alpha = 0.5, density = True, label = 'Unpaired')
        ax1.set_xlim([-0.05, 1.05])
        ax1.set_xlabel('icSHAPE score', fontsize = 14)
        ax1.set_ylabel('# of bases (normalized)', fontsize = 14)
        ax1.legend(fontsize = 14)
        # Plot the probability density function of the fitted Beta distribution
        ax2 = ax1.twinx()
        xindex = np.linspace(0, 1, num = 100)
        ax2.plot(xindex, beta.pdf(xindex, *self.parameters['paired_beta']), color = Utils.VRED)
        ax2.plot(xindex, beta.pdf(xindex, *self.parameters['unpaired_beta']), color = Utils.VBLUE)
        ax2.set_ylabel('Probability density', fontsize = 14)
        plt.show()
    
    
    def compute_log_likelihood(self, ensemble):
        bp_vector = ensemble.get_bpp_per_base()
        log_likelihood = 0.0
        for i in range(ensemble.rna.N):
            icSHAPE_score = ensemble.rna.icSHAPE.scores[i]
                    
            # Give an infinitesimal offset to the boundary value
            if not np.isnan(icSHAPE_score):
                if icSHAPE_score < Utils.EPS:
                    icSHAPE_score = Utils.EPS
                if icSHAPE_score > 1.0 - Utils.EPS:
                    icSHAPE_score = 1.0 - Utils.EPS 
                
                # Compute log-likelihood for base i
                base_log_likelihood = 0.0
                # Precision control by cutting off an infinitesimal
                if bp_vector[i] > Utils.EPS:
                    # The probability of the icSHAPE score from paired base
                    prob_paired = beta.pdf(icSHAPE_score, *self.parameters['paired_beta'])
                    base_log_likelihood += bp_vector[i] * prob_paired
                # Precision control by cutting off an infinitesimal
                if 1.0 - bp_vector[i] > Utils.EPS:
                    # The probability of the icSHAPE score from unpaired base
                    prob_unpaired = beta.pdf(icSHAPE_score, *self.parameters['unpaired_beta'])
                    base_log_likelihood += (1.0 - bp_vector[i]) * prob_unpaired

                log_likelihood += np.log(base_log_likelihood)
            
        return log_likelihood