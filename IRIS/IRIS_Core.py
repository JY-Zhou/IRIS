import sys
import numpy as np
from scipy.stats import expon, norm
import itertools
from sklearn.cluster import AgglomerativeClustering
from sklearn.linear_model import Lasso
from scipy.sparse import csc_matrix
from scipy.optimize import lsq_linear
import matplotlib as mpl


import RNA

from IRIS.Utils import *
from IRIS.Info_RNA import *

    
class Scoring:
    
    def __init__(self, rna, ostream = sys.stdout):
        self.rna = rna
        self.ostream = ostream
        self.PARIS_support = np.zeros((self.rna.N, self.rna.N))
        
    
    def compute_PARIS_support(self):
        print('=== Convert PARIS reads to PARIS support matrix ===', file = self.ostream)
        arm_len = []
        for (ll, lr), (rl, rr) in self.rna.PARIS.blocks:
            arm_len.append(lr-ll)
            arm_len.append(rr-rl)
            
        # The 3-sigma rule of Normal distribution
        sd = np.mean(arm_len)
        sd /= 6
        print('--- The standard deviation is %.3f' % sd, file = self.ostream)
            
        for (ll, lr), (rl, rr) in self.rna.PARIS.blocks:
            # The center of intervals as the mean value
            l_mu = (ll + lr) / 2
            r_mu = (rl + rr) / 2
            # Regions for counting read coverage
            ll = int(max(0, l_mu - 3*sd))
            lr = int(min(self.rna.N, l_mu + 3*sd))
            rl = int(max(0, r_mu - 3*sd))
            rr = int(min(self.rna.N, r_mu + 3*sd))
            l_vec = np.array([_ for _ in range(0, self.rna.N)])
            r_vec = np.array([_ for _ in range(0, self.rna.N)]) 
            # Transform via Normal distribution 
            l_support = norm.pdf(l_vec, l_mu, sd)
            r_support = norm.pdf(r_vec, r_mu, sd)
            # Normaling by the peak of the distribution
            l_support /= np.max(l_support)
            r_support /= np.max(r_support)
            # Sum up to the matrix
            self.PARIS_support += np.outer(l_support, r_support)
            self.PARIS_support += np.outer(r_support, l_support)
            
        # Trim infinitesimal values as 0
        self.PARIS_support[self.PARIS_support < Utils.EPS] = 0    
        # The logarithm transformation
        self.PARIS_support = np.log(self.PARIS_support+1)
        
        print('=== Finished! ===\n\n', file = self.ostream)
    
    
    def plot_PARIS_support(self, xticks):      
        # Draw the heatmap of triangle matrix
        fig, axs = plt.subplots(1, 1, figsize = (8, 4), constrained_layout=True)
        ax = axs
        Utils.plot_triangle_matrix(ax, self.PARIS_support, Utils.black2red_cmap)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, fontsize = 14)

        # Draw the color bar legend
        PARIS_support_cmap_norm = mpl.colors.Normalize(0, np.max(self.PARIS_support))
        cbar1 = plt.colorbar(cm.ScalarMappable(PARIS_support_cmap_norm, Utils.black2red_cmap), aspect = 20, ax = axs, location = 'left')
        cbar1.set_label(label = 'PARIS support', fontsize = 16)

        plt.show()
        
        
        
class Generating:
    
    def __init__(self, rna, PARIS_support, ostream = sys.stdout):
        self.rna = rna
        self.PARIS_support = PARIS_support
        self.ostream = ostream
        self.candidate_structures = []
    
    
    def __check_stem(self, left_part, right_part, k):
        if len(left_part) != k or len(right_part) != k:
            return False
        p, q = 0, k-1
        while p < k:
            if not left_part[p] + right_part[q] in Utils.VALID_BP:
                return False
            p, q = p + 1, q - 1
        return True
    
    
    def __collect_stems(self, fraction, k_range):
        positive_scores = self.PARIS_support[np.nonzero(self.PARIS_support)]
        threshold = np.quantile(positive_scores, fraction)
        print('--- The %.3f fraction of PARIS support is %.3f' % (fraction, threshold), file = self.ostream)
        
        self.stems = []
        for k in k_range:
            total_n_stems, passed_n_stems = 0, 0
            # Scan the forward sequence
            for l in range(self.rna.N):
                # Scan the reverse sequence
                for r in range(self.rna.N):
                    # Reserve the length for forming loops
                    if r - l >= k + Utils.MIN_LOOP_LEN:
                        left_part = self.rna.seq[l:l+k]
                        right_part = self.rna.seq[r:r+k]
                        if self.__check_stem(left_part, right_part, k):
                            support = np.mean(self.PARIS_support[l:l+k, r:r+k])
                            if support > threshold:
                                self.stems.append((l, r, k, support))
                                passed_n_stems += 1
                            total_n_stems += 1
            print('--- %d %d-mers can perfectly match (%d have high PARIS supports)' 
                  % (total_n_stems, k, passed_n_stems), file = self.ostream)
        print('--- %d stems are collected with high PARIS supports\n' % len(self.stems), file = self.ostream)
        
        
    def __remove_redundant_stems(self):
        clean_stems = []
        # Sort stems by length (from long to short)
        stems = sorted(self.stems, key = lambda x: - x[2])
        for stm in stems:
            l, r, k, support = stm
            passed = True
            # Check if this stem is covered by collected stems
            for l_, r_, k_, support_ in clean_stems:
                if l >= l_ and l+k <= l_+k_ and r >= r_ and r+k <= r_+k_:
                    passed = False
            if passed:
                clean_stems.append(stm)
        print('--- %d stems are retained after eliminating redundant stems' % len(clean_stems), file = self.ostream)
        self.stems = clean_stems
        
        
    def __get_hard_constraints(self):
        # Single stems as hard constraints
        all_combination_indices = [(i,) for i in range(len(self.stems))]
        # Double stems as hard constraints
        all_combination_indices.extend(list(itertools.combinations([i for i in range(len(self.stems))], 2)))

        print('--- %d combinations of stems need to check if is compatible' % len(all_combination_indices), file = self.ostream)
        self.hard_constraints = []
        for combination in all_combination_indices:
            # Representing dot-bracket notation use 0, 1, -1, where 0:., 1:(, -1:)
            structure_vector = np.zeros(self.rna.N, dtype = np.byte) 
            compatible = True
            for i in combination:
                l, r, k, support = self.stems[i]
                
                # Check if this stems overlaps with included stem
                if np.sum(np.abs(structure_vector[l:l+k])) + np.sum(np.abs(structure_vector[r:r+k])) != 0:
                    compatible = False; break
                
                # Ensure pseudoknot-free by the following rules:
                # 1) the prefix summation of the inner vector enclosed by the stem should be always non-negative
                close_stem_checker = 0
                for i in range(l+k, r):
                    close_stem_checker += structure_vector[i]
                    if close_stem_checker < 0:
                        compatible = False; break
                # 2) the summation of the inner vector enclosed by the stem should be zero
                if close_stem_checker != 0:
                    compatible = False; break
                
                if compatible:
                    structure_vector[l:l+k] = 1
                    structure_vector[r:r+k] = -1
                else:
                    break
            
            # Add the compatible stem as new hard constraint
            if compatible:
                structure_string = ['.' for _ in range(self.rna.N)]
                for i in range(self.rna.N):
                    if structure_vector[i] == 1:
                        structure_string[i] = '('
                    elif structure_vector[i] == -1:
                        structure_string[i] = ')'
                structure_string = ''.join(structure_string)
                self.hard_constraints.append(structure_string)
                
        print('--- %d compatible stems are reserved as hard constraints\n' % len(self.hard_constraints), file = self.ostream)
        
    
    def __proceed_to_localmin(self):
        print('--- Perform the constrained RNA folding algorithm', file = self.ostream)
        self.localmin_structures = []
        for i in range(len(self.hard_constraints)):
            if (i+1) % (len(self.hard_constraints) // 10) == 0:
                print('--- %d structures completed ...' % (i+1), file = self.ostream)
            # Thanks for ViennaRNA package!
            fc = RNA.fold_compound(self.rna.seq)
            fc.hc_add_from_db(self.hard_constraints[i], RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
            self.localmin_structures.append(fc.mfe()[0])
    
    
    def __cluster_structures(self, C):
        # Prepare base pair distance
        print('--- Prepare base pair distance matrix', file = self.ostream)
        M = len(self.localmin_structures)
        distance_matrix = np.zeros((M, M))
        for i in range(M):
            for j in range(i, M):
                # Thanks for the ViennaRNA package again!
                distance = RNA.bp_distance(self.localmin_structures[i], self.localmin_structures[j])
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
        
        # Perform clustering
        print('--- Perform agglomerative clustering', file = self.ostream)
        clustering_model = AgglomerativeClustering(n_clusters = C, affinity = 'precomputed', linkage = 'average')
        clustering_model.fit(distance_matrix)
        cluster = {}
        for i in range(len(clustering_model.labels_)):
            label = clustering_model.labels_[i]
            if label in cluster:
                cluster[label].append(i)
            else:
                cluster[label] = [i]
                
        # Retain the structure with lowest free energy in each cluster
        self.cluster_leaders = []
        for label in cluster:
            min_energy, min_structure = np.inf, ''
            for i in cluster[label]:
                energy = RNA.eval_structure_simple(self.rna.seq, self.localmin_structures[i])
                if energy < min_energy:
                    min_energy = energy
                    min_structure = self.localmin_structures[i]
            self.cluster_leaders.append(min_structure) 
    
    
    def generate(self, C, fraction, k_range):
        print('=== Generate candidate structures ===', file = self.ostream)
        print('+++ %d structures expected' % C, file = self.ostream)
        print('+++ Stems with PARIS support higher than %.3f of all PARIS support are collected' % fraction, file = self.ostream)
        print('+++ The range of length k for stems is ', k_range, file = self.ostream)
        self.__collect_stems(fraction, k_range)
        self.__remove_redundant_stems()
        self.__get_hard_constraints()
        self.__proceed_to_localmin()
        if len(self.localmin_structures) < C:
            self.candidate_structures = self.localmin_structures
        else:        
            self.__cluster_structures(C)
            self.candidate_structures = self.cluster_leaders
        print('--- Eventually %d structures are generated as candidates' % len(self.candidate_structures), file = self.ostream)
        print('=== Finished! === \n\n', file = self.ostream)
        
        
        
class Picking:
    
    def __init__(self, rna, PARIS_support, candidate_structures, ostream = sys.stdout):
        self.rna = rna
        self.PARIS_support = PARIS_support
        self.candidate_structures = candidate_structures
        self.ensemble = None
        self.ostream = ostream
        
        # Preprocess the tensor of all candidate structures
        dummy_ensemble = Structure_Ensemble(self.rna)
        dummy_ensemble.load(self.candidate_structures, np.ones(len(self.candidate_structures)))
        self.X = dummy_ensemble.get_tensor()
        # Reshape the base pairing matrices of candidate structures
        self.X_sparse = csc_matrix(self.X.reshape((self.rna.N ** 2, len(self.candidate_structures))))
        # Reshape the matrix of PARIS support
        self.Y = self.PARIS_support.reshape((self.rna.N ** 2))
        
        
    def __LASSO_filter(self, K):
        print('--- Eliminate redundant structures using LASSO regression', file = self.ostream)
        # Not necessary to remove redundancy
        if len(self.candidate_structures) <= 50:
            return [c for c in range(len(self.candidate_structures))]
                
        gamma = 1
        while gamma > Utils.EPS:
            model = Lasso(alpha = gamma, positive = True)
            model.fit(self.X_sparse, self.Y)
            nz_indices = np.nonzero(model.coef_)[0]
            print('--- Try gamma = %.2e, and get %d non-zero coefficients' % (gamma, len(nz_indices)), file = self.ostream)
            if len(nz_indices) > K:
                break
            gamma /= 2
            
        if gamma <= Utils.EPS:
            raise ValueError('Error: LASSO regression failed.')
        
        print('--- %d non-zero coefficients are retained (end up with gamma = %.2e)' % (len(nz_indices), gamma), file = self.ostream)
        return nz_indices
        
        
    def pick(self, K):
        print('=== Pick the optimal ensemble with %d structures ===' % K, file = self.ostream)
    
        # Preprocess PARIS support and free enery of each candidate structure
        fc = RNA.fold_compound(self.rna.seq, RNA.md()); fc.pf()
        max_logprob, min_logprob = -np.inf, np.inf
        max_support, min_support = -np.inf, np.inf
        for c in range(len(self.candidate_structures)):
            # the mean PARIS support of base pairs in the structure
            support = np.mean(self.PARIS_support[self.X[:,:,c] > 0]) 
            max_support = max(max_support, support)
            min_support = min(min_support, support)
            # the log-probability from the Boltzmann distribution
            logprob = np.log(fc.pr_energy(RNA.eval_structure_simple(self.rna.seq, self.candidate_structures[c])))
            max_logprob = max(max_logprob, logprob)
            min_logprob = min(min_logprob, logprob)
        
        # Compute the parameter lambda for the exponential distribution as the likelihood
        lamb = (max_logprob - min_logprob) / (max_support - min_support)
        print('--- The lambda for the exponential distribution is %.3f' % lamb, file = self.ostream)
        
        # Get indices of structures with non-zero coefficients in LASSO regression
        nonredundant_indices = self.__LASSO_filter(K)
        
        # Enumerate all the K combinations 
        combination_indices = list(itertools.combinations([_ for _ in range(len(nonredundant_indices))], K))
        candidate_ensembles = []
        for indices in combination_indices:
            subset_indices = [nonredundant_indices[i] for i in indices]

            # Compute the proportion alpha using Linear regression
            cur_alpha = lsq_linear(self.X[:,:,subset_indices].reshape((self.rna.N ** 2, K)), self.Y, bounds = (0, np.inf)).x
            
            cur_ensemble = Structure_Ensemble(self.rna)
            cur_ensemble.load([self.candidate_structures[i] for i in subset_indices], np.array(cur_alpha))
            candidate_ensembles.append(cur_ensemble)
        print('--- %d candidate ensembles are generated' % len(candidate_ensembles), file = self.ostream)
            
        # Pick the optimal ensemble using a Bayesian model
        max_log_posterior = -np.inf
        for t in range(len(candidate_ensembles)):
            log_prior = np.log(fc.pr_energy(candidate_ensembles[t].get_mean_free_enery()))
            log_likelihood = expon.logpdf(max_support - candidate_ensembles[t].get_mean_PARIS_support(self.PARIS_support), scale = 1.0/lamb)
            if log_prior + log_likelihood > max_log_posterior:
                max_log_posterior = log_prior + log_likelihood
                self.ensemble = candidate_ensembles[t]
                
        print('=== Finished ===\n\n', file = self.ostream)