import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import rgb2hex, ListedColormap
import matplotlib as mpl
from PIL import Image

from sklearn.manifold import MDS
import scipy.spatial

import RNA


class Utils:
    
    # Important constants
    EPS = np.finfo(np.float32).eps

    MIN_LOOP_LEN = 4
    VALID_BP = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']

    VBLUE = '#1ba1e2' #{27,161,226}
    VGREEN = '#339933' #{51,153,51}
    VRED = '#e51400' #{229,20,0}
    VORANGE = '#f09609' #{240,150,9}
    VPINK = '#d80073' #{216,0,115}
    VCYAN = '#00aba9' #{0,171,169}
    VLGRAY = '#666666'

    # Create a color map for PARIS support
    log_step = (10 - np.logspace(0, 1, num = 256, base = 10)[::-1]) / 20
    black2red_cmap = ListedColormap(cm.get_cmap('hot')(log_step))
    
    white2red_cmap = ListedColormap([[1-i/256*27/256, 1-i/256*236/256, 1-i/256, 1] for i in range(256)])

    

    @staticmethod
    def get_pairs(structure):
        pairs = []
        # Count normal pairs
        ptable = RNA.ptable(structure)
        for i in range(1, len(ptable)):
            if ptable[i] > 0 and i < ptable[i]:
                pairs.append((i-1, ptable[i]-1))
        # Count level-1 pseudoknot pairs
        ss = structure.replace('(','.').replace(')','.').replace('[', '(').replace(']', ')')
        ptable = RNA.ptable(ss)
        for i in range(1, len(ptable)):
            if ptable[i] > 0 and i < ptable[i]:
                pairs.append((i-1, ptable[i]-1))
        # Count level-2 pseudoknot pairs
        ss = structure.replace('(','.').replace(')','.').replace('{', '(').replace('}', ')')
        ptable = RNA.ptable(ss)
        for i in range(1, len(ptable)):
            if ptable[i] > 0 and i < ptable[i]:
                pairs.append((i-1, ptable[i]-1))
        return pairs
    

    @staticmethod
    def plot_triangle_matrix(ax, matrix, cmap):
        x, y, c, N = [], [], [], matrix.shape[0]
        for i in range(N):
            for j in range(i+1, N):
                x.append((i+j)/2+1)
                y.append((j-i)/2)
                c.append(matrix[i, j])
        ax.scatter(x, y, s = 1, c = c, cmap = cmap)
        ax.plot([0, N/2, N, 0], [0, N/2, 0, 0], c = Utils.VLGRAY, linestyle = '--', linewidth = 1)
        ax.axis('equal')
        ax.set(frame_on = False)
        ax.axes.get_yaxis().set_visible(False)    
        
        
    @staticmethod
    def emphasize(score, baseline, operator):
        if operator == '>=':
            if score >= baseline:
                return '\x1b[31m%.1f\x1b[0m' % score
        if operator == '<=':
            if score <= baseline:
                return '\x1b[31m%.1f\x1b[0m' % score
        return '%.1f' % score
        
        
    @staticmethod
    def plot_structure_landscape(sequence, structures, PARIS_support, art_style):
        collection = []
        for label in structures:
            collection.extend(list(zip([label for _ in range(len(structures[label]))], structures[label])))

        # Get the base pairing distance matrix
        print('--- Computing the base-pairing distance matrix')
        bp_dist_mat = np.array([[RNA.bp_distance(x[1], y[1]) for x in collection] for y in collection])
        print('--- Base pair distance matrix: \n', bp_dist_mat)

        # Perform multidimensional scaling
        print('--- Embedding to 2D using multidimensional embedding (MDS)')
        embedded_res = MDS(n_components=2, dissimilarity='precomputed', random_state = 29).fit_transform(bp_dist_mat)

        # Divide the landscape using Voronoi diagram
        print('--- Partition using Voronoi diagram and color by free energy')
        bound_v = np.max(np.abs(embedded_res))
        dummy_v = bound_v * 10
        dummy_points = [[dummy_v, dummy_v], [-dummy_v, -dummy_v], [-dummy_v, dummy_v], [dummy_v, -dummy_v]]
        vor = scipy.spatial.Voronoi(np.append(embedded_res, dummy_points, axis = 0))

        embedded_structures = {label:[] for label in structures}
        for i in range(len(collection)):
            label = collection[i][0]
            embedded_structures[label].append(embedded_res[i])

        # Prepare colors for free energy
        collected_energy = [RNA.eval_structure_simple(sequence, collection[i][1]) for i in range(len(collection))]
        energy_cmap_norm = mpl.colors.Normalize(vmin = np.min(collected_energy), vmax = 0, clip=True)
        energy_cmap = cm.ScalarMappable(cmap = cm.Blues_r, norm = energy_cmap_norm)
        energy_cmap.set_array([])

        # Prepare colors for PARIS support
        supports = []
        for i in range(len(collection)):
            supports.append(np.mean([PARIS_support[l, r] for l, r in Utils.get_pairs(collection[i][1])]))
        PARIS_cmap_norm = mpl.colors.Normalize(vmin = 0, vmax = np.max(supports))
        PARIS_cmap = cm.ScalarMappable(cmap = cm.Greens, norm = PARIS_cmap_norm)
        PARIS_cmap.set_array([])


        score = [collected_energy, supports]
        colormap = [energy_cmap, PARIS_cmap]
        title = ['Free energy (kcal/mol)', 'PARIS support']
        ticks = [[-40, -30, -20, -10, 0], [0, 0.5, 1, 1.5]]
        for f in [0, 1]:
            # Plot the Voronoi diagram
            fig = scipy.spatial.voronoi_plot_2d(vor, show_points = False, show_vertices = False, line_width = 1)
            fig.set_size_inches((8, 8))
            pad = 5
            bound_l = np.min([x[0] for x in embedded_res]) - pad
            bound_r = np.max([x[0] for x in embedded_res]) + pad
            bound_d = np.min([x[1] for x in embedded_res]) - pad
            bound_u = np.max([x[1] for x in embedded_res]) + pad
            plt.xlim([bound_l, bound_r]), plt.ylim([bound_d, bound_u])
            plt.xticks([]), plt.yticks([])
            
            # Plot embedded points denoting structures
            for label in structures:
                points = embedded_structures[label]
                plt.scatter([p[0] for p in points], [p[1] for p in points], alpha = 1, label = label, \
                            c = art_style[label]['color'], \
                            marker = art_style[label]['marker'], \
                            s = art_style[label]['scale'] * 20, \
                            edgecolor = art_style[label]['edgecolor'], \
                            zorder = art_style[label]['zorder'] + 2)

            # Color each area occupied by a structure
            for i in range(len(vor.point_region) - len(dummy_points)):
                reg = vor.regions[vor.point_region[i]]
                polygon = []
                for v in reg:
                    polygon.append(vor.vertices[v])
                plt.fill(*zip(*polygon), alpha = 1, zorder = 0, color = colormap[f].to_rgba(score[f][i]))

            plt.legend(loc = 'upper left', fontsize = 12)
            cbar = plt.colorbar(colormap[f], orientation='horizontal', cax = fig.add_axes([0.17,0.1,0.7,0.02]), ticks = ticks[f])
            cbar.set_label(title[f], fontsize = 16)
            plt.show()
    

    
class Structure_Plotter:

    VARNA = 'java -cp ./lib/VARNAv3-92.jar fr.orsay.lri.varna.applications.VARNAcmd '
    
    def __init__(self, sequence, PARIS_support):
        self.sequence = sequence
        self.PARIS_support = PARIS_support
    
    def __decorate_pair(self, structure):
        decor_pairs = []
        PARIS_support_cmap_norm = mpl.colors.Normalize(0, np.max(self.PARIS_support))
        for l, r in Utils.get_pairs(structure):
            support = self.PARIS_support[l, r]
            thickness = 4
            color = rgb2hex(Utils.black2red_cmap(PARIS_support_cmap_norm(support)))
            decor_pairs.append((l+1, r+1, thickness, color))
        return decor_pairs
    
    
    def __trim_white_background(self, figure_path):
        img = Image.open(figure_path)
        img = img.convert("RGBA")
        pixdata = img.load()
        width, height = img.size
        for y in range(height):
            for x in range(width):
                if pixdata[x, y] == (255, 255, 255, 255):
                    pixdata[x, y] = (255, 255, 255, 0)
        img.save(figure_path + 'clean.png',  "PNG")
    
        
    def plot(self, structure, output_file, resolution = 10):
        cmd = Structure_Plotter.VARNA 
        cmd += ' -algorithm line -drawBase True -baseOutline "#000000" -baseName "#000000" -bpStyle none -bpIncrement 1 '
        cmd +=  '-resolution %d ' % resolution
        cmd += '-sequenceDBN "%s" ' % self.sequence
        cmd += '-structureDBN "%s" ' % ''.join(['.' for _ in range(len(self.sequence))])
        cmd += '-o "%s" ' % output_file
        pair_cmd = ''
        for pars in self.__decorate_pair(structure):
            pair_cmd += '(%d,%d):thickness=%d,color=%s;' % pars
        cmd += '-auxBPs "%s" ' % pair_cmd

        os.system(cmd)
        
        self.__trim_white_background(output_file)