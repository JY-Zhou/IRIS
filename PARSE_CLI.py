import argparse
from PARSE.Info_RNA import *
from PARSE.PARSE_Core import *
from Bio import SeqIO


VERSION = 'v1.0'
def run_PARSE(args):
    # Load reference sequence and PARIS reads as an Info_RNA object
    print('--- Loading inputs ...')
    print('  --- RNA sequence : %s' % args.input_fasta)
    print('  --- PARIS reads: %s' % args.input_bam)
    print('')
    for ref_rna in SeqIO.parse(args.input_fasta, 'fasta'):
        concerning_rna = Info_RNA(ref_rna.id, ref_rna.seq)
    concerning_rna.PARIS = Info_PARIS(args.input_bam)
    
    # Perform scoring of PARSE
    PARSE_scoring = Scoring(concerning_rna)
    PARSE_scoring.compute_PARIS_support()
    
    PARSE_generating = Generating(concerning_rna, PARSE_scoring.PARIS_support)
    PARSE_generating.generate(args.C, args.F, args.krange)
    
    PARSE_picking = Picking(concerning_rna, PARSE_scoring.PARIS_support, PARSE_generating.candidate_structures)
    PARSE_picking.pick(args.K)
    
    print('--- The predicted ensemble is saved in %s' % args.output)
    PARSE_picking.ensemble.save(args.output)
    print('')
    

parser = argparse.ArgumentParser(prog = 'python PARSE_CLI.py', description='PARSE %s: A method for predicting PARIS-based RNA secondary structure ensembles' % VERSION)
parser.add_argument('-i', metavar = 'FASTA file', action = 'store', \
                    type=str, nargs = '?', required = True,\
                    help='the input RNA sequence in FASTA format', dest = 'input_fasta')

parser.add_argument('-r', metavar = 'BAM file', action = 'store', \
                    type=str, nargs = '?', required = True,  \
                    help='the input mapped PARIS reads in BAM format', dest = 'input_bam')

parser.add_argument('-o', metavar = 'output file', action = 'store', \
                    type=str, nargs = '?', required = True, \
                    help='the output file of the predicted ensemble', dest = 'output')

parser.add_argument('-K', metavar = 'int', action = 'store', \
                    type=int, nargs = '?', required = True,\
                    help='the expected number of representative structures in the ensemble', dest = 'K')

parser.add_argument('-C', metavar = 'int', action = 'store', \
                    type=int, nargs = '?', default = 100,\
                    help='the number of candidate structures generated (default: 100)', dest = 'C')

parser.add_argument('-F', metavar = 'float', action = 'store', \
                    type=float, nargs = '?', default = 0.75, \
                    help='the fraction threshold for filtering stems (default: 0.75)', dest = 'F')

parser.add_argument('--krange', metavar = '', action = 'store', \
                    type=int, nargs = '+', default = [4, 5, 6, 7],\
                    help='the range of length k of stems (default: 4 5 6 7)', dest = 'krange')

args = parser.parse_args()
run_PARSE(args)
