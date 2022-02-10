#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=PARIS_STAR_index
#SBATCH --output=PARIS_STAR_index_%J.log

source /rhome/jianyu/.bashrc
WORKDIR=~/bigdata/EnsembleRNA/proc_PARIS

#REFFILE=~/bigdata/EnsembleRNA/ref_test/all.fa
#INDEXNAME=test_index
#SAPAR=11
#CHRBIT=10

#REFFILE=~/bigdata/EnsembleRNA/ref_test/ZIKV/ZIKV_59.fa
#INDEXNAME=zikv_index
#SAPAR=6
#CHRBIT=6

#REFFILE=~/bigdata/EnsembleRNA/ref_train/bpRNA_Homo.S._curated_nonsim.fa
#INDEXNAME=train_index
#SAPAR=6
#CHRBIT=7

REFFILE=~/bigdata/EnsembleRNA/ref_known/Rfam_H.S._nonsim.fa
INDEXNAME=known_index
SAPAR=7
CHRBIT=7

STAR --runMode genomeGenerate\
     --genomeDir ${WORKDIR}/${INDEXNAME}\
     --genomeFastaFiles ${REFFILE}\
     --genomeSAindexNbases ${SAPAR}\
     --genomeChrBinNbits ${CHRBIT}