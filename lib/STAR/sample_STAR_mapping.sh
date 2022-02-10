#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=PARIS_STAR_align
#SBATCH --output=PARIS_STAR_align_%J.log

source /rhome/jianyu/.bashrc
WORKDIR=~/bigdata/EnsembleRNA/proc_PARIS

SEQDIR=~/bigdata/EnsembleRNA/data_PARIS/huh7
#INDEXNAME=train_index
#ALIGNNAME=huh7_train_align
INDEXNAME=test_index
ALIGNNAME=huh7_test_align
#INDEXNAME=known_index
#ALIGNNAME=huh7_known_align
#INDEXNAME=zikv_index
#ALIGNNAME=huh7_zikv_align

# please sumbit 2 individual jobs 59_24_ and 59_72_ for huh7
#BATCHNAME=59_24_
#BATCHNAME=59_72_
BATCHNAME=766_24_
#BATCHNAME=766_72_


#SEQDIR=~/bigdata/EnsembleRNA/data_PARIS/hek293t
#INDEXNAME=train_index
#ALIGNNAME=hek293t_train_align
#INDEXNAME=test_index
#ALIGNNAME=hek293t_test_align
#INDEXNAME=known_index
#ALIGNNAME=hek293t_known_align

# please sumbit 3 individual jobs SRR2814763.1_ ~ SRR2814765.1_ for hek293t
#BATCHNAME=SRR2814763.1_
#BATCHNAME=SRR2814764.1_
#BATCHNAME=SRR2814765.1_

echo ${BATCHNAME}

STAR --runThreadN 16 \
     --genomeDir ${WORKDIR}/${INDEXNAME} \
     --readFilesIn ${SEQDIR}/${BATCHNAME}trim_nodup.fastq \
     --outFileNamePrefix ${WORKDIR}/${ALIGNNAME}/${BATCHNAME} \
     --outSAMtype BAM SortedByCoordinate \
     --outSJfilterReads Unique \
     --outSJfilterOverhangMin 5 5 5 5 \
     --outSJfilterCountUniqueMin 1 1 1 1 \
     --outSJfilterCountTotalMin 1 1 1 1 \
     --outSJfilterDistToOtherSJmin 0 0 0 0 \
     --scoreGapNoncan 0 \
     --scoreGapGCAG 0 \
     --scoreGapATAC 0 \
     --limitBAMsortRAM 10000000000 \
     --alignIntronMin 10
     #--twopassMode Basic \
     #--limitSjdbInsertNsj 4000000 \
     #--limitOutSJcollapsed 10000000 \
     #--limitIObufferSize 300000000