#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-00:00:00
in=1o_TA2450nanoWGS2lr42locus
out=2o_TA2450nanoWGS2lr42locus.parse

perl ~/scripts/sam/samparser.minimap2.pl \
	--input ${in}.parse.sam \
	--identical 12000 \
	--mismatches 20 100 --tail 95 100 \
	--gap 1000 --mappingscore 30 \
	> ${out}.sam

# sam2fq
module load SAMtools
samtools fastq ${out}.sam > ${out}.fq

