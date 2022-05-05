#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3G
#SBATCH --time=1-00:00:00
ref=/bulk/liu3zhen/research/projects/Lr42/Lr42_nanopore/1_ref/1o_conserved.Lr42+gene.mmi
fq=/bulk/guifanglin/raw/TA/TA2450-nanopore/guppy4.2.2/TA2450wgs.np.gp422.fq
out=1o_TA2450nanoWGS2lr42locus
minimap2 -ax map-ont -t 24 $ref $fq > ${out}.sam 

perl  ~/scripts/sam/samparser.minimap2.pl \
	--input ${out}.sam \
	--identical 3000 \
	--mismatches 20 100 --tail 90 100 \
	--gap 1000 --mappingscore 30 \
	> ${out}.parse.sam

