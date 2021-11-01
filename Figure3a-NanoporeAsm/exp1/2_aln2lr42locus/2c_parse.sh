#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-00:00:00
ref=/bulk/liu3zhen/research/projects/Lr42/Lr42_nanopore/1_ref/1o_conserved.Lr42+gene.mmi
out=1o_TA2450nanoWGS2lr42locus
perl ~/scripts/sam/samparser.minimap2.pl \
	--input ${out}.parse.sam \
	--identical 8000 \
	--mismatches 16 100 --tail 95 100 \
	--gap 100 --mappingscore 30 \
	> ${out}.parse.sam

