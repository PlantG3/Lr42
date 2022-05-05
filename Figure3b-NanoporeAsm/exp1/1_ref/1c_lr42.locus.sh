#!/bin/bash
awk '$1 == "Chr1D" && $2>8000000 && $3<10000000' /bulk/guifanglin/project/wheat/14-NLRs/4-genome_walking/4-cgrd/cgrdo_TA2450_TA10132.segments.txt > TA2450_TA10132.CGRD.partial.txt
bedtools getfasta -bed 1i_conserved.Lr42.locus.bed -fi /bulk/guifanglin/project/wheat/0_reference/organized/Aetv4.0/Aetv4.0.full.fasta -fo 1o_conserved.Lr42.locus.fasta
cat 1o_conserved.Lr42.locus.fasta TA2450_Lr42.gene.locus.fas > 1o_conserved.Lr42+gene.fasta
minimap2 -x map-ont -d 1o_conserved.Lr42+gene.mmi 1o_conserved.Lr42+gene.fasta

