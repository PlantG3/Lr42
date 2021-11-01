#!/bin/bash
bedtools getfasta -bed 1i_conserved.Lr42.locus.bed -fi /bulk/guifanglin/project/wheat/0_reference/organized/Aetv4.0/Aetv4.0.full.fasta -fo 1o_conserved.Lr42.locus.fasta
#cat 1o_conserved.Lr42.locus.fasta TA2450_Lr42.gene.locus.fas > 1o_conserved.Lr42+gene.fasta
minimap2 -x map-ont -d 1o_conserved.Lr42.locus.mmi 1o_conserved.Lr42.locus.fasta

