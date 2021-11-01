#!/bin/bash
# get sequences from two assemblies
perl fasta.format.pl --seqprefix exp1_ ../../exp1/3_flye/assembly.fasta > exp1.ctg.fasta
perl fasta.format.pl --seqprefix exp2_ ../../exp2/3_flye/assembly.fasta > exp2.ctg.fasta
cat exp1.ctg.fasta exp2.ctg.fasta > all.fasta

# manually examine the overlap and make a reorganizing table 1i_merge.info
perl fasta.reorganiz.pl --fasta all.fasta --table 1i_merge.info > Lr42locus.draft.fasta

# nucmer to compare
nucmer Lr42locus.draft.fasta ../../exp1/1_ref/1o_conserved.Lr42.locus.fasta

