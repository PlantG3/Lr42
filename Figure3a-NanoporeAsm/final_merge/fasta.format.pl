#!/usr/bin/perl -w
#
# ====================================================================================================
# File: fasta.format.pl
# Author: Sanzhen Liu
# Date: 12/1/2016
# ====================================================================================================

use strict;
use warnings;
use Getopt::Long;

my ($seq, $size, $seq_name, %seqhash, %sizehash);
my ($sort, $bpperline, $removeN, $based, $seq_prefix);
my $help;

# default
$sort = "n";
$based = "l";
$bpperline = 80;
$seq_prefix = "";

sub prompt {
    print <<EOF;
    Usage: perl fasta.format.pl <Input Fasta Files> [options]
    --bpperline: bp per line ($bpperline)
    --sort     : sort sequences (i=increasing; d=decreasing; n=no sorting,default)
    --based    : sequence name (n) or sequence length (l, default)
    --rmN      : remove sequences with N if specified
    --seqprefix: prefix to add to each sequence name ($seq_prefix)
    --help     : help information
EOF
exit;
}
# read the parameters:
&GetOptions("bpperline=i" => \$bpperline,
            "sort=s"      => \$sort,
			"based=s"     => \$based,
			"rmN"         => \$removeN,
			"seqprefix=s" => \$seq_prefix,
			"help|h"      => \$help) || &prompt;

if ($help) { &prompt; }
if (@ARGV<1) { &prompt; }

if ($sort ne "i" and $sort ne "d" and $sort ne "n") {
	print "Only three values can be assigned for --sort: i (increasing), d (decreasing), and n (no)\n";
	exit;
}

###
### main
###

my @ori_seqnames;

open(IN, $ARGV[0]) || die;
while (<IN>) {
	chomp;
	if (/^>(\S+)/) {
		if (defined $seq_name) {
			push(@ori_seqnames, $seq_name);
			$seqhash{$seq_name} = $seq;
			$sizehash{$seq_name} = $size;
		}
    	$seq_name = $1;
		$size = 0;
		$seq = "";
 	 } else {
		$seq .= $_;
		$size += length($_);
	}
}
# last element:
push(@ori_seqnames, $seq_name);
$seqhash{$seq_name} = $seq;
$sizehash{$seq_name} = $size;
close IN;

my @allseqnames;
if ($sort eq "i") {
	if ($based eq "l") {
		@allseqnames = sort {$sizehash{$a} <=> $sizehash{$b}} keys %sizehash;
	} elsif ($based eq "n") {
		@allseqnames = sort {$a cmp $b} keys %sizehash;
	} 
} elsif ($sort eq "d") {
	if ($based eq "l") {
		@allseqnames = sort {$sizehash{$b} <=> $sizehash{$a}} keys %sizehash;
	} elsif ($based eq "n") {
		@allseqnames = sort {$b cmp $a} keys %sizehash;
	}
} else {
	@allseqnames = keys %sizehash;
	@allseqnames = @ori_seqnames;
}

# sequence output
foreach my $eachseq (@allseqnames) {
	if (defined $removeN and $seqhash{$eachseq} =~ /N/) {
		print STDERR "$eachseq was removed due to the presence of N\n";
	} else {
		printf("%s%s%s\n", ">", $seq_prefix, $eachseq);
		&format_print($seqhash{$eachseq}, $bpperline);
	}
}

###
### function for formatted output:
###
sub format_print {
	my ($inseq, $formatlen) = @_;
	while (my $chunk = substr($inseq, 0, $formatlen, "")) {
		print "$chunk\n";
	}
}

