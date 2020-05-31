#!/usr/bin/perl
use strict;
use warnings;

my %hash;
open(IN,"../Data/Deltas_1SNPin50kb_GeneName2.tsv");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  $hash{$li[1]}=1;
}
close(IN);

open(OUT,">../Data/DEGs_edgeR_dAF.tsv");
open(IN,"../Data/Allgenes_edgeR.tsv");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  if($_=~/^TranscriptID/){
    print OUT "$_\tdAF\n";
  }else{
    if(defined $hash{$li[0]}){
      print OUT "$_\tYes\n";
    }else{
      print OUT "$_\tNo\n"
    }
  }
}
close(IN);
