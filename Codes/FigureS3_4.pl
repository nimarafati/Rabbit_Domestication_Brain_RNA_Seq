#!/usr/bin/perl
use strict;
use warnings;

my %test;
#remove redundant descriptions for a gene
open(OUT,">../Data/Deltas_1SNPin50kb_GeneName2.tsv");
open(IN,"../Data/Deltas_1SNPin50kb_GeneName.tsv");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  $test{$li[2]}="$li[3]\t$li[4]";
}
foreach my $key (sort {$a cmp $b} keys %test){
  print OUT "$key\t$test{$key}\n";
}
close(IN);
close(OUT);

my %hash;
my %hasht;
open(IN,"../Data/raw_data/Ens76_Orthologs.txt");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  next if($_=~/^Ensembl/);
  if(defined $li[2]){
    if(defined $hash{$li[0]} and $hash{$li[0]} ne $li[2]){
      $hash{$li[0]}="$hash{$li[0]}\,$li[2]";
    }else{
      $hash{$li[0]}=$li[2];
    }
    $hasht{$li[0]}=$li[3];
  }
}
close(IN);

open(OUT,">../Data/Deltas_1SNPin50kb_GeneName2_ortho.tsv");
open(IN,"../Data/Deltas_1SNPin50kb_GeneName2.tsv");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  if(defined $hash{$li[0]}){
    print OUT "$_\t$hash{$li[0]}\t$hasht{$li[0]}\n";
  }else{
    print OUT "$_\n";
  }
}
close(IN);
