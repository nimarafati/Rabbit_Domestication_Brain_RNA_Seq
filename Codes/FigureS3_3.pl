#!/usr/bin/perl
use strict;
use warnings;

my %gene;
my %tid;
my %st;
my %en;
open(IN,"../Data/raw_data/Ens76_Rabbit_genename_ID_chr_start_end_strand.txt");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  my $chr=$li[3];
  $gene{$chr}{$li[1]}=$li[0];
  $tid{$chr}{$li[1]}=$li[2];
  $st{$chr}{$li[1]}=$li[4];
  $en{$chr}{$li[1]}=$li[5];
}
close(IN);

open(OUT,">../Data/Deltas_1SNPin50kb_GeneName.tsv");
open(IN,"../Data/Deltas_1SNPin50kb_2.txt");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  my $chr;
  if($li[0]=~/chr(\S+)/){
    $chr=$1;
  }
  foreach my $key (sort {$a cmp $b} keys %{$tid{$chr}}){
    if($li[1] > $st{$chr}{$key}-100000 and $li[1] < $en{$chr}{$key}+100000){
      print OUT "$_\t$key\t$tid{$chr}{$key}\t$gene{$chr}{$key}\t$st{$chr}{$key}\t$en{$chr}{$key}\n";
    }
  }
}
close(IN);
