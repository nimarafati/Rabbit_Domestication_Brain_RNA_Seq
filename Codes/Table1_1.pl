#!/usr/bin/perl
#perl Table1_1.pl Amy Hyp Hipp ParTemp
#For GREAT analysis

use strict;
use warnings;

my %hash;
open(IN,"../Data/raw_data/Ens76_Orthologs.txt");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  next unless($_=~/ortholog_one2one/);
  $hash{$li[1]}=$li[2];
}
close(IN);

my %chr;
my %st;
my %en;
open(IN,"../Data/raw_data/Ens76_HumanGeneCoordinates.txt");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  $chr{$li[0]}=$li[2];
  $st{$li[0]}=$li[3];
  $en{$li[0]}=$li[4];
}
close(IN);

for my $reg (@ARGV){
  open(OUT,">../Data/GREAT_$reg\_DEGs_Ens76HumanGeneCoordinates.bed");
  open(IN,"../Data/DEGs_$reg\_edgeR.tsv");
  while(<IN>){
    chomp;
    my @li=split(/\s+/,$_);
    next if($_=~/^logFC/);
    if(defined $hash{$li[0]}){
      my $gene=$hash{$li[0]};
      if(defined $chr{$gene} and $chr{$gene}!~/([A-Z])/){
        print OUT "chr$chr{$gene}\t$st{$gene}\t$en{$gene}\n"
      }
    }
  }
  close(IN);
  close(OUT);
  open(OUT,">../Data/GREAT_$reg\_background_Ens76HumanGeneCoordinates.bed");
  open(IN,"../Data/Allgenes_$reg\_edgeR.tsv");
  while(<IN>){
    chomp;
    next if($_=~/^logFC/);
    my @li=split(/\s+/,$_);
    if(defined $hash{$li[0]}){
      my $gene=$hash{$li[0]};
      if(defined $chr{$gene} and $chr{$gene}!~/([A-Z])/){
        print OUT "chr$chr{$gene}\t$st{$gene}\t$en{$gene}\n"
      }
    }
  }
  close(IN);
  close(OUT);
}
