#!/usr/bin/perl
#perl TableS2_1.pl Amy Hyp Hipp ParTemp
#For DAVID analysis

use strict;
use warnings;

my %hash;
open(IN,"../Data/raw_data/Ens76_Orthologs.txt");
while(<IN>){
  chomp;
  next if($_=~/^Ensembl/);
  my @li=split(/\s+/,$_);
  next unless defined($li[3]);
  if($li[3]=~/one2one/){
    $hash{$li[1]}=$li[2];
  }
}
close(IN);

my %hash2;
for my $reg (@ARGV){
  open(OUT,">../Data/DAVID_$reg\_background_HumanEns76.txt");
  open(IN,"../Data/Allgenes_$reg\_edgeR.tsv");
  while(<IN>){
    chomp;
    next if($_=~/^Transcript/);
    my @li=split(/\s+/,$_);
    if(defined $hash{$li[0]}){
      print OUT "$hash{$li[0]}\n";
    }
  }
  close(IN);
  close(OUT);
  open(OUT,">../Data/DAVID_$reg\_DEGs_HumanEns76.txt");
  open(IN,"../Data/DEGs_$reg\_edgeR.tsv");
  while(<IN>){
    chomp;
    next if($_=~/^logFC/);
    my @li=split(/\s+/,$_);
    if(defined $hash2{$li[0]}){
      $hash2{$li[0]}++;
    }else{
      $hash2{$li[0]}=1;
    }
    if(defined $hash{$li[0]}){
      print OUT "$hash{$li[0]}\n";
    }
  }
  close(IN);
  close(OUT);
}
