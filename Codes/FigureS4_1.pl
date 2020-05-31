#!/usr/bin/perl
use strict;
use warnings;

my %tss;
open(IN,"../Data/raw_data/Oryctolagus_cuniculus.OryCun2.0.76.gtf");
while(<IN>){
  chomp;
  my @li=split(/\t/,$_);
  next unless($li[2] eq "gene");
  my $chr="chr$li[0]";
  my $geneid;
  if($li[8]=~/gene_id \"(\S+)\"\;/){
    $geneid=$1;
  }
  if($li[6] eq "+"){
    $tss{$chr}{$geneid}=$li[3];
  }else{
    $tss{$chr}{$geneid}=$li[4];
  }
}

open(OUT,">../Data/All_dAF_TSS100kb.txt");
print OUT "Chr\tPos\tdAF\tBin\tGenes\n";
open(IN,"../Data/raw_data/All_dAF_0.5over");
while(<IN>){
  chomp;
  my @li=split(/\t/,$_);
  next if($_=~/^CHR/);
  my $dafr=int(10*$li[2]+0.5)/10;
  my $dafr2;
  if($dafr > $li[2]){
    $dafr2=$dafr;
    $dafr=$dafr-0.05;
  }else{
    $dafr2=$dafr+0.05;
  }
  print OUT "$li[0]\t$li[1]\t$li[2]\t$dafr\-$dafr2\t";
  my $count=0;
  my @temp;
  foreach my $key (keys %{$tss{$li[0]}}){
    if($tss{$li[0]}{$key}-100000 < $li[1] and $tss{$li[0]}{$key}+100000 > $li[1]){
      push @temp,"$key\($tss{$li[0]}{$key}\)";
      $count++;
    }
  }
  if($count > 0){
    my $temp=join(",",@temp);
    print OUT "$temp\n";
  }else{
    print OUT "-\n";
  }
}
