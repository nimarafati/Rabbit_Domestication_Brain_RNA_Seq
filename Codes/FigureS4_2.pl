#!/usr/bin/perl
use strict;
use warnings;

my %hash;
open(IN,"../Data/raw_data/Ens76_GeneCoordinates.txt");
while(<IN>){
  chomp;
  my @li=split(/\t/,$_);
  $hash{$li[4]}=$li[3];
}
close(IN);

my %degs;
my @regs=("Amy","Hyp","Hipp","ParTemp");
foreach my $reg (@regs){
  open(IN,"../Data/DEGs_$reg\_edgeR.tsv");
  while(<IN>){
    chomp;
    my @li=split(/\t/,$_);
    next if($_=~/^TranscriptID/);
    my $deg=$hash{$li[0]};
    $degs{$reg}{$deg}=1;
  }
}

open(OUT,">../Data/All_dAF_TSS100kb_DEGs.txt");
print OUT "Chr\tPos\tdAF\tBin\tGenes\tDEG_Amy\tDEG_Hyp\tDEG_Hipp\tDEG_ParTemp\n";
open(IN,"../Data/All_dAF_TSS100kb.txt");
while(<IN>){
  chomp;
  my @li=split(/\t/,$_);
  next if($_=~/^Chr/);
  next if($li[2] eq "gene");
  print OUT "$_";
  if($li[4]=~/ENSOCUG/){
    my @li2=split(/,/,$li[4]);
    foreach my $reg (@regs){
      my $temp2="No";
      foreach my $gene (@li2){
        my $temp;
        if($gene =~ /(\S+)\((\S+)\)/){
          $temp=$1;
        }
        $temp2="Yes" if(defined $degs{$reg}{$temp});
      }
      print OUT "\t$temp2";
    }
    print OUT "\n";
  }else{
    print OUT "\t-\t-\t-\t-\n";
  }
}
