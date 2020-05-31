#!/usr/bin/perl
#perl Table2.pl Amy Hyp Hipp ParTemp
use strict;
use warnings;
#sudo perl -MCPAN -e 'install Statistics::R'
use Statistics::R;

my %hash;
open(IN,"../Data/raw_data/All-Sweeps.bed");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  $hash{$li[0]}{$li[1]}=$li[2];
}
close(IN);
my %chr;
my %st;
my %en;
my %gene;
open(IN,"../Data/raw_data/Ens76_GeneCoordinates.txt");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  $chr{$li[4]}=$li[0];
  $st{$li[4]}=$li[1];
  $en{$li[4]}=$li[2];
  $gene{$li[4]}=$li[3];
}
close(IN);
print "Region\tDEGs_Sweep\tDEGs_NonSweep\tNonDEGs_Sweep\tNonDEGs_NonSweep\n";
for my $reg (@ARGV){
  my %co; my %gene;
  $co{de}{sw}=0; $co{de}{nsw}=0;
  $co{nde}{sw}=0; $co{nde}{nsw}=0;
  open(IN,"../Data/DEGs_$reg\_edgeR.tsv");
  while(<IN>){
    chomp;
    next if($_=~/^Transcript/);
    my @li=split(/\s+/,$_);
    $gene{$li[0]}=1;
    my $co=0;
    if(defined $chr{$li[0]}){
      my $chr=$chr{$li[0]};
      if(defined $hash{$chr}){
        foreach my $key (sort {$a <=> $b} keys %{$hash{$chr}}){
          next if($en{$li[0]}<$key-100000 or $st{$li[0]}>$hash{$chr}{$key}+100000);
          $co++;
          print "$reg\t$li[0]\n"
        }
      }
    }
    if($co>0){
      $co{de}{sw}++;
    }else{
      $co{de}{nsw}++;
    }
  }
  close(IN);

  open(IN,"../Data/Allgenes_$reg\_edgeR.tsv");
  while(<IN>){
    chomp;
    next if($_=~/^Transcript/);
    my @li=split(/\s+/,$_);
    next if(defined $gene{$li[0]});
    my $co=0;
    if(defined $chr{$li[0]}){
      my $chr=$chr{$li[0]};
      if(defined $hash{$chr}){
        foreach my $key (sort {$a <=> $b} keys %{$hash{$chr}}){
          next if($en{$li[0]}<$key-100000 or $st{$li[0]}>$hash{$chr}{$key}+100000);
          $co++;
        }
      }
    }
    if($co>0){
      $co{nde}{sw}++;
    }else{
      $co{nde}{nsw}++;
    }
  }
  close(IN);
  print "$reg\t$co{de}{sw}\t$co{de}{nsw}\t$co{nde}{sw}\t$co{nde}{nsw}\n";

  my $R = Statistics::R->new() ;
  $R->startR ;
  $R->set('a1',$co{de}{sw});
  $R->set('a2',$co{de}{nsw});
  $R->set('b1',$co{nde}{sw});
  $R->set('b2',$co{nde}{nsw});
  $R->send(q`fisher.test(matrix(c(a1,a2,b1,b2), nrow=2, byrow=T))`);
  my $res = $R->read;
  $R->stopR();
  print "$res\n";
}
