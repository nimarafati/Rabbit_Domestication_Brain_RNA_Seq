#!/usr/bin/perl
#perl TableS2_2.pl Amy Hyp Hipp ParTemp
#For summarizing the DAVID results

use strict;
use warnings;

open(OUT,">../Data/DAVID_res_sig.tsv");
print OUT "Region\tFeature\tTerm\tNum_DEGs\tNum_bg\tFold_enrich\tFDR\n";
for my $reg (@ARGV){
  open(IN,"../Data/DAVID_res_$reg.txt");
  my %hash;
  while(<IN>){
    chomp;
    next if($_=~/^Category/);
    my @li=split(/\t/,$_);
    next if($li[12] > 0.05);
    $hash{$li[0]}{$li[12]}{$li[1]}="$li[1]\t$li[2]\t$li[7]\t$li[9]\t$li[12]";
  }
  foreach my $key (sort {$a cmp $b} keys %hash){
    foreach my $key2 (sort {$a <=> $b} keys %{$hash{$key}}){
      foreach my $key3 (sort {$a cmp $b} keys %{$hash{$key}{$key2}}){
        print OUT "$reg\t$key\t$hash{$key}{$key2}{$key3}\n";
      }
    }
  }
  close(IN);
}