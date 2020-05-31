#!/usr/bin/perl
#perl Table1_2.pl Amy Hyp Hipp ParTemp
#For summarizing the GREAT results

use strict;
use warnings;

open(OUT,">../Data/GREAT_res_sig.tsv");
print OUT "Region\tFeature\tID\tTerm\tNum_DEGs\tNum_bg\tFold_enrich\tFDR\n";
for my $reg (@ARGV){
  open(IN,"../Data/GREAT_res_$reg.tsv");
  my %hash;
  while(<IN>){
    chomp;
    next if($_=~/^#/);
    my @li=split(/\t/,$_);
    next if($li[5] > 0.05 or $li[6] > 0.05);
    $hash{$li[0]}{$li[6]}{$li[1]}="$li[2]\t$li[13]\t$li[14]\t$li[7]\t$li[6]";
  }
  foreach my $key (sort {$a cmp $b} keys %hash){
    foreach my $key2 (sort {$a <=> $b} keys %{$hash{$key}}){
      foreach my $key3 (sort {$a cmp $b} keys %{$hash{$key}{$key2}}){
        print OUT "$reg\t$key\t$key3\t$hash{$key}{$key2}{$key3}\n";
      }
    }
  }
  close(IN);
}
