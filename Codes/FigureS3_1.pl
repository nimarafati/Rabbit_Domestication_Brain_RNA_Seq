#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my %hash2;
my $last=0;
my $last2;
my $t=0;
my $chr="chr1";
open(IN,"../Data/raw_data/Deltas_vs_OryCun2.liftover.29mCons.mapped.chrUnFixed.sorted.DeltasOver0.8.bed_NOT_ANY_OTHER_EXON.bed");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  if($chr eq "$li[0]" and $li[1] - $last < 50000){
    $last2="$last2\t$li[1]";
    $hash2{$t}=$last2;
  }else{
    $t++;
    $last2=$li[1];
    $hash2{$t}=$last2;
    $hash{$t}=$chr;
  }
  $last=$li[1];
  $chr=$li[0];
}
close(IN);

open(OUT,">../Data/Deltas_1SNPin50kb.txt");
foreach my $key (sort {$a<=>$b} keys %hash2){
  print OUT "$key\t$hash{$key}\t$hash2{$key}\n";
}
