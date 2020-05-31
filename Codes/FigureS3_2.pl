#!/usr/bin/perl
use strict;
use warnings;

my $tmp2=100000;
my $t=0;
open(OUT,">../Data/Deltas_1SNPin50kb_2.txt");
open(IN,"../Data/Deltas_1SNPin50kb.txt");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  if($li[-1]-$li[2]<50000){
    unless (defined $li[4]){
      print OUT "$li[1]\t$li[2]\n";
    }else{
      my $tmp2;my $tmp3;my $t=0;
      my $tmp=($li[-1]+$li[2])/2;
      for (my $a=2; $a<=$#li; $a++){
        if($t==0){
          $tmp2=abs($tmp-$li[$a]);
          $tmp3=$li[$a];
        }elsif(abs($tmp-$li[$a])<$tmp2){
          $tmp3=$li[$a];
          $tmp2=abs($tmp-$li[$a]);
          $t++;
        }
      }
      print OUT "$li[1]\t$tmp3\n";
    }
  }else{
    if(defined $li[4]){
      if($li[4] eq $li[-1]){
        print OUT "$li[1]\t$li[3]\n";
      }else{
        print OUT "$li[1]\t$li[2]\n$li[1]\t$li[-1]\n";
      }
    }else{
      print OUT "$li[1]\t$li[2]\n";
    }
  }
}
close(IN);
