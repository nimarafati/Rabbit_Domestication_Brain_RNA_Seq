#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my %end;
my $t=0;
open(IN,"../Data/raw_data/All_dAF_0.5over");
while(<IN>){
  chomp;
  $t++;
  my @li=split(/\s+/,$_);
  if($t>1){
    $hash{$li[0]}{$li[1]}=$li[2];
    $end{$li[0]}=$li[1];
  }
}
close(IN);

my %daf;
my %num;
foreach my $chr (sort {$a cmp $b} keys %hash){
  OUT:for (my $n=0;$n<$end{$chr};$n+=5000){
    my $daf;
    my $ndaf;
    foreach my $pos (sort {$a<=>$b} keys %{$hash{$chr}}){
      next if($pos < $n);
      next OUT if($pos > $n + 10000);
      $daf+=$hash{$chr}{$pos};
      $ndaf++;
      $daf{$chr}{$n}=$daf/$ndaf;
      $num{$chr}{$n}=$ndaf;
    }
  }
}

open(OUT,">../Data/dAF_window_10kb.txt");
print OUT "CHR\tSTART\tEND\tNUM_SNP\tAVE_dAF\n";
foreach my $chr (sort {$a cmp $b} keys %daf){
  foreach my $pos (sort {$a<=>$b} keys %{$daf{$chr}}){
    my $tmp=$pos+10000;
    print OUT "$chr\t$pos\t$tmp\t$num{$chr}{$pos}\t$daf{$chr}{$pos}\n" if(defined $daf{$chr}{$pos});
    print OUT "$chr\t$pos\t$tmp\t$num{$chr}{$pos}\tNA\n" unless(defined $daf{$chr}{$pos});
  }
}
