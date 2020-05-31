#!/usr/bin/perl
#perl TableS3_2.pl Amy Hyp Hipp ParTemp
use strict;
use warnings;
use Statistics::R;

for (my $a=0; $a<=$#ARGV; $a++){
  my %hash;
  $hash{de}{sw} = 0; $hash{de}{nsw} = 0;
  $hash{nde}{sw} = 0; $hash{nde}{nsw} = 0;
  my $reg=$ARGV[$a];
  my $n=$a+12;
  open(IN,"../Data/dAF_DEGs_results.tsv");
  while(<IN>){
    chomp;
    my @li=split(/\t/,$_);
    next if($li[0]=~/TranscriptID/);
    if($li[$n] ne "NA" and $li[$n] == 1){
      if($li[11] > 2){
        $hash{de}{sw}++ if(defined $hash{de}{sw});
        $hash{de}{sw}=1 unless(defined $hash{de}{sw});
        print "$reg\t$li[0]\n";
      }else{
        $hash{de}{nsw}++ if(defined $hash{de}{nsw});
        $hash{de}{nsw}=1 unless(defined $hash{de}{nsw});
      }
    }elsif($li[$n] ne "NA" and $li[$n] == 0){
      if($li[11] > 2){
        $hash{nde}{sw}++ if(defined $hash{nde}{sw});
        $hash{nde}{sw}=1 unless(defined $hash{nde}{sw});
      }else{
        $hash{nde}{nsw}++ if(defined $hash{nde}{nsw});
        $hash{nde}{nsw}=1 unless(defined $hash{nde}{nsw});
      }
    }
  }
  close(IN);

  print "$reg\t$hash{de}{sw}\t$hash{de}{nsw}\t$hash{nde}{sw}\t$hash{nde}{nsw}\n";
  my $R = Statistics::R->new() ;
  $R->startR ;
  $R->set('a1',$hash{de}{sw});
  $R->set('a2',$hash{de}{nsw});
  $R->set('b1',$hash{nde}{sw});
  $R->set('b2',$hash{nde}{nsw});
  $R->send(q`fisher.test(matrix(c(a1,a2,b1,b2), nrow=2, byrow=T))`);
  my $res = $R->read;
  $R->stopR();
  print "$res\n";
}