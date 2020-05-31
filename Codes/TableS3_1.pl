#!/usr/bin/perl
#perl TableS3_1.pl Amy Hyp Hipp ParTemp
use strict;
use warnings;

my %hash;
my $t=0;
open(IN,"../Data/raw_data/All_dAF_0.5over");
while(<IN>){
  chomp;
  $t++;
  next if($t<2);
  my @li=split(/\t/,$_);
  my $chr;
  if($li[0]=~/chr(\S+)/){$chr=$1};
  $hash{$chr}{$li[1]}=$li[2];
  $t++;
}
close(IN);

my %degs;
for (my $a=0; $a<=$#ARGV; $a++){
  my $reg=$ARGV[$a];
  open(IN, "../Data/DEGs_$reg\_edgeR.tsv");
  while(<IN>){
  chomp;
    next if($_=~/^TranscriptID/);
    my @li=split(/\t/,$_);
    $degs{$reg}{$li[0]}=1;
  }
}

my %deg;
my %deg2;
$t=0;
open(IN,"../Data/Allgenes_edgeR.tsv");
while(<IN>){
  chomp;
  $t++;
  next if($t<2);
  my @li=split(/\s+/,$_);
  for (my $a=0; $a<=$#ARGV; $a++){
    my $reg=$ARGV[$a];
    my $a1 = $a * 4 + 1;
    my $a2 = ($a+1) * 4;
    if($li[$a1] eq "NA" or $li[$a2] eq "NA"){
      $deg{$reg}{$li[0]}="NA";
    }elsif(defined $degs{$reg}{$li[0]}){
      $deg{$reg}{$li[0]}=1;
    }else{
      $deg{$reg}{$li[0]}=0;
    }
  }
  if($deg{Amy}{$li[0]} eq "1" and $deg{Hyp}{$li[0]} eq "1" and $deg{Hipp}{$li[0]} eq "1" and $deg{ParTemp}{$li[0]} eq "1"){
    $deg2{$li[0]}=1;
  }else{
    $deg2{$li[0]}=0;
  }
}

my %hash2;
$t=0;
open(IN,"../Data/raw_data/Ens76_Rabbit_genename_ID_chr_start_end_strand.txt");
while(<IN>){
  chomp;
  $t++;
  next if($t<2);
  my @li=split(/\s+/,$_);
  if(defined $li[6]){
    $hash2{$li[2]}{name}=$li[0] if($li[0]=~/(\S+)/);
    $hash2{$li[2]}{name}="NA" unless($li[0]=~/(\S+)/);
    $hash2{$li[2]}{chr}=$li[3];
    $hash2{$li[2]}{start}=$li[4];
    $hash2{$li[2]}{end}=$li[5];
    $hash2{$li[2]}{strand}=$li[6];
  }else{
    $hash2{$li[1]}{name}="NA";
    $hash2{$li[1]}{chr}=$li[2];
    $hash2{$li[1]}{start}=$li[3];
    $hash2{$li[1]}{end}=$li[4];
    $hash2{$li[1]}{strand}=$li[5];
  }
}

open(OUT,">../Data/dAF_DEGs_results.tsv");
print OUT "TranscriptID\tGenename\tChr\tStart\tEnd\tStrand\tAve_gene_dAF\tNum_gene_dAF\tAve_up100kb_dAF\tNum_up100kb_dAF\tupdown50kb_dAF0.9\tupdown100kb_dAF0.9\tDEGs_Amy\tDEGs_Hyp\tDEGs_Hipp\tDEGs_ParTemp\tDEGs_All\n";
foreach my $key (keys %deg2){
  if(defined $hash2{$key}{chr}){
    my $chr=$hash2{$key}{chr};
    my ($sum_gdaf,$av_gdaf,$sg,$n_gdaf,$sum_udaf,$av_udaf,$su,$n_udaf,$gen_len,$updown,$updown_10)=(0,0,0,0,0,0,0,0,0,0,0);
    foreach my $snp (sort{$a <=> $b} keys %{$hash{$chr}}){
      $gen_len=$hash2{$key}{end}-$hash2{$key}{start};
      if($snp >= $hash2{$key}{start} and $snp <= $hash2{$key}{end}){
        $sum_gdaf+=$hash{$chr}{$snp};
        $sg++;
      }
      if($snp > $hash2{$key}{start}-50000 and $snp < $hash2{$key}{end}+50000 and $hash{$chr}{$snp}>0.9){
        $updown+=1;
      }
      if($snp > $hash2{$key}{start}-100000 and $snp < $hash2{$key}{end}+100000 and $hash{$chr}{$snp}>0.9){
        $updown_10+=1;
      }
      if($hash2{$key}{strand} == 1){
        if($snp > $hash2{$key}{start}-100000 and $snp <= $hash2{$key}{start}){
          $sum_udaf+=$hash{$chr}{$snp};
          $su++;
        }
      }else{
        if($snp >= $hash2{$key}{end} and $snp < $hash2{$key}{end}+100000){
          $sum_udaf+=$hash{$chr}{$snp};
          $su++;
        }
      }
    }
    if($sg>0){
      $av_gdaf=$sum_gdaf/$sg;
      $n_gdaf=$sg/$gen_len*10000;
    }else{
      $av_gdaf="NA";
      $n_gdaf=0;
    }
    if($su>0){
      $av_udaf=$sum_udaf/$su;
      $n_udaf=$su/10;
    }else{
      $av_udaf="NA";
      $n_udaf=0;
    }
    print OUT "$key\t$hash2{$key}{name}\t$hash2{$key}{chr}\t$hash2{$key}{start}\t$hash2{$key}{end}\t$hash2{$key}{strand}\t$av_gdaf\t$n_gdaf\t$av_udaf\t$n_udaf\t$updown\t$updown_10\t$deg{Amy}{$key}\t$deg{Hyp}{$key}\t$deg{Hipp}{$key}\t$deg{ParTemp}{$key}\t$deg2{$key}\n";
  }else{
    print "$key\n";
  }
}
