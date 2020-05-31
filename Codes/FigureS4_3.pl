#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my %hash3;
my @regs=("Amy","Hyp","Hipp","ParTemp");
foreach my $reg (@regs){
  for (my $a=50;$a<=95;$a=$a+5){
    my $a2=$a/100;
    my $b=$a2+0.05;
    my $temp="$a2\-$b";
    $hash{$temp}{$reg}=0;
    $hash3{$reg}{"Yes"}=0;
    $hash3{$reg}{"No"}=0;
  }
}

my %hash2;
my %count;
#my $count;
open(IN,"../Data/All_dAF_TSS100kb_DEGs.txt");
while(<IN>){
  chomp;
  my @li=split(/\t/,$_);
  next if($li[4]=~/-/ or $li[0] eq "Chr");
  my $daf=$li[3];
  if($daf eq "1-1.05"){
    $daf="0.95-1";
  }
  $count{$daf}++ if(defined $count{$daf});
  $count{$daf}=1 unless(defined $count{$daf});
#  $count++;
  for($a=5;$a<=8;$a++){
    my $t=$a-5;
    if($li[$a] eq "Yes"){
      $hash{$daf}{$regs[$t]}++;
      $hash3{$regs[$t]}{"Yes"}++;
    }else{
      $hash2{$daf}{$regs[$t]}++;
      $hash3{$regs[$t]}{"No"}++;
    }
  }
}

open(OUT,">../Data/All_dAF_TSS100kb_DEGs_res.txt");
print OUT "Bin\tTotal\tRegion\tDEGs\tNum\tTotal_Reg_deg\n";
foreach my $key (sort {$a cmp $b} keys %hash){
  foreach my $key2 (keys %{$hash{$key}}){
    my $tmp1 = $hash3{$key2}{Yes};
    my $tmp2 = $hash3{$key2}{No};
    print OUT "$key\t$count{$key}\t$key2\tDEGs\t$hash{$key}{$key2}\t$tmp1\n";
    print OUT "$key\t$count{$key}\t$key2\tnonDEGs\t$hash2{$key}{$key2}\t$tmp2\n";
  }
}
