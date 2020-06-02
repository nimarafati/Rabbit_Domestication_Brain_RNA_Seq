#!/usr/bin/perl
#perl DatabaseS1.pl Amy　Hyp Hipp ParTemp
#For summmarizing DEGs with information aboyt up/down-regulation

use warnings;
use strict;

my %sweep;
my %swe_name;
open(IN,"../Data/raw_data/All-Sweeps.bed");
while(<IN>){
  chomp;
  my @li=split(/\s+/,$_);
  $sweep{$li[0]}{$li[1]}=$li[2];
  $swe_name{$li[0]}{$li[1]}=$li[3]
}
close(IN);

my %daf;
open(IN,"../Data/dAF_DEGs_results.tsv");
while(<IN>){
chomp;
	next if($_=~/^TranscriptID/);
	my @li=split(/\t/,$_);
	$daf{$li[0]}=$li[11] if($li[11]>2);
}

my %chr;
my %start;
my %end;
my %gene;
my %geneid;
open(IN,"../Data/raw_data/Ens76_GeneCoordinates.txt");
while(<IN>){
chomp;
	my @li=split(/\t/,$_);
	$gene{$li[4]}=$li[5];
	$geneid{$li[4]}=$li[3];
	$chr{$li[4]}=$li[0];
	$start{$li[4]}=$li[1];
	$end{$li[4]}=$li[2];
}
close(IN);

my %rab;
open(IN,"../Data/raw_data/Ens76_Orthologs.txt");
while(<IN>){
chomp;
	my @li=split(/\t/,$_);
	if(defined $li[2]){
		$rab{$li[1]}=$li[2];
	}
}
close(IN);

my %human;
open(IN,"../Data/raw_data/Ens76_HumanGeneName.txt");
while(<IN>){
chomp;
	my @li=split(/\t/,$_);
	$human{$li[0]}=$li[1];
}
close(IN);

open(OUT,">../Data/UpDown_All.tsv");
print OUT "Brain area\tUp/Down\tGeneID\tTranscriptIDt\tName\tLogFC\tFDR\tAdjacent sweep region\tNumber of SNPs with dAF > 0.9\n";
for (my $n=0; $n<=$#ARGV; $n++){
	my $reg=$ARGV[$n];
	my %hash;
	my %wildup;
	my %domup;
	my ($w,$d)=(0,0);
	open(IN,"../Data/DEGs_$reg\_edgeR.tsv");
	while(<IN>){
	chomp;
		next if($_=~/^TranscriptID/);
		my @li=split(/\t/,$_);
		if($li[7] > 0){
			$wildup{$li[0]}=1;
			$w++;
		}else{
			$domup{$li[0]}=1;
			$d++;
		}
		$hash{$li[0]}="$li[7]\t$li[10]";
	}
	print "wildup:$w\tdomesticup:$d\n";

	foreach my $key(sort {$geneid{$a} cmp $geneid{$b}} keys %domup){
		my $name;
		my $nameg=$geneid{$key};
		if($gene{$key}!~/ENSOCUT/){
			$name=$gene{$key};
		}else{
			if(defined $rab{$key}){
				my $tmp=$rab{$key};
				if(defined $human{$tmp}){
					$name=$human{$tmp};
				}else{
					$name="-";
				}
			}else{
				$name="-";
			}
		}

		my ($sweep, $daf, $t) = ("-", "-", 0); 
		my ($chr, $start, $end) = ($chr{$key}, $start{$key}, $end{$key}); 
		if(defined $sweep{$chr}){
			foreach my $key2 (sort {$a <=> $b} keys %{$sweep{$chr}}){
				if($end > $key2 - 100000 and $start < $sweep{$chr}{$key2} + 100000){
					$sweep = $swe_name{$chr}{$key2} if($t==0);
					$sweep = $sweep.", ".$swe_name{$chr}{$key2} if($t>0);
					$t++;
				}
			}
		}
		$daf = $daf{$key} if(defined $daf{$key});
		print OUT "$reg\tUp\t$nameg\t$key\t$name\t$hash{$key}\t$sweep\t$daf\n";
	}
	foreach my $key(sort {$wildup{$a} cmp $wildup{$b}} keys %wildup){
		my $name;
		my $nameg=$geneid{$key};
		if($gene{$key}!~/ENSOCUT/){
			$name=$gene{$key};
		}else{
			if(defined $rab{$key}){
				my $tmp=$rab{$key};
				if(defined $human{$tmp}){
					$name=$human{$tmp};
				}else{
					$name="-";
				}
			}else{
				$name="-";
			}
		}

		my ($sweep, $daf, $t) = ("-", "-", 0); 
		my ($chr, $start, $end) = ($chr{$key}, $start{$key}, $end{$key}); 
		if(defined $sweep{$chr}){
			foreach my $key2 (sort {$a <=> $b} keys %{$sweep{$chr}}){
				if($end > $key2 - 100000 and $start < $sweep{$chr}{$key2} + 100000){
					$sweep = $swe_name{$chr}{$key2} if($t==0);
					$sweep = $sweep.", ".$swe_name{$chr}{$key2} if($t>0);
					$t++;
				}
			}
		}
		$daf = $daf{$key} if(defined $daf{$key});
		print OUT "$reg\tDown\t$nameg\t$key\t$name\t$hash{$key}\t$sweep\t$daf\n";
	}
}