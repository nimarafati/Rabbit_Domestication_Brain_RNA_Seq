#!/usr/bin/perl
#perl Figure1.pl ../Data/DEGs_4regions_edgeR.tsv
#This script have to be run after DatabseS2.pl

use warnings;
use strict;

my %gene;
my %geneid;
open(IN,"../Data/raw_data/Ens76_GeneCoordinates.txt");
while(<IN>){
chomp;
	my @li=split(/\t/,$_);
	$gene{$li[4]}=$li[5];
	$geneid{$li[4]}=$li[3];
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

my $fileout=substr($ARGV[0], 0, -4);
open(OUT,">$fileout\_name.tsv");
print OUT "TranscripID\tGeneID\tName\n";
open(IN,$ARGV[0]);
while(<IN>){
chomp;
	next if($_=~/^TranscriptID/);
	my @li=split(/\t/,$_);
	for (my $a=0; $a<=$#li; $a++){
		next if($li[$a]!~/^ENSOCUT/);
		my $name;
		my $nameg=$geneid{$li[$a]};
		if($gene{$li[$a]}!~/ENSOCUT/){
			$name=$gene{$li[$a]};
		}else{
			if(defined $rab{$li[$a]}){
				my $tmp=$rab{$li[$a]};
				if(defined $human{$tmp}){
					$name=$human{$tmp};
				}else{
					$name="-";
				}
			}else{
				$name="-";
			}
		}
		print OUT "$li[$a]\t$geneid{$li[$a]}\t$name\n";
		last ;
	}
}