#!/usr/bin/perl
#perl DatabaseS3.pl Amy Hyp Hipp ParTemp
#For adding gene names for DAVID outputs

use warnings;
use strict;

my %human;
open(IN,"../Data/raw_data/Ens76_HumanGeneName.txt");
while(<IN>){
chomp;
	next if($_=~/^Ensembl/);
	my @li=split(/\t/,$_);
	$human{$li[0]}=$li[1];
}
close(IN);

for (my $a=0; $a<=$#ARGV; $a++){
	my $reg=$ARGV[$a];
	open(OUT,">../Data/DAVID_res_$reg\_sig.tsv");
	open(IN,"../Data/DAVID_res_$reg.txt");
    while(<IN>){
    chomp;
		my @li=split(/\t/,$_);
		next if($_=~/^Category/ or $li[12] > 0.05);
		my @gene=split(/, /,$li[5]);
		for (my $a=0;$a<=$#gene;$a++){
			my $id = $gene[$a];
			if(defined $human{$id}){
				$gene[$a]=$human{$id};
			}
		}
		$li[5] = join(", ",@gene);
		my $tmp = join("\t",@li);
		print OUT "$tmp\n";
    }
    close(IN);
}