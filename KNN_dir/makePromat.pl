#!usr/bin/perl
use warnings;
use strict;

#read from the clustalW2 outputs, generate the PID score matrix accordingly
#served as a preparation part of PWM prediction

die "$0 requires argument.\n $0 <input folder> <fasta file name> <output folder> <matrix file name>\n" if $#ARGV!=3;

my $inputdir=$ARGV[0];
my $fafile=$ARGV[1];
my $outputdir=$ARGV[2];
my $mafile=$ARGV[3];

my @arr;

#get familiy members of MF
my @MF;
open(LIST,"$inputdir/$fafile")||die("Failed to open $inputdir/$fafile file!");
while(<LIST>){
	chomp;
	if(/^\>(\w+)/){
		push @MF, $1;
	}
}

unless(-e $outputdir){
	eval {system("mkdir -p $outputdir")};
	print "Error: when make dir $outputdir: $@\n" if $@;
}

#call clustalw2
system("clustalw2 -INFILE=$inputdir/$fafile -ALIGN -TYPE=PROTEIN -OUTFILE=$outputdir/tmp.aln >>$outputdir/alignment.txt");

open(FILE, "$outputdir/alignment.txt")||die("Failed to open $outputdir/alignment.txt!");
while(<FILE>){
	chomp;
	if(/^Sequences \((\d+):(\d+)\) Aligned. Score:  (\d+)/){
		$arr[$1-1][$2-1]=$3;
		$arr[$2-1][$1-1]=$3;
	}
}
close(FILE);

for (0..$#arr){
	$arr[$_][$_]=100;
}

my $dnd=(split(/\./, $fafile))[0];
unlink "$outputdir/tmp.aln","$outputdir/alignment.txt","$inputdir/$dnd.dnd";

my ($i,$j);
open(FILE,">$outputdir/$mafile")||die("Failed to open $mafile!");
for (@MF) {
	print FILE "\t$_";
}
print FILE "\n";
for $i(0..$#arr){
	print FILE "$MF[$i]";
	for $j(0..$#{$arr[$i]}){
		print FILE "\t$arr[$i][$j]";
	}
	print FILE "\n";
}
close(FILE);



