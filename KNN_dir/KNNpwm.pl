#!usr/bin/perl
use warnings;
use strict;

use Parallel::Simple qw( prun );
use List::Util qw(shuffle);


sub pwmLength{
	### (pwm dir string, reference of array) -> length of each PWMs (not including the header line) in the dir

	my ($pwmdir, $listRef)=@_;

	my $matLen=2000;	#arbitray large number
	my $fileNum=0;
	foreach (@$listRef){
		open(FILE,"$pwmdir/$_.mat")||die("Failed to open $pwmdir/$_.mat\n!");
		my @file=<FILE>;
		$fileNum++;
		close FILE;

		if ($#file != $matLen) {
			if ($fileNum==1) {
				$matLen=$#file;
			}else{
				die("All PWMs should have same length: $matLen, but $_.mat has length $#file!\n");
			}
			
		}
	}

	return $matLen;	
}



sub pwmMinLength{
	### (pwm dir string, reference of array) -> minimum length of all PWMs (not including the header line) in the dir

	my ($pwmdir, $listRef)=@_;

	my $minLen=2000;	#arbitray large number
	foreach (@$listRef){
		open(FILE,"$pwmdir/$_.mat")||die("Failed to open $pwmdir/$_.mat\n!");
		my @file=<FILE>;
		close FILE;

		if ($#file < $minLen) {
			$minLen=$#file;
		}
	}

	return $minLen;	#return the min length of the PWM
}



sub loadPromat{
	### (protein similarity file name) -> reference of 2D-hash containing PIDs $PIDs{$pi}{$pj}

	my ($file)=@_;
	my %distance;
	my @proteins;
	open(FILE,"$file")||die("Failed to open $file!");
	while(my $line=<FILE>){
		$line=~s/(\s+)$//;
		if($line=~/^\w+/){
			my @ele=split(/\s+/, $line);
			print "ERROR: #proteins is different from #matrix_columns-1" if $#ele!=$#proteins+1;
			foreach(my $i=1;$i<=$#ele;$i++){
				$distance{$ele[0]}{$proteins[$i-1]}=$ele[$i];
			}
		}elsif($line=~/^\s+\w+/){	#header line
			$line=~s/^(\s+)//;
			@proteins=split(/\s+/, $line);
		}
	}
	close FILE;

	# foreach my $x (@proteins) {
	# 	print "$x\n";
	# }
	# print "\n\n------\n";

	my %PIDs;
	for my $pi(@proteins) {	#for my $pi(keys %distance) {	#
		my $sum=0;
		for my $pj(keys %{$distance{$pi}}){
			if ($pi ne $pj) {
				$PIDs{$pi}{$pj}=$distance{$pi}{$pj};
			}
		}
	}

	return \%PIDs;
}


sub predictPWM{
	### (string of name for target PWM that to be predicted, reference of 1D array containing weights, minum length of all PWMs, string of PWM directory) -> reference of 2D array for predicted PWM pwm[position][base]

	my ($targetPWM, $wVectorR, $matLen, $pwmdir)=@_;

	my @predictPWM;
	for (my $i = 0; $i < $matLen; $i++) {
		for (my $j = 0; $j < 4; $j++) {
			$predictPWM[$i][$j]=0;
		}
	}

	for my $pwm(keys %{$wVectorR}) {
		open(FILE,"$pwmdir/$pwm.mat")||die("Failed to open $pwmdir/$pwm.mat\n!");
		my $linum=0;
		while(my $line=<FILE>){
			next unless $line=~/^\d+/;
			
			$line=~s/(\s+)$//;
			my @ele=split(/\s+/, $line);
			for (my $i = 0; $i <= $#ele; $i++) {
				$predictPWM[$linum][$i]+=$ele[$i]*$wVectorR->{$pwm};
			}
			$linum++;
		}
		close FILE;
	}

	# (optional) final normalization to make a distribution
	for (my $i = 0; $i < $matLen; $i++) {
		my $tmpsum=0;
		for (my $j = 0; $j < 4; $j++) {
			$tmpsum+=$predictPWM[$i][$j];
		}
		for (my $j = 0; $j < 4; $j++) {
			$predictPWM[$i][$j]/=$tmpsum;
		}
	}

	return \@predictPWM;
}



sub printPredictPWM{
	### (reference of PWM, output file name with full path) ->none

	my ($pwmR, $pwmFile)=@_;
	open(FILE,">$pwmFile")||die("Failed to open $pwmFile!");
	print FILE "A C G U\n";
	for(@$pwmR){
		for my $j(@$_){
			print FILE "$j ";
		}
		print FILE "\n";
	}
	close FILE;
}



sub getCosine{
	### get cosine similarity between predicted PWM and the true PWM

	my ($pwmFile1, $pwmFile2, $baseNum)=@_;

	my $EXE='./cosineSim-pairwise.py'; 

	my $string=`python $EXE $pwmFile1 $pwmFile2 $baseNum`;
	my @cosVal;
	if ($string=~/maxCos=(\S+), minCos=(\S+), avgCos=(\S+), overallCos=(\S+), median=(\S+)/) {		
		push @cosVal, ($1, $2, $3, $4, $5);
	}else{
		die "ERROR output from python script!\n";
	}

	return \@cosVal;
}



sub wrapPredictPWM{
	my ($targetPWM, $wVectorR, $matLen, $pwmDir, $outputDir)=@_;

	my $pwmR=predictPWM($targetPWM, $wVectorR, $matLen, $pwmDir);

	eval{printPredictPWM($pwmR, "$outputDir/$targetPWM.mat")};
	print "Error: $@\n" if $@;

	# (optional) get cosine similarity for predicted PWM vs true PWM
	# my $baseNum=3;
	# my $csvFile="$outputDir/cosineSim-base$baseNum.csv";
	# my $cosValR=getCosine("$pwmDir/$targetPWM.mat", "$outputDir/$targetPWM.mat", $baseNum);
	# open(FILE,">>$csvFile")||die("Failed to add to $csvFile!");
	# print FILE "$targetPWM,$cosValR->[0],$cosValR->[1],$cosValR->[2],$cosValR->[3],$cosValR->[4]\n";
	# close FILE;
}



sub topKweight{
	my ($k, $targetPWM, $PIDsR)=@_;

	unless(exists $PIDsR->{$targetPWM}){
		die("$targetPWM is not a key in the hash!\n");
	}

	my @sortedPWMs=sort {$PIDsR->{$targetPWM}{$b} <=> $PIDsR->{$targetPWM}{$a}} keys %{$PIDsR->{$targetPWM}};

	my $weightVectorR;
	my $tmpsum=0;
	for (my $i = 0; $i < $k; $i++) {
		$weightVectorR->{$sortedPWMs[$i]}=$PIDsR->{$targetPWM}{$sortedPWMs[$i]};
		$tmpsum+=$PIDsR->{$targetPWM}{$sortedPWMs[$i]};
	}

	for (keys %$weightVectorR){
		$weightVectorR->{$_}/=$tmpsum;
	}

	return $weightVectorR;
}



sub arraySymDiff{
	### (ref of array A, ref of array B) -> ref of array containing elements in either A or B but not both (i.e., symmetric set diff between A and B)

	my ($a, $b) = @_;

	my %count;
	my @diff;

	foreach (@$a, @$b) {
		$count{$_}++;
	}

	foreach (keys %count) {
		if ($count{$_} == 1) {
			push @diff, $_;
		}
	}

	return \@diff;
}



sub foldPredictPWMs{
	### (ref of array of folds start index, ref of array of shuffled protein names, ref of a 2d hash containing subset of PIDs, PWM mat length, pwm dir, output dir, number of nearest neighbors, which fold, how many folds in total) ->none
	
	my ($foldsStartIdxR, $proteinNamesShuffledR, $PIDsR, $matLen, $matdir, $outdir, $k, $f, $n) = @_;

	my $startIdx=$foldsStartIdxR->[$f];
	my $endIdx= $f==$n-1 ? $#{$proteinNamesShuffledR} : $foldsStartIdxR->[$f+1]-1;

	my @testNames = @{$proteinNamesShuffledR}[$startIdx..$endIdx];
	my $trainNamesR = arraySymDiff($proteinNamesShuffledR, \@testNames);

	my %subPIDs;
	foreach my $testName (@testNames) {
		foreach my $trainName (@$trainNamesR) {
			$subPIDs{$testName}{$trainName} = $PIDsR->{$testName}{$trainName};
		}
	}

	# my @run;
	foreach my $testName (@testNames) {
		my $weightVecR = &topKweight($k, $testName, \%subPIDs);
		&wrapPredictPWM($testName, $weightVecR, $matLen, $matdir, $outdir);
		# push @run, [\&wrapPredictPWM, $testName, $weightVecR, $matLen, $matdir, $outdir];
	}
	# prun(@run) or die(Parallel::Simple::errplus());
}



############# Main ##############
# Note: PWMs should be the same length (can be trimmed based on information content using STAMP etc.):
#assume the PWM is in the format like:
#--------headerline----------
#pos1_A	pos1_C	pos1_G	pos1_U
#pos2_A	pos2_C	pos2_G	pos2_U
#...
#each cell contains a frequency value 

die "$0 requires 4 argument.\n $0 <homolog PWM dir with PWMs having the same length> <protein similarity matrix file derived from ClustalW> <optK of KNN> <output dir for predicted PWMs>\n" if $#ARGV!=3;

my $matdir = $ARGV[0]; #homolog PWM folder. The PWMs should have the same length
my $protFile=$ARGV[1];#protein similarity matrix file. Obtained from ClustalW
my $k=$ARGV[2];	#k of knn, choose by user. The optimalK should be choosen from cross-validation
my $outdir=$ARGV[3]; #folder for the predicted PWMs

my $PIDsR=loadPromat($protFile);

my @proteinNames=keys %$PIDsR;

my $matLen=pwmLength($matdir, \@proteinNames);

unless(-e $outdir){
	eval {system("mkdir -p $outdir")};
	print "Error: when make dir $outdir: $@\n" if $@;
}

# n-fold; set $n_fold to number of proteins for leave one out prediction
my $n_fold=$#proteinNames+1; 

# shuffle the protein array order
my @proteinNamesShuffled = shuffle @proteinNames;

# set the start and end index for each fold
my $foldSize=int(($#proteinNames+1)/$n_fold);
my @foldSizes;
foreach (1..$n_fold) {
	push @foldSizes, $foldSize;
}
my $res = $#proteinNames+1 - $foldSize*$n_fold;
for (my $i = 0; $i < $res; $i++) {
	$foldSizes[$i]++;
}
my @foldsStartIdx;
$foldsStartIdx[0]=0;
for (my $i = 1; $i < $n_fold; $i++) {
	$foldsStartIdx[$i] = $foldsStartIdx[$i-1] + $foldSizes[$i-1];
}

my @run;
foreach my $f (0..$n_fold-1) {
	# &foldPredictPWMs(\@foldsStartIdx, \@proteinNamesShuffled, $PIDsR, $matLen, $matdir, $outdir, $k, $f, $n_fold);
	push @run, [\&foldPredictPWMs, \@foldsStartIdx, \@proteinNamesShuffled, $PIDsR, $matLen, $matdir, $outdir, $k, $f, $n_fold];
}
prun(@run) or die(Parallel::Simple::errplus());
