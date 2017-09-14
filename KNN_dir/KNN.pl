#!usr/bin/perl
use warnings;
use strict;

use Parallel::Simple qw( prun );

#die "$0 requires argument.\n nonpredict_ccV2-general.pl <family> <mode cc:1 or cc2:2>\n eg: perl nonpredict_ccV2.pl MF0005 1\n" if $#ARGV!=3;

#my $pwmName=$ARGV[1];
# my $outdir = $ARGV[3]; #'/ubc/cs/research/connections/people/syang11/work/coevo/mat_out';
# unless(-e $outdir){
# 	eval {system("mkdir -p $outdir")};
# 	print "Error: when make dir $outdir: $@\n" if $@;
# }

# my $EXE='/ubc/cs/research/connections/people/syang11/work/coevo/cosineSim-pairwise.py';


# (pwm dir string, reference of array) -> minimum length of all PWMs (not including the header line) in the dir
sub pwmMinLength{
	my ($pwmdir, $listRef)=@_;

	my $minLen=20;
	foreach (@$listRef){
		my $string=`wc -l $pwmdir/$_.mat`;
		if ($string=~/^(\d+)\s+/) {
			if ($1<$minLen) {
				$minLen=$1-1;
			}
		}else{
			print "ERROR: output from wc is illegal: $string\n";
		}
	}

	return $minLen;
}

# (protein similarity file name) -> reference of 2D-hash containing normalized weights $weight{$pi}{$pj}
sub loadPromat{
	my ($file)=@_;
	my %distance;
	my @proteins;
	open(FILE,"$file")||die("Failed to open $file!");
	while(my $line=<FILE>){
		$line=~s/(\s+)$//;
		if($line=~/^(\S+)/){
			my @ele=split(/\s+/, $line);
			print "ERROR: #proteins is different from #matrix_columns-1" if $#ele!=$#proteins+1;
			foreach(my $i=1;$i<=$#ele;$i++){
				$distance{$ele[0]}{$proteins[$i-1]}=$ele[$i];
			}
		}elsif($line=~/\S+/){	#header line
			$line=~s/^(\s+)//;
			@proteins=split(/\s+/, $line);
		}
	}
	close FILE;

	my %weight;
	for my $pi(keys %distance) {
		my $sum=0;
		for my $pj(keys %{$distance{$pi}}){
			if ($pi ne $pj) {
				$weight{$pi}{$pj}=$distance{$pi}{$pj};
				$sum+=$distance{$pi}{$pj};
			}
		}
		for my $pj(keys %{$weight{$pi}}){
			$weight{$pi}{$pj}/=$sum;
		}
	}

	return \%weight;
}

# (string of name for target PWM that to be predicted, reference of 1D array containing weights, minum length of all PWMs, string of PWM directory) -> reference of 2D array for predicted PWM pwm[position][base]
sub predictPWM{
	my ($targetPWM, $wVectorR, $minLen, $pwmdir)=@_;

	my @predictPWM;
	for (my $i = 0; $i < $minLen; $i++) {
		for (my $j = 0; $j < 4; $j++) {
			$predictPWM[$i][$j]=0;
		}
	}

	for my $pwm(keys %{$wVectorR}) {
		open(FILE,"$pwmdir/$pwm.mat")||die("Failed to open $pwmdir/$pwm.mat\n!");
		my $linum=0;
		while(my $line=<FILE>){
			next unless $line=~/^\d+/;
			next unless $linum<$minLen;
			$line=~s/(\s+)$//;
			my @ele=split(/\s+/, $line);
			for (my $i = 0; $i <= $#ele; $i++) {
				$predictPWM[$linum][$i]+=$ele[$i]*$wVectorR->{$pwm};
			}
			$linum++;
		}
		#debug
		# if($targetPWM=~/RNCMPT00239/){
		# 	print "$linum, $minLen, $#predictPWM\n";
		# }
		close FILE;
	}

	return \@predictPWM;
}

# (reference of PWM, output file name with full path) ->none
sub printPredictPWM{
	my ($pwmR, $pwmFile)=@_;
	open(FILE,">$pwmFile")||die("Failed to open $pwmFile!");
	print FILE "AC RNCMPT00007 ID V\$Arnt NA Arnt_KH\n";
	for(@$pwmR){
		for my $j(@$_){
			print FILE "$j ";
		}
		print FILE "\n";
	}
	close FILE;
}

sub wrapPredictPWM{
	my ($targetPWM, $wVectorR, $minLen, $pwmDir, $outputDir, $baseNum)=@_;

	my $pwmR=predictPWM($targetPWM, $wVectorR, $minLen, $pwmDir);

	eval{printPredictPWM($pwmR, "$outputDir/$targetPWM.mat")};
	print "Error: $@\n" if $@;
}

sub topKweight{
	my ($k, $targetPWM, $weightR)=@_;

	unless(exists $weightR->{$targetPWM}){
		die("$targetPWM is not a key in the hash!\n");
	}

	my @sortedPWMs=sort {$weightR->{$targetPWM}{$b} <=> $weightR->{$targetPWM}{$a}} keys %{$weightR->{$targetPWM}};

	my $wVectorR;
	my $tmpsum=0;
	for (my $i = 0; $i < $k; $i++) {
		$wVectorR->{$sortedPWMs[$i]}=$weightR->{$targetPWM}{$sortedPWMs[$i]};
		$tmpsum+=$weightR->{$targetPWM}{$sortedPWMs[$i]};
	}

	for (keys %$wVectorR){
		$wVectorR->{$_}/=$tmpsum;
	}

	return $wVectorR;
}

my $matdir = $ARGV[0]; #homolog PWM folder
my $protFile=$ARGV[1];#protein similarity matrix file. Obtained from ClustalW
my $baseNum=$ARGV[2];
my $k=$ARGV[3];	#k of knn, choose by user. The optimalK should be choosen from cross-validation

my $weightR=loadPromat($protFile);
my @tmp=keys %$weightR;
my $minMatLen=pwmMinLength($matdir, \@tmp);

my $testDir="./test";
unless(-e $testDir){
	eval {system("mkdir -p $testDir")};
	print "Error: when make dir $testDir: $@\n" if $@;
}

my @run;
foreach my $pi(@tmp){
	my $wVecR=topKweight($k, $pi, $weightR);
	push @run, [\&wrapPredictPWM, $pi, $wVecR, $minMatLen, $matdir, $testDir, $baseNum];
}
prun(@run) or die(Parallel::Simple::errplus());


