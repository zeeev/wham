#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

truth.pl -t INR|INV|DUP|DEL -w 25 truth.bed wham.vcf 

Options:

-t <FLAG>   --required-- type of structural variant in truth.bed file
-C <INT>    --optional-- number of softClipped positions in pileup [10000]
-E <INT>    --optional-- number of reads supporting endpoint       [0]
-N <INT>    --optional-- number of reads in consensus              [0]
-TP <FLAG>  --optional-- print true positives
-FP <FLAG>  --optional-- print false positives 
-w <INT>    --optional-- the number of bases on either side of the true event 
                         to count as acceptable [25]

Description:

This script takes \"truth\" bed file adds slop around the breakpoints 
(pos1 & pos2) and then intersects WHAM calls, independent of SV type. 

The truth files used in our experiments can be found with: Chaisson et al 2014 
(http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/)

Output:

If the true positive or false positive flag is toggled the hits will be printed
to SDOUT.

By default the summary statisitcs are printed to SDERR.

program type truth-window FP TP FDR
WHAM INDEL 25 18546 5239 0.779735127180996


";

my ($help);

my $WINDOW = 25;
my $TYPE;
my $GENOTYPE; 
my $FP;
my $TP;
my $EP = 0;
my $PU = 0;
my $SU = 0;
my $NC = 0;
my $CU = 10000;
my $RD = 10000;

my $opt_success = GetOptions('help'     => \$help,
			     "window=s" => \$WINDOW,
			     'GENOTYPE' => \$GENOTYPE,
			     "TYPE=s"   => \$TYPE,
			     "PU=s"     => \$PU,
			     "SU=s"     => \$SU,
			     "CU=s"     => \$CU,
			     "EP=s"     => \$EP,
			     "NC=s"     => \$NC,
			     "RD=s"     => \$RD,
			     'TP'       => \$TP,
			     'FP'       => \$FP
    );

die $usage if $help || ! $opt_success || ! defined $TYPE;

my $file  = shift;
my $fileb = shift;
die $usage unless $file;
open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

#chr1    726766  726851  insertion       85      (GAATG)n:FULL
my %truth;
my $ntruth = 0;

while (<$IN>) {

    $ntruth++;    
    chomp;

    $_ =~ s/chr//g;
    my @l = split /\s+/, $_;
    
    for(my $i = $l[1] - $WINDOW; $i < $l[1] + $WINDOW; $i++){	
	$truth{$l[0]}{$i} = $_;
    }
    for(my $i = $l[2] - $WINDOW; $i < $l[2] + $WINDOW; $i++){
	$truth{$l[0]}{$i} = $_;
    }
}

close $IN;

open(my $INB, '<', $fileb) or die "Can't open $fileb for reading\n$!\n";

my $tp    = 0;
my $fp    = 0;

my %HIT;
my $UNI;
my %bad;

$bad{'hs37d5'}     = 1;
$bad{'MT'}         = 1;
$bad{'GL0002191'}  = 1;
$bad{'GL0002241'}  = 1;
$bad{'GL0002211'}  = 1;
$bad{'GL0001921'}  = 1;
$bad{'GL0001931'}  = 1;
$bad{'GL0001951'}  = 1;
$bad{'GL0001991'}  = 1;
$bad{'GL0002041'}  = 1;
$bad{'GL0002051'}  = 1;
$bad{'GL0002071'}  = 1;
$bad{'GL0002081'}  = 1;
$bad{'GL0002121'}  = 1;
$bad{'GL0002141'}  = 1;
$bad{'GL0002161'}  = 1;
$bad{'GL0002171'}  = 1;
$bad{'GL0002201'}  = 1;
$bad{'GL0002251'}  = 1;
$bad{'GL0002261'}  = 1;
$bad{'GL0002281'}  = 1;
$bad{'GL0002321'}  = 1;
$bad{'GL0002341'}  = 1;
$bad{'GL0002351'}  = 1;
$bad{'GL0002371'}  = 1;
$bad{'GL0002391'}  = 1;
$bad{'GL0002401'}  = 1;
$bad{'GL000191.1'} = 1;
$bad{'GL000192.1'} = 1;
$bad{'GL000193.1'} = 1;
$bad{'GL000194.1'} = 1;
$bad{'GL000195.1'} = 1;
$bad{'GL000198.1'} = 1;
$bad{'GL000199.1'} = 1;
$bad{'GL000202.1'} = 1;
$bad{'GL000203.1'} = 1;
$bad{'GL000204.1'} = 1;
$bad{'GL000205.1'} = 1;
$bad{'GL000206.1'} = 1;
$bad{'GL000207.1'} = 1;
$bad{'GL000208.1'} = 1;
$bad{'GL000210.1'} = 1;
$bad{'GL000211.1'} = 1;
$bad{'GL000212.1'} = 1;
$bad{'GL000213.1'} = 1;
$bad{'GL000214.1'} = 1;
$bad{'GL000215.1'} = 1;
$bad{'GL000216.1'} = 1;
$bad{'GL000217.1'} = 1;
$bad{'GL000218.1'} = 1;
$bad{'GL000219.1'} = 1;
$bad{'GL000220.1'} = 1;
$bad{'GL000221.1'} = 1;
$bad{'GL000222.1'} = 1;
$bad{'GL000224.1'} = 1;
$bad{'GL000225.1'} = 1;
$bad{'GL000226.1'} = 1;
$bad{'GL000227.1'} = 1;
$bad{'GL000228.1'} = 1;
$bad{'GL000229.1'} = 1;
$bad{'GL000230.1'} = 1;
$bad{'GL000231.1'} = 1;
$bad{'GL000232.1'} = 1;
$bad{'GL000234.1'} = 1;
$bad{'GL000235.1'} = 1;
$bad{'GL000236.1'} = 1;
$bad{'GL000237.1'} = 1;
$bad{'GL000238.1'} = 1;
$bad{'GL000239.1'} = 1;
$bad{'GL000240.1'} = 1;
$bad{'GL000241.1'} = 1;
$bad{'GL000242.1'} = 1;
$bad{'GL000243.1'} = 1;
$bad{'GL000244.1'} = 1;
$bad{'GL000247.1'} = 1;
$bad{'GL000248.1'} = 1;

SV: while (<$INB>) {
    chomp;
    next SV if $_ =~ /^#/;

    my @l = split /\t/, $_;
   
    next SV if defined $bad{$l[0]};

    my %info = map {split /;|=/} $l[7];

#    next SV if $info{"MQF"} > 0.05;

    my @alt;
    my $endSeqid = "NAN";
    my $endPos   = "NAN";
    my $endCount = 0;

    my $number = () = $l[4] =~ /N/gi;

    if($info{"BE"} ne '.'){
	@alt = split /,/, $info{"BE"};
	$endSeqid = $alt[0];
	$endPos   = $alt[1];
	$endCount = $alt[2];
    } 
    else{	
	next SV if $info{"SI"} < 0.5 ;
	next SV if $info{"PU"} < 6;
	next SV if $info{"NC"} < 6 && $info{"SU"} < 1;
    }
   
#    if($info{"SVLEN"} ne '.'){
#	next SV if $info{"SVLEN"} < 50;
#    }

    next SV if ($number / $info{"NC"}) > 1.5 ;
    next SV if $number > 18;
    next SV if $info{"CU"} > 1.40 * $info{"RD"};

    next SV if $endCount   < $EP;
    next SV if $info{"NC"} < $NC;
    next SV if $info{"CU"} > $CU;
    next SV if $info{"RD"} > $RD;
    
    next SV if defined $bad{$endSeqid};
    
    if(defined $truth{$l[0]}{$l[1]} || defined $truth{$endSeqid}{$endPos}){

	my $skey = $l[0];
	my $pkey = $l[1];

	if(defined $truth{$endSeqid}{$endPos}){
	    $skey = $endSeqid;
	    $pkey = $endPos;
	}
	if(defined $HIT{$truth{$skey}{$pkey}}){
	    next SV;
	}
	if($TP){
	    print "TP:\t$_\n";
#	    print  $truth{$skey}{$pkey} , "\n";
	}
	$tp ++;

	my $flagg = 0;
	if($l[9] =~ /1\/1/){
            $flagg = 1;
	}
	if($GENOTYPE && $TP){
	    print "\t$flagg\n";
	}
	else{
	    if($TP){
		print "\n";
	    }
	}
	if($skey ne -1){
	    $HIT{$truth{$skey}{$pkey}} = 1;
	}
    }
    else{
	if($FP){
	    print "FP:\t", $_, "\n";
	}
	$fp++;
    }
}



close $INB;

my $FDR = $fp / ($tp + $fp);

print STDERR "WHAM\t$TYPE\t$WINDOW\t$fp\t$tp\t$FDR\n";


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

