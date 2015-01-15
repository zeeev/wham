#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

truth_lumpy.pl -t INR|INV|DUP|DEL -w 25 truth.bed lumpy.bedbe 

Description:

-t <FLAG>   --required-- type of structural variant in truth.bed file
-w <INT>    --optional-- the number of bases on either side of the true event to count as acceptable [25]

Description:

This script takes \"truth\" bed file adds slop around the breakpoints
(pos1 & pos2) and then intersects LUMPY calls, independent of SV type. True
positives are determined by MAX.

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
my $TYPE;
my $TP;
my $WINDOW = 25;
my $DEPTH  = 10;
my $opt_success = GetOptions('help'     => \$help,
			     "TP"       => \$TP,
			     "window=s" => \$WINDOW,
			     "depth=s"  => \$DEPTH,
			     "TYPE=s"   => \$TYPE );

die $usage if $help || ! $opt_success || ! defined $TYPE ;

my $file  = shift;
my $fileb = shift;
die $usage unless $file;
open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

#1       7851301 7851302 1       7851202 7851203 1::DUP0444::1   255     +       +

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

my %bad;

$bad{'hs37d5'}     = 1;
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

my %HIT;

SV: while (<$INB>) {

#1       7017630 7017662 1       8017622 8017664 1       4.15492e-27     +       -       TYPE:DELETION   IDS:1,9;2,3     STRANDS:+-,12   MAX:1:7017643;1:8017664 95:1:7017641-7017644;1:8017664-8017664
    chomp;
    $_ =~ s/\.//g;
    next if $_ =~ /^#/;
    my @l = split /\t/, $_;
    
    next SV if defined $bad{$l[0]} ;
    next SV if defined $bad{$l[5]} ;

#MAX:1:7017643;1:8017664
    my @max = split /:|;/, $l[13];

    shift @max;
    
    if(defined $truth{$max[0]}{$max[1]} || defined $truth{$max[2]}{$max[3]}){
	
	my $seqid = 0;
	my $key   = 1;
	
	if(defined $truth{$max[0]}{$max[1]}){
	    $key   = 1;
	    $seqid = 0;
	}
	if(defined $truth{$max[2]}{$max[3]}){
	    $seqid = 2;
	    $key   = 3;
	}
	
	if(defined $HIT{$truth{$max[$seqid]}{$max[$key]}}){
	    next SV;
	}
	
	if($TP){
#	    print  $truth{$max[$seqid]}{$max[$key]}, "\n";
	    print "$_\n";
	}
	
	$tp ++;
	
	$HIT{$truth{$max[$seqid]}{$max[$key]}} = 1;
	
    }
    else{
#	print STDERR $_ , "\n";
	$fp++;
    }
}



close $INB;

my $FDR = $fp / ($tp + $fp);

print STDERR "LUMPY\t$TYPE\t$WINDOW\t$fp\t$tp\t$FDR\n";


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

