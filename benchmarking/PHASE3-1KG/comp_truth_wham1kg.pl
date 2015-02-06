#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

comp_truth_wham.pl truth.vcf wham.vcf

Description:

No really, how the hell do you use this thing!

";


my ($help);
my $FP;
my $opt_success = GetOptions('help'    => \$help,
			     'FP'      => \$FP,
			      );

die $usage if $help || ! $opt_success;

my $file  = shift;
my $fileb = shift;
die $usage unless $file && $fileb;
open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

my %truth;
my %bad;

$bad{'MT'}         = 1;
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
$bad{'chrUn_gl000245'} = 1;
$bad{'chrUn_gl000246'} = 1;
$bad{'chrUn_gl000247'} = 1;
$bad{'chrUn_gl000248'} = 1;
$bad{'chr6_cox_hap2'} = 1;
$bad{'chr6_dbb_hap3'} = 1;
$bad{'chr6_mann_hap4'} = 1;
$bad{'chr6_mcf_hap5'} = 1; 
$bad{'chr6_qbl_hap6'} = 1;
$bad{'chr6_ssto_hap7'} = 1;
$bad{'chr11_gl000202_random'} = 1;
$bad{'chr17_gl000203_random'} = 1;
$bad{'chr17_gl000204_random'} = 1;
$bad{'chr17_gl000205_random'} = 1;
$bad{'chr17_gl000206_random'} = 1;
$bad{'chr18_gl000207_random'} = 1;
$bad{'chr19_gl000208_random'} = 1;
$bad{'chr19_gl000209_random'} = 1;
$bad{'chr1_gl000191_random'}  = 1;
$bad{'chr1_gl000192_random'}  = 1; 
$bad{'chr21_gl000210_random'} = 1;
$bad{'chr4_gl000193_random'} = 1;
$bad{'chr4_gl000194_random'} = 1;
$bad{'chr7_gl000195_random'} = 1;
$bad{'chr8_gl000196_random'} = 1;
$bad{'chr8_gl000197_random'} = 1;
$bad{'chr9_gl000198_random'} = 1;
$bad{'chr9_gl000199_random'} = 1;
$bad{'chr9_gl000200_random'} = 1;
$bad{'chr9_gl000201_random'} = 1;

my $n = 0;

VAR: while (<$IN>) {
    chomp;
    next VAR if $_ =~ /^#/;
    my @l = split /\t/, $_;

    $n++;

    for(my $i = $l[3] - 25; $i < $l[3] + 25; $i++){
	$truth{$l[0]}{$i} = $_;
    }
    for(my $i = $l[4] - 25; $i < $l[4] + 25; $i++){
	$truth{$l[0]}{$i} = $_;
    }
}

close $IN;

print STDERR  "n vars in truth set: ", $n, "\n";

open(my $INB, '<', $fileb) or die "Can't open $fileb for reading\n$!\n";

my %HIT;
my %UNI;
my $tp ;
my $fp ; 


VAR: while (<$INB>) {
    chomp;
    next VAR if $_ =~ /^#/;
    
    my @l = split /\t/, $_;

    my %INFO = map {split /;|=/ } $l[7];

#    next VAR if $INFO{MQF} > 0.1;
#    next VAR if $INFO{MQ}  < 50 ; 

    my $endSeqid = -1;
    my $endPos   = -1;
    my $count    =  0;

    if($INFO{BE} ne '.'){
	($endSeqid, $endPos, $count) = split /,/, $INFO{BE};
    }

    next VAR if defined $bad{$l[0]};
    next VAR if defined $bad{$endSeqid};

    $l[0] =~ s/chr//g;
    $endSeqid =~ s/chr//g;

    if(! defined $truth{$l[0]}{$l[1]} && ! defined $truth{$endSeqid}{$endPos}){
	$fp++;

	if($FP){
            print "FP:";
	    print "\t", $_;
            print "\n";
	}
	next VAR;
    }

    if(defined $truth{$l[0]}{$l[1]}){
	next VAR if defined $UNI{ $truth{$l[0]}{$l[1]} } ;
    }
    if(defined $truth{$endSeqid}{$endPos}){
	next VAR if defined $UNI{$truth{$endSeqid}{$endPos}};
    }
    if (defined $truth{$l[0]}{$l[1]}){

	$UNI{$truth{$l[0]}{$l[1]}}++;

	$tp++;

	my @truth = split /\t/, $truth{$l[0]}{$l[1]} ;
        if(! $FP){
	    print $truth[2];
	    print "\t", $_;
	    print "\n";
	}
	next VAR;
    }
    
    if (defined $truth{$endSeqid}{$endPos}){

        $UNI{$truth{$endSeqid}{$endPos}}++;
	$HIT{$endSeqid}{$endPos}++;
	$tp++;

	my @truth = split /\t/, $truth{$endSeqid}{$endPos} ;
        if(! $FP){
	    print $truth[2];
	    print "\t", $_;
	    print "\n";
	    next VAR;
	}
    }

}
close $INB;
print STDERR "tp: $tp fp:$fp\n";

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------


