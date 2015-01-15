#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

truth.pl -t INR|INV|DUP|DEL -w 25 truth.events wham.vcf 

Description:

-t <STRING> --required-- type defined by SVsim.
-w <INT>    --optional-- the number of bases on either side of the true event to count as acceptable [25]
-d <INT>    -- depth  -- the depth of bam file used for simulation [10]

";

my ($help);
my $TYPE;
my $WINDOW = 25;
my $DEPTH  = 10;
my $GENOTYPE; 
my $opt_success = GetOptions('help'     => \$help,
			     "window=s" => \$WINDOW,
			     "depth=s"  => \$DEPTH,
			     'GENOTYPE' => \$GENOTYPE,
			     "type=s"   => \$TYPE );

die $usage if $help || ! $opt_success || ! $TYPE;

my $file  = shift;
my $fileb = shift;
die $usage unless $file;
open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

#1       7851301 7851302 1       7851202 7851203 1::DUP0444::1   255     +       +

my %truth;
my $ntruth = 0;

while (<$IN>) {

    $ntruth++;

    $_ =~ s/\.//g;

    chomp;

    my @l = split /\s+/, $_;

    if($TYPE eq "DEL"){
	
	#DEL0022 DEL     50000   8       107462299       8       107462299       +       8       107512298       +
    
	my %load ;
	$load{"ID"}    = $l[0];
	$load{"TYPE"}  = $l[1];
	$load{"LENGH"} = $l[2];
	$load{"CHR"}   = $l[3];
	$load{"START"} = $l[4];
	$load{"END"}   = $l[9];
	$load{"LINE"}   = $_;

	for(my $i = $l[4] - $WINDOW; $i < $l[4] + $WINDOW; $i++){
	    $truth{$l[3]}{$i} = \%load;
	}

	for(my $i = $l[6] - $WINDOW; $i < $l[6] + $WINDOW; $i++){
            $truth{$l[3]}{$i} = \%load;
	}
	
	for(my $i = $l[9] - $WINDOW; $i < $l[9] + $WINDOW; $i++){
	    $truth{$l[3]}{$i} = \%load;
	}
    }
    if($TYPE eq "INR"){
#    INR0001 INR     500000  15      43942052        Y       16875461        +       Y       17375460        +
#INR0502 INR     50      2       127792228       1       6589207 +       1       6589256 +
	my %load;

	$load{"ID"}     = $l[0];
        $load{"TYPE"}   = $l[1];
        $load{"LENGH"}  = $l[2];
        $load{"CHRL"}   = $l[3];
        $load{"CHRR"}   = $l[5];
	$load{"START"}  = $l[4];
	$load{"STARTR"} = $l[6];
        $load{"ENDR"}   = $l[9];
	$load{"LINE"}   = $_;

	#right at the spot

	for(my $i = $l[4] - $WINDOW; $i < $l[4] + $WINDOW; $i++){
	    $truth{$l[3]}{$i} = \%load;
        }

	# other chr left

	for(my $i = $l[6] - $WINDOW; $i < $l[6] + $WINDOW ; $i++){
            $truth{$l[5]}{$i} = \%load;
        }

	#other chr upstream
	for(my $i = $l[9] - $WINDOW; $i < $l[9] + $WINDOW ; $i++){
            $truth{$l[8]}{$i} = \%load;
        }

    }
    if($TYPE eq "INV" || $TYPE eq "DUP"){
    #INV0001 INV     500000  15      43942052        15      43942052        +       15      44442051        +
    #DUP0001 DUP     500000  15      43942052        15      43942052        +       15      44442051        +
	my %load;

        $load{"ID"}   =  $l[0];
        $load{"TYPE"}  = $l[1];
	$load{"LENGH"} = $l[2];
	$load{"CHR"}   = $l[3];
        $load{"START"}  = $l[4];
	$load{"STARTR"} = $l[6];
        $load{"ENDR"}   = $l[9];
	$load{"LINE"}   = $_;

	for(my $i = $l[4] - $WINDOW; $i < $l[4] + $WINDOW; $i++){
            $truth{$l[3]}{$i} = \%load;
	}
        for(my $i = $l[6] - $WINDOW ; $i < $l[6] + $WINDOW ; $i++){
	    $truth{$l[3]}{$i} = \%load;
        }

	for(my $i = $l[9] - $WINDOW; $i < $l[9] + $WINDOW; $i++){
	    $truth{$l[3]}{$i} = \%load;
	}

    }
}

close $IN;

open(my $INB, '<', $fileb) or die "Can't open $fileb for reading\n$!\n";

my $tp    = 0;
my $fp    = 0;

my %HIT;
my $UNI;


SV: while (<$INB>) {
    chomp;
    $_ =~ s/\.//g;
    next SV if $_ =~ /^#/;
    my @l = split /\t/, $_;
    next SV if $_ =~ /BE=;/;

    my $end = $1 if $_ =~ /BE=(.*?);/;
    $end = $1 if $_ =~ /BE=(.*?)\t/;
    
    my $endSeqid = -1;
    my $endPos   = -1;


    if($_ !~ /BE=\.;/){
	my @alt = split /,/, $end;
	$endSeqid = $alt[0];
	$endPos = $alt[1];
    }

    if(defined $truth{$l[0]}{$l[1]} || defined $truth{$endSeqid}{$endPos}){

	my $skey = $l[0];
	my $pkey = $l[1];

	if(defined $truth{$endSeqid}{$endPos}){
	    $skey = $endSeqid;
	    $pkey = $endPos;
	}


	if(defined $HIT{$truth{$skey}{$pkey}{"ID"}}){
	    next SV;
	}


	print  $truth{$skey}{$pkey}{"TYPE"}, "\t";
	print  $truth{$skey}{$pkey}{"LENGH"}, "\t";
	print $WINDOW, "\t";
	print  $DEPTH, "\t";
	print  "1", "\t";
	print "WHAM", "\t";
	print  $truth{$skey}{$pkey}{"ID"};

	$tp ++;

	my $flagg = 0;
	if($l[9] =~ /1\/1/){
            $flagg = 1;
	}
	if($GENOTYPE){
	    print "\t$flagg\n";
	}
	else{
	    print "\n";
	}

#        print STDERR "$TYPE\t$DEPTH\tWHAM\t$flagg" , "\n";

	if($skey ne -1){
	    $HIT{$truth{$skey}{$pkey}{"ID"}} = 1;
	}
    }
    else{
#	print STDERR $_, "\n";
	$fp++;
    }
}



close $INB;

my $FDR = $fp / ($tp + $fp);

print STDERR "$TYPE\t$DEPTH\t$WINDOW\tWHAM\t$fp\t$tp\t$FDR\n";


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

