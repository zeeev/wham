#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

truth.pl -t INR|INV|DUP|DEL -w 25 truth.events lumpy.bedbe 

Description:

-t <STRING> --required-- type defined by SVsim.
-w <INT>    --optional-- the number of bases on either side of the true event to count as acceptable [25]
-d <INT>    -- depth  -- the depth of bam file used for simulation [10]

";

my ($help);
my $TYPE;
my $WINDOW = 25;
my $DEPTH  = 10;
my $opt_success = GetOptions('help'     => \$help,
			     "window=s" => \$WINDOW,
			     "depth=s"  => \$DEPTH,
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

    my @l = split /\t/, $_;

    if($TYPE eq "DEL"){
	
	#DEL0022 DEL     50000   8       107462299       8       107462299       +       8       107512298       +
    
	my %load ;
	$load{"ID"}    = $l[0];
	$load{"TYPE"}  = $l[1];
	$load{"LENGH"} = $l[2];
	$load{"CHR"}   = $l[3];
	$load{"START"} = $l[4];
	$load{"END"}   = $l[9];

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

	my %load;

	$load{"ID"}     = $l[0];
        $load{"TYPE"}   = $l[1];
        $load{"LENGH"}  = $l[2];
        $load{"CHRL"}   = $l[3];
        $load{"CHRR"}   = $l[5];
	$load{"START"}  = $l[4];
	$load{"STARTR"} = $l[6];
        $load{"ENDR"}   = $l[9];

	#right at the spot

	for(my $i = $l[4] - $WINDOW; $i < $l[4] + $WINDOW; $i++){
	    $truth{$l[3]}{$i} = \%load;
        }

	for(my $i = $l[4] - $WINDOW - int($l[2] / 2); $i < $l[4] + $WINDOW - int($l[2] / 2); $i++){
            $truth{$l[3]}{$i} = \%load;
        }

	for(my $i = $l[4] - $WINDOW + int($l[2] / 2); $i < $l[4] + $WINDOW + int($l[2] / 2); $i++){
            $truth{$l[3]}{$i} = \%load;
	}


	# upstream

	for(my $i = $l[4] - $WINDOW + $l[2]; $i < $l[4] + $WINDOW + $l[2]; $i++){
            $truth{$l[3]}{$i} = \%load;
        }

	#downstream

	for(my $i = $l[4] - $WINDOW - $l[2]; $i < $l[4] + $WINDOW - $l[2]; $i++){
            $truth{$l[3]}{$i} = \%load;
        }


	# other chr left

	for(my $i = $l[6] - $WINDOW; $i < $l[6] + $WINDOW ; $i++){
            $truth{$l[5]}{$i} = \%load;
        }

	#other chr upstream

	for(my $i = $l[6] - $WINDOW + $l[2]; $i < $l[6] + $WINDOW + $l[2]; $i++){
            $truth{$l[5]}{$i} = \%load;
        }

	for(my $i = $l[6] - $WINDOW - $l[2]; $i < $l[6] + $WINDOW - $l[2]; $i++){
            $truth{$l[5]}{$i} = \%load;
        }

	

	for(my $i = $l[9] - $WINDOW; $i < $l[9] + $WINDOW; $i++){
            $truth{$l[5]}{$i} = \%load;
        }


	for(my $i = $l[9] - $WINDOW - $l[2]; $i < $l[9] + $WINDOW - $l[2]; $i++){
            $truth{$l[5]}{$i} = \%load;
        }

	for(my $i = $l[9] - $WINDOW + $l[2]; $i < $l[9] + $WINDOW + $l[2]; $i++){
            $truth{$l[5]}{$i} = \%load;
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

SV: while (<$INB>) {

#1       7017630 7017662 1       8017622 8017664 1       4.15492e-27     +       -       TYPE:DELETION   IDS:1,9;2,3     STRANDS:+-,12   MAX:1:7017643;1:8017664 95:1:7017641-7017644;1:8017664-8017664
    chomp;
    $_ =~ s/\.//g;
    next if $_ =~ /^#/;
    my @l = split /\t/, $_;

#MAX:1:7017643;1:8017664

    my @max = split /:|;/, $l[13];

    shift @max;

    if(defined $truth{$max[0]}{$max[1]} || defined $truth{$max[2]}{$max[3]}){

        my $seqid = 0;
        my $key   = 1;

        if(defined $truth{$max[0]}{$max[1]}){
            $seqid = 0;
            $key   = 1;
        }
        if(defined $truth{$max[2]}{$max[3]}){
            $seqid = 2;
            $key   = 3;
        }

	if(defined $HIT{$truth{$max[$seqid]}{$max[$key]}}){
            next SV;
        }
    	
	print  $truth{$max[$seqid]}{$max[$key]}{"TYPE"}, "\t";
	print  $truth{$max[$seqid]}{$max[$key]}{"LENGH"}, "\t";
	print  $WINDOW, "\t";
	print  $DEPTH, "\t";
	print  "1", "\t";
	print  "LUMPY", "\t";
	print  $truth{$max[$seqid]}{$max[$key]}{"ID"}, "\n";

	$tp ++;

	$HIT{$truth{$max[$seqid]}{$max[$key]}} = 1;
#	$HIT{$truth{$l[$seqid]}{$l[$key]}{"ID"}} = 1;
    
    }
    else{
#	print STDERR $_ , "\n";
	$fp++;
    }
}



close $INB;

my $FDR = $fp / ($tp + $fp);

print STDERR "$TYPE\t$DEPTH\t$WINDOW\tLUMPY\t$fp\t$tp\t$FDR\n";


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

