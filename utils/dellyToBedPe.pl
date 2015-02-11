#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw( max sum );

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

perl dellyToBedPe.pl --append --file delly.vcf > delly.bed 2> delly.err

Description:

Converts a delly VCF into the BEDPE format.

Options:

file       - <STRING> - required - filename

-h,--help       - <FLAG>   - optional - print help statement and die
-a,--append     - <FLAG>   - optional - concatenate WHAM name to SV annotations (bcbio compatible)
-b,--buffer     - <INT>    - optional - add basepair to both sides of the SV call [0]

Info:

-a -- append -- if wham calls are classified this addeds the annotation
-b -- buffer -- add x bp of slop in both directions to both coordinates


";


my ($help);
my $FILE;
my $APPEND        = 0;
my $BUFFER        = 0;
my $ANNOTATED_FLAG = 1;

my $opt_success = GetOptions('help'         => \$help,
			     'file=s'       => \$FILE,
			     'buffer=s'     => \$BUFFER,
			     'append'       => \$APPEND
		      );

die $usage if $help || ! $opt_success;

if(! defined $FILE){
    print STDERR "\nFATAL: file was not provided.\n ";
    die $usage;
}


processLines();

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub processLines{

    open(my $FH, "<", $FILE) || die "FATAL: could not open $FILE for reading\n";    

    my $svCount = 0;
    
    VCF: while(my $line = <$FH>){
	chomp $line;
	next VCF if $line =~ /^#/;
	next VCF if $line =~ /LowQual/;
	
	my @vcfLine = split /\t/, $line;


	$svCount++;
	
	my $fivePrimeChr  = $vcfLine[0];
	
	$line =~ /CHR2=(.*?);/;
	
	my $threePrimeChr = $1; 

	#one based to zerobased
	my $startPos0 = $vcfLine[1]   - 1 ;
	my $startPos1 = $vcfLine[1]   - 1 ;
	$line =~ /;END=(.*?);/;
	my $endPos0   = $1  - 1 ;
	my $endPos1   = $1  - 1 ;

	$line =~ /;CIPOS=(.*?);/;

	my @CIPOS = split /,/, $1;

	$line =~ /;CIEND=(.*?);/;

	my @CIEND = split /,/, $1;
	    
	$startPos0 -= $BUFFER + $CIPOS[0] ;
	$startPos1 += $BUFFER + $CIPOS[1] ;
	$endPos0   -= $BUFFER + $CIEND[0] ;
	$endPos1   += $BUFFER + $CIEND[1] ;

	#wham reports both positional orders 
	
	if($fivePrimeChr eq $threePrimeChr){
	    if($startPos0 > $endPos1){
		my $b0 = $startPos0;
		my $b1 = $startPos1;
		$startPos0 = $endPos0;
		$startPos1 = $endPos1;
		$endPos0 = $startPos0;
		$endPos1 = $startPos1;
	    }
	}

	my $bedline .= "$fivePrimeChr"    ;
	   $bedline .= "\t$startPos0"     ;
           $bedline .= "\t$startPos1"     ;
	   $bedline .= "\t$threePrimeChr" ;
	   $bedline .= "\t$endPos0"       ;
	   $bedline .= "\t$endPos1"       ;

	$line =~ /SVTYPE=(.*?);/;
	my $SVTYPE = $1;
	$line =~ /SVLEN=(.*?);/;
	my $SVLEN =  $1;
	

	if($ANNOTATED_FLAG){
	
	    if($APPEND){
		$bedline .= "\t$svCount:delly\_$SVTYPE:$SVLEN";
	    }
	    else{
		$bedline .= "\t$svCount:$SVTYPE:$SVLEN";
	    }
	    
	 

	   $bedline .= "\t0";

	}
	else{
	    $bedline .= "\t$svCount:.:$SVLEN\t.";
	}
	
	$bedline .= "\t.";
	$bedline .= "\t.";
	$bedline .= "\t$vcfLine[6]";
	$bedline .= "\t$vcfLine[7]";

	print "$bedline\n";
	
    }
    close $FH;
}

