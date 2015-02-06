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

perl whamToBedPe.pl --append --file wham.vcf > wham.bed 2> wham.err

Description:

Converts a classified or unclassified WHAM VCF into the BEDPE format.
If the WC and WP fields are not present in the VCF header whamToBed does not
provide the type of structural variant.

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
my $PAIRED        = 0;
my $BUFFER        = 0;
my $SINGLE_BUFFER = 0;

my $opt_success = GetOptions('help'         => \$help,
			     'file=s'       => \$FILE,
			     'buffer=s'     => \$BUFFER,
			     'append'       => \$APPEND,
			     'singletons=s' => \$SINGLE_BUFFER,
			     'paired'       => \$PAIRED	     
		      );

die $usage if $help || ! $opt_success;

if(! defined $FILE){
    print STDERR "\nFATAL: file was not provided.\n ";
    die $usage;
}

my $ANNOTATED_FLAG = 0;

checkForAnnotations();

if($ANNOTATED_FLAG){
    print STDERR "INFO: BED will include SV type and annotation score\n";
}
else{
    print STDERR "INFO: BED will NOT include SV type and annotation score\n";
    print STDERR "      Classify WHAM ouput to get these features in the bed\n";
}

processLines();

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

#Checking if the WHAM vcf has been annotated by looking for the WC field in 
#the VCF file.

sub checkForAnnotations {
    
    open(my $FH, "<", $FILE) || die "FATAL: could not open $FILE for reading\n";

    my $headerFlag = 1;
    while(my $line = <$FH>){
	if($line =~ /ID=WC/){
	    $ANNOTATED_FLAG = 1;
	    last;
	}
	if($line !~ /^#/){
	    last;
	}
    }
    close $FH;
}

#-----------------------------------------------------------------------------

sub processLines{

    open(my $FH, "<", $FILE) || die "FATAL: could not open $FILE for reading\n";    

    my $svCount = 0;
    
    VCF: while(my $line = <$FH>){
	chomp $line;
	next VCF if $line =~ /^#/;
	
	my @vcfLine = split /\t/, $line;
	
	my %info = map{ split /=|;/} $vcfLine[7];

	#skipping unpaired breakpoints
	next VCF if $info{"BE"} eq '.';

	$svCount++;

	my @otherBreak = split /,/, $info{"BE"};
	
	my $fivePrimeChr  = $vcfLine[0];
	my $threePrimeChr = $otherBreak[0]; 

	#one based to zerobased
	my $startPos0 = $vcfLine[1]    - 1 ;
	my $startPos1 = $vcfLine[1]    - 1 ;
	my $endPos0   = $otherBreak[1] - 1 ;
	my $endPos1   = $otherBreak[1] - 1 ;

	$startPos0 -= $BUFFER;
	$startPos1 += $BUFFER;
	$endPos0   -= $BUFFER;
	$endPos1   += $BUFFER;

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

	if($ANNOTATED_FLAG){

	    if($APPEND){
		$bedline .= "\t$svCount:wham\_$info{WC}:$info{SVLEN}";
	    }
	    else{
		$bedline .= "\t$svCount:$info{WC}:$info{SVLEN}";
	    }
	    
	    my @probs = split /,/, $info{WP};
	    
	    my $maxP = max @probs;

	    $bedline .= "\t$maxP";

	}
	else{
	    $bedline .= "\t$svCount:.:$info{SVLEN}\t.";
	}
	
	$bedline .= "\t.";
	$bedline .= "\t.";
	$bedline .= "\t$vcfLine[6]";
	$bedline .= "\t$vcfLine[7]";

	print "$bedline\n";
	
    }
    close $FH;
}

