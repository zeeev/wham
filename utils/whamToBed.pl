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

perl whamToBed.pl --append --file wham.vcf > wham.bed 2> wham.err

Description:

Converts a classified or unclassified WHAM VCF into the BED format.
If the WC and WP fields are not present in the VCF header whamToBed does not
provide the type of structural variant.

Options:

file       - <STRING> - required - filename

-h,--help       - <FLAG>   - optional - print help statement and die
-a,--append     - <FLAG>   - optional - concatenate WHAM name to SV annotations (bcbio compatible)
-p,--paired     - <FLAG>   - optional - output paired breakpoints only (increase specificity)
-b,--buffer     - <INT>    - optional - add basepair to both sides of the SV call [0]
-s,--singletons - <INT>    - optional - add basepair to only unpaired breakpoints [0]

Info:

paired - This option filters out sites were there was no split read support for 
         the other side of the breakpoint. 

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

    VCF: while(my $line = <$FH>){
	chomp $line;
	next VCF if $line =~ /^#/;

	my @vcfLine = split /\t/, $line;
	
	my %info = map{ split /=|;/} $vcfLine[7];

	my $startPos = $vcfLine[1] - 1 ;
	my $endPos = $vcfLine[1]   - 1 ;

	if($info{'BE'} ne '.'){
	    $startPos -= $BUFFER;
	    $endPos   += $BUFFER;
	}

	if($info{"BE"} eq '.'){
	    next VCF if $PAIRED;
	    $startPos -= $SINGLE_BUFFER;
	    $endPos   += $SINGLE_BUFFER;
	}
	else{
	    my @end = split /,/, $info{BE};
	    if($end[0] ne $vcfLine[0]){
		print STDERR "WARNING: translocations cannot be represented 
                in bed, skipping: $line\n";
		next VCF;
	    }
	    else{
		$endPos = $end[1];
	    }
	}

	if($startPos > $endPos){
	    my $temp = $startPos;
	    $startPos = $endPos;
	    $endPos = $temp;
	}
	
	my $bedLine = "$vcfLine[0]\t$startPos";
	$bedLine .= "\t$endPos";

	if($ANNOTATED_FLAG){
	    if($APPEND){
		$bedLine .= "\twham\_$info{WC}";
	    }
	    else{
		$bedLine .= "\t$info{WC}";
	    }
	    
	    my @probs = split /,/, $info{WP};
	    
	    my $maxP = max @probs;
	    my $sumP = sum @probs;

	    my $normP = $maxP / $sumP;

	    $bedLine .= "\t$normP";

	}
	else{
	    $bedLine .= "\t$vcfLine[0]:$startPos:$endPos\t\.";
	}

	print "$bedLine\n";
	
    }
    close $FH;
}

