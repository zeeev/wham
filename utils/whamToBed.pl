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

file   - <STRING> - required - filename
help   - <FLAG>   - optional - print help statement and die
append - <FLAG>   - optional - concatenate WHAM name to SV annotations (bcbio compatible)
paired - <FLAG>   - optional - output paired breakpoint only

";


my ($help);
my $FILE;
my $APPEND = 0;
my $PAIRED = 1;

my $opt_success = GetOptions('help'      => \$help,
			     'file=s'    => \$FILE,
			     'append'    => \$APPEND,
			     'paired'    => \$PAIRED
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

	my $bedLine = "$vcfLine[0]\t$vcfLine[1]";
	
	my $endPos = $vcfLine[1];

	if($info{BE} eq 'nan'){
	    next VCF if $PAIRED;
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
	    $bedLine .= "\t$vcfLine[0]:$vcfLine[1]:$endPos\t\.";
	}

	print "$bedLine\n";
	
    }
    close $FH;
}

