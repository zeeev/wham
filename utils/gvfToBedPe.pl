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

perl gvfToBedPe.pl --buffer 10 --file stuff.gvf > stuff.bedpe 2> stuff.err

Description:

Converts a GVF to bedpe

Options:

file       - <STRING> - required - filename

-h,--help       - <FLAG>   - optional - print help statement and die
-b,--buffer     - <INT>    - optional - add basepair to both sides of the SV call [0]

Info:

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
			     'buffer=s'     => \$BUFFER
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

    GVF: while(my $line = <$FH>){
	chomp $line;
	next GVF if $line =~ /^#/;

	$svCount++;
	
	my @gvfLine = split /\t/, $line;
	
	my $svlen = $gvfLine[4] - $gvfLine[3];
	my $name = "$svCount:$gvfLine[2]:$svlen";
	
	my $bedline;
	
	my $start0 =   $gvfLine[3] - 1 - $BUFFER;
	my $start1 =   $gvfLine[3] - 1 + $BUFFER;

	my $end0   = $gvfLine[4] - 1 - $BUFFER;
	my $end1   = $gvfLine[4] - 1 + $BUFFER;

	$bedline .= "$gvfLine[0]";
	$bedline .= "\t$start0";
	$bedline .= "\t$start1";
	$bedline .= "\t$gvfLine[0]";
	$bedline .= "\t$end0";
	$bedline .= "\t$end1";
	$bedline .= "\t$name";
	$bedline .= "\t.";
	$bedline .= "\t.";
	$bedline .= "\t.";
	$bedline .= "\t$gvfLine[8]\n";
	
	print $bedline;
    }
    close $FH;
}

