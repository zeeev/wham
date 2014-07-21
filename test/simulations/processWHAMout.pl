#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

processWHAMout.pl wham-TYPE-truepos.out

Description:


";


my ($help);
my $opt_success = GetOptions('help'    => \$help,
			      );

die $usage if $help || ! $opt_success;
#DELETION-_832969.wham.out
my $file = shift;

$file =~ /_([0-9]+)_.*\.wham/;

my $pos = $1;

print STDERR "INFO: position: $pos\n";

die $usage unless $file;
open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

LINE: while (<$IN>) {
    chomp;
    next LINE if $_ =~ /^\#/;
    
    $_ =~ /LRT=(.*?)\s+/;
    my $lrt = $1;

    my @l = split /\t/, $_;

    my $positive = 0;

    if(abs($l[1] - $pos) <= 500){
	$positive = 1;
    }
    print $l[1], "\t", $lrt, "\t", $positive, "\t", $file, "\t", $pos, "\n";
	
}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub {

}

