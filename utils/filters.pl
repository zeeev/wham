#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

filters.pl -w <STRING>

Description:

This script provides a set of heuristic filters designed to reduce false positives

Options:

--wham, -w    <STRING> - WHAM file to be filtered
--depth,-d    <BOOL>   - To remove sites where depth is 2 times the average across indivduals 

";

my ($help);
my $wham ;
my $depth;
my $opt_success = GetOptions('help'    => \$help ,
			     'wham=s'  => \$wham ,
			     'depth'   => \$depth,
    );

die $usage if $help || ! $opt_success || ! $wham;


open (my $IN, '<', $wham) or die "Can't open $wham for reading\n$!\n";

my $depth_sum       = 0; 
my $number_of_sites = 0;
my $avg_depth       = 0;

LINE: while (<$IN>) {
    chomp;
    next LINE if $_ =~ /^#/;
    $number_of_sites ++;

    my @l = split /\s+/, $_;

    for(my $i = 10; $i < $#l ; $i++){
	my @gt = split /:/, $l[$i];
	$depth_sum += $gt[-1];
	$number_of_sites ++;
    }
    last LINE if $number_of_sites > 100_000;
}

$avg_depth = $depth_sum / $number_of_sites;

close $IN;

print STDERR "number of sites assayed                         : ", $number_of_sites, "\n";
print STDERR "average depth across individuals                : ", $avg_depth      , "\n";
print STDERR "now filtering sites where average depth is above: ", $avg_depth * 10;

open (my $INB, '<', $wham) or die "Can't open $wham for reading\n$!\n";

LINEB: while(<$INB>){
    chomp;
    if($_ =~ /^#/){
	print $_;
	print "\n";
	next LINEB;
    }
    
    my $site_depth_sum = 0;
    my $number_of_cov  = 0;
    
    my @l = split /\s+/, $_;

    for(my $i = 10; $i < $#l ; $i++){
        my @gt = split /:/, $l[$i];
        $site_depth_sum += $gt[-1];
	if($gt[0] ne './.'){
	    $number_of_cov++;
	}
    }

    $_ =~ /LRT=(.*?);/;
    my $flag = "PASS"; 

    while( $flag  eq "PASS"){

	if($1 < 1){
	    $flag = "LOW LRT";
	    last;
	}

	if($number_of_cov < 2){
	    $flag = "LOW GC";
	    last;
	}
	my $asite = $site_depth_sum / $number_of_cov;
	
	if($asite > $avg_depth * 10){

	    $flag = "HIGH DEPTH: $asite";   
	    last;
	}
	last;
    }

    if($flag ne "PASS"){
	print STDERR "$flag : $l[0]-$l[1]:$l[7]\n";    
	next LINEB;
    }

    print $_;
    print "\n";

}


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------


