#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

quick_human.pl <How the hell do you use this thing>

Description:

No really, how the hell do you use this thing!

";


my ($help);
my $opt_success = GetOptions('help'    => \$help,
			      );

die $usage if $help || ! $opt_success;

open (my $IN, '<', "human.hg19.build.txt") or die "Can't open human.hg19.build.txt for reading\n$!\n";


while (<$IN>) {
    chomp;
    my @l = split /\s+/, $_;
    for(my $i = 0; $i < $l[1]; $i += 10_000_000){
	my $end = $i + 10_000_000;

	print "WHAM-BAM -r $l[0]:$i-$end -t \$TARG -b \$BACK > tmp-$l[0]-$i-$end.wham.vcf\n" 
    }
}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub {

}

