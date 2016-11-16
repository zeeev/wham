#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "

Synopsis:

cat whamg.vcf | perl filtWhamG.pl > whamg.filt.vcf

Description:

   This set of filters is for high coverage genomes.
   Below is a list of filters.

   1. Removes SVs less than 50Bp or greater than 1Mb.

   2. Removes <DEL>s that don't have \"TF\" greater than 3.
      - This will remove real small deletions.

   3. Removes <DUP>s that do not have \"U\" greater than 3.

   4. Removes <INV>s that do not have \"V\" greater than 3.

   5. Removes SVs that have a max(CW) < 0.2.  These calls
      are poorly categorized. The dispersion of the weight
      suggests they are just poor mapping regions.

  6.  Removes SVs that have more than 0.2 CW going to a translocation.
      Calls with high cross seqid mapping are often false positives.

  7.  Requires at least one individual joint called has high support SP > 7.
      This removes calls that have low support across many individual.

";

my ($help);
my $opt_success = GetOptions('help'    => \$help,
			      );

die $usage if $help || ! $opt_success;

my %lookup;

 $lookup{"DEL"} = 0;
 $lookup{"DUP"} = 1;
 $lookup{"INV"} = 2;
 $lookup{"INS"} = 3;
 $lookup{"BND"} = 4;


while (<STDIN>) {
    chomp;
    if($_ =~ /^\#/){
        print "$_\n";
    }
    else{

        my @line = split /\t/, $_;

        my $len = $1 if $_ =~ /SVLEN=(.*?);/;
        my $pass = 1;
        $pass = 0 if abs($1) < 50 || abs($1) > 1000000;

        my $type = $1 if $_ =~ /;SVTYPE=(.*?);/;

        my $tf = $1 if $_ =~ /;TF=(.*?);/;
        my $u  = $1 if $_ =~ /;U=(.*?);/;
        my $v  = $1 if $_ =~ /;V=(.*?)\tGT/;
        if($type eq "DEL" && $tf < 4){
            $pass = 0;
        }
        if($type eq "DUP" && $u < 4){
            $pass = 0;
        }
        if($type eq "INV" && $v < 4){
            $pass = 0;
        }

        my @cw = split ",", $1 if $_ =~ /;CW=(.*?);/;
        $pass = 0 if $cw[$lookup{$type}] < 0.2;
        #these are translocation-ish
        $pass = 0 if $cw[-1] > 0.2;

        my $enough = 0;

        my $n = 0;

        for(my $i = 9; $i < scalar @line ; $i++){
            my @gt = split /:/, $line[$i];

            if($gt[-1] > 0){
                $n++;
            }

            if($gt[-1] > 7){
                $enough = 1;
            }
        }

        my $min = $n * 1.5;

        $_ =~ /A=(.*?);/;

        $pass = 0 if $1 < $min;

        $pass = 0 if $enough == 0;

        print "$_\n" if $pass == 1;
    }
}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
