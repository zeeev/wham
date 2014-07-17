#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Fasta;
use Text::Wrap;
$Text::Wrap::columns = 72;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

generateSV.pl --number 10 --type DELETION --original 

Description:

Options:
    -n,number    number of sequences to make
    -t,type      type of SV [INSERTION|DELETION|DUPLICATION]
    -o,original  the template fasta


";


my ($help);
my $number = 10;
my $type       ;
my $original   ;
my $opt_success = GetOptions('help'             => \$help,
                             'type=s'           => \$type,
			     'number=s'         => \$number,
                             'original=s'	=> \$original,
		      );


die $usage unless $number && $original;


my $db = Bio::DB::Fasta->new($original);

my ($seqid, @s) = $db->ids();
my $seq         = $db->seq($seqid);

my $SEQL = length($seq) - 1;
my @SEQS;

map{push @SEQS, $seq} (1..10);

DELETION();

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub DELETION{

    my $trueL = int(rand(9000) + 1000);
    my $trueP = int(rand($SEQL - 1000) + 1000);

    my $index = 0;
		 
    my $uni = sprintf "%08X\n", rand(0xffffffff);
    
    foreach my $seq (@SEQS){
	$index += 1;
	my $start  = $trueP + int(rand(100) - 100);
	my $end    = $trueP + $trueL + int(rand(100) - 100);
	my $length = $end - $start; 
	substr($seq, $start, $end-$start) = "";
	print ">Sample-$index-DELETION-$trueP-$start-$end-$length-$uni";
	print wrap('','', "$seq\n");
    }
		 
}
    
