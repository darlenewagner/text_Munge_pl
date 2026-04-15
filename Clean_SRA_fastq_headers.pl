#!/usr/bin/perl
use strict;
use Getopt::Long;
use warnings;

my $strand = '1:N:0:1';

my $singleEnd = 0;

my $result = GetOptions(
    'strand:s' => \$strand,
    ) or die "Error in command line option argument\n";

#print $result, "\n";

if(($strand !~ /^(1|2)\:/) || ($strand !~ /^SINGLE/))
{
    die "Error: The first integer in '--strand' argument must be 1 or 2\n";
}

if($strand !~ /^SINGLE/)
{
    $singleEnd = 1;
}


while(<STDIN>)
{
    if(($. % 4) == 1)
     {
	my $ws = index($_, ' '); 
	my $id = substr($_, 0, $ws);
	if($singleEnd == 0){
	    print $id, " ", $strand, "\n";
	}
     }
    elsif(($. % 4) == 3)
     {
	 my $first = substr($_, 0, 1);
	 print $first, "\n";
     }
    else
     {
	print;
     }
}
