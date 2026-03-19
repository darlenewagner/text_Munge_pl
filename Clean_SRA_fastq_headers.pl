#!/usr/bin/perl
use strict;
use Getopt::Long;
use warnings;

my $strand = '1:N:0:1';

GetOptions(
    'strand|s|=s' => \$strand,
    ) or die "Error in command line option argument\n";

if($strand !~ /^(1|2)\:/)
{
    die "Error: The first integer in '--strand' argument must be 1 or 2\n";
}


while(<STDIN>)
{
    if(($. % 4) == 1)
     {
	my $ws = index($_, ' '); 
	my $id = substr($_, 0, $ws);
	print $id, " ", $strand, "\n";
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
