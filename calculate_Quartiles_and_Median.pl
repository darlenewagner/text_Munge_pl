use strict;
use Statistics::Lite;

## split fasta header at input string, $trimmer, to shorten header in fasta input file, $fasta


while(<FASTA>)
{
    if($_ =~ /^>/)
    {
	my @newHeader = split $trimmer, $_;
	print $newHeader[0], "\n";
    }
    else
    {
       print $_;
    }
}

close(FASTA)

