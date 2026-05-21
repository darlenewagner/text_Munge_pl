use strict;

## split fasta header at input string, $trimmer, to shorten header in fasta input file, $fasta

my $fasta = $ARGV[0];

open(FASTA, $fasta) || die "Can't find fasta-formatted input file, $fasta $!";

while(<FASTA>)
{
    if($_ =~ /^>/)
    {
	print $_;
    }
    else
    {
	my $cleanedLine = $_;
	$cleanedLine =~ tr/WRYS/NNNN/;
       print $cleanedLine;
    }
}

close(FASTA)

