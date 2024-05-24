use strict;

## For output from vamb/bin/concatenate.py, count contigs for each strain, 'S'
## Append number after 'C' in fasta header, reset count for each new strain
## Uses an example of lookbehind assertion in the split function


my $fasta = $ARGV[0];

open(FASTA, $fasta) || die "Can't find fasta-formatted input file, $fasta $!";

my $strain = '';
my $prev_strain = '';
my $contigCount = 1;

while(<FASTA>)
{
    if($_ =~ /^>/)
    {
        $prev_strain = $strain;
	$strain = substr($_, 1, index($_, 'C') - 1);
	if($prev_strain ne $strain)
	 {
	     $contigCount = 1;
	 }
	my @newHeader = split(/(?<=\dC)/, $_);
	print $newHeader[0], $contigCount, "\n";
	$contigCount++;
    }
    else
    {
       print $_;
    }
}

close(FASTA)

