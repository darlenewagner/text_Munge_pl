use strict;
## Assume input BlastN output file is tab-delimited with 'virus' or '10239' in 17th column

## Alphanumeric Genbank or RefSeq expected
my $accession = 'EU234061';

## Two or more words expected, the first one strictly alphabetical
my $virusName = '';

## For searching when $accession or $virusName are not supplied
my $roguesGallery = 'HIV-1';

my @multi_segment = ();
my $i = 0;
my @line = ();

while(<STDIN>)  ## Read BlastN file, use mode according to which string is nonempty
  {
      if(($accession =~ /[A-Z][A-Z][A-Z0-9]+/) || ($accession =~ /^N(C|T|X|Z)_[A-Z0-9]+/))
        { 
	    @line = split(/\t/, $_);
	    if($line[16] =~ /virus|10239/)
	     {
		 my $newLine = $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[14]."\t".$line[17];
		 if($newLine =~ m/$accession/)
		 {
		     print $newLine;
		 }
	     }
	 }   
      elsif($virusName =~ /^[A-Z]/)
        {
	    @line = split(/,\s/, $_);
	}
  }

