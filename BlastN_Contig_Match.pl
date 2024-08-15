use strict;
use Getopt::Long;
## Assume input BlastN output file is tab-delimited with 'virus' or '10239' in 17th column

## Alphanumeric Genbank or RefSeq expected
my $accession = '';

## Two or more words expected, the first one strictly alphabetical
my $virusName = '';

$accession = $ARGV[0];

my $help = '';
my $verbose = '';
my $advanced = '';

GetOptions(
    'help'     => \$help,
    'verbose'  => \$verbose,
    'advanced' => \$advanced,
    );



if($accession =~ /^[A-Z][a-z]+\s+[A-Za-z0-9]+\.?\s+[A-Za-z0-9]+/)
  {
    $virusName = $ARGV[0];
     if($verbose)
        {
	    print "Search by species/strain name:\n";
	    print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
  }
elsif($accession =~ /^[A-Z][A-Z][A-Z0-9]+\.?[0-9]$/)
  {
    if($verbose)
      {
	  print "Search by GenBank accession:\n";
	  print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
      }
  }
elsif($accession =~ /^N(C|T|X|Z)_[0-9]+\.?[0-9]$/)
  {
     if($verbose)
       {
	  print "Search by RefSeq accession:\n";
	  print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
       }
  }
elsif(!$advanced)
  {
    die "Search string too short and/or lacks sufficient information. $!"
  }

## For searching when $accession or $virusName are not supplied
my $roguesGallery = qr/HIV\-(1|2)\sisolate\s.+,\scomplete\sgenome/;

my @multi_segment = ();
my $i = 0;
my @line = ();

while(<STDIN>)  ## Read BlastN file, use mode according to which string is nonempty
  {
      if(($accession =~ /^[A-Z][A-Z][A-Z0-9]+\.?[0-9]$/) || ($accession =~ /^N(C|T|X|Z)_[0-9]+\.?[0-9]$/))
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
      elsif($virusName =~ /^[A-Z][a-z]+\s+[A-Za-z0-9]+\.?\s+[A-Za-z0-9]+/)
        {
	    @line = split(/\t/, $_);
	    if($line[16] =~ /virus|10239/)
	    {
		my $newLine = $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[14]."\t".$line[17];
		if($newLine =~ m/$virusName/)
		{
                   print $newLine;     
		} 
	    }
     	}
      elsif(($_ =~ $roguesGallery))
        {
	    @line = split(/\t/, $_);
	    if($line[16] =~ /virus|10239/)
	    {
		my $newLine = $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[14]."\t".$line[17];
		if($newLine =~ $roguesGallery)
		{
                    print $newLine;
		}
	    }
	}
  }

