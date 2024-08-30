#!/usr/bin/perl
use strict;
use Getopt::Long;
# Prepare metagenomic profile data for presentation in .html via Krona/2.7
# Create krona-readable metagenomic profiles in seven-column:
#name    taxID   taxRank genomeSize      numReads        numUniqueReads  abundance
# or four-column format:
#name    taxID    taxRank    abundance

my $infile = $ARGV[0];

my $help = '';
my $verbose = '';
my $seven = '';
my $simple = '';

GetOptions(
    'help|h|'  => \$help,
    'verbose|v|'  => \$verbose,
    'sevenColumn|c|' => \$seven,
    'simpleTax|s|' => \$simple,
    ) or die "Error in command line option arguments\n";

if($help)
   {
       print "\nUsage: CAMI_format_to_Krona.pl --verbose --sevenColumn --simpleTax  MetAlgin_or_MetaPhlAn.format.file\n";
       print "\n<------------------------------------------- Options -------------------------------------------->\n";
       print "--verbose,                Show original column names in STDERR\n";
       print "--sevenColumn,            Output in seven-column format (four-column default)\n";
       print "--simpleTax,              Remove redundant strain information\n";
       print "--help,                   Display this help message\n\n";
       exit;
   }


my @line = ();

open(INFILE, $infile) || die "Can't find CAMI-formatted input file, $infile $!";

if($seven)
  {
     print "name\ttaxID\ttaxRank\tgenomeSize\tnumReads\tnumUniqueReads\tabundance\n";
  }
else
  {
      print "name\ttaxID\ttaxRank\tabundance\n";
  }    

while(<INFILE>)
  {
     if($_ !~ /^(@|#|\s)/)
     {
	 #print $_;
	 @line = split(/\t/, $_);
	 my @taxPath = split(/\|/, $line[3]);
         my $pathlen = scalar @taxPath - 1;
	 chomp $line[4];
	 if($simple)
	 {
	     my @cleanTax = split(/\./, $line[0]);
	     
	     if($taxPath[$pathlen] !~ /\sunknown\sstrain/)
	     {	 
		 ## suppress print if strain data is redundant
		 if($line[1] ne 'strain')
	         {
      	           print $taxPath[$pathlen], "\t", $cleanTax[0], "\t", $line[1], "\t", $line[4], "\n";
	         }
	     }
	   }
	 else
	   {
	       print $taxPath[$pathlen], "\t", $line[0], "\t", $line[1], "\t", $line[4], "\n";
	   }
     }
  }

