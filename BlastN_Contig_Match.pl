#!/usr/bin/perl
use strict;
use Getopt::Long;
## Input BlastN file is tab-delimited with 'virus' or '10239' in 17th column
## Output format option from original BlastN command had 18 columns:
## -outfmt '6 qseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore sgi sacc staxids sblastnames stitle'

## Alphanumeric Genbank or RefSeq expected
my $accession = '';

## Two or more words expected, the first one strictly alphabetical
## Is empty string until assigned after $accession fails format test
my $virusName = '';

## For searching when $accession or $virusName are not supplied
my $highConsequenceStr = qr/(HIV\-(1|2)\sisolate\s.+,\scomplete\sgenome|Human\simmunodeficiency\svirus\s(1|2)\sproviral)/;
my @HighConsequenceArray = ();
my @SortingArray = ();

$accession = $ARGV[0];

my $help = '';
my $verbose = '';
my $sortIt = '';
my $advanced = '';

GetOptions(
    'help|h|'  => \$help,
    'verbose'  => \$verbose,
    'sort'     => \$sortIt,
    'highConsequence' => \$advanced,
    ) or die "Error in command line option arguments\n";

if($help)
   {
       print "\nUsage: BlastN_Contig_Match.pl [ACCESSION or \"VIRUS NAME STR.\"] --verbose --highConsequence < BlastN-Output.tsv\n";
       print "\n<------------------------------------------------- Options -------------------------------------------------->\n";
       print "ACCESSION or VIRUS NAME STR#, Search key (required)\n";
       print "--verbose,                    Output column headers (optional)\n";
       print "--sort,                       Sort output by MatchID (NCBI Accession)\n";
       print "--highConsequence,            Search for HIV-1 or other high-consequence pathogens (optional)\n";
       print "--help,                       Display this help message\n\n";
       exit;
   } 

if($accession =~ /^[A-Z][a-z]+\s+([A-Za-z][A-Za-z0-9]?[0-9]?\s)?([A-Za-z0-9]+(-|_)?[A-Za-z0-9]+|sp|str|strain)\.?\s+([A-Z]\s)?[A-Za-z0-9]+(-|_)?[A-Za-z0-9]+/)
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
      die "Search string too short and/or lacks sufficient information --$!";
  }




my %Sortable = {};
my $i = 0;
my @line = ();

while(<STDIN>)  ## Read BlastN file, use mode according to which string is nonempty
  {
      if(($accession =~ /^[A-Z][A-Z][A-Z0-9]+\.?[0-9]$/) || ($accession =~ /^N(C|T|X|Z)_[0-9]+\.?[0-9]$/))
        { 
	    @line = split(/\t/, $_);
	    if($line[16] =~ /virus|10239/)
	    {
		chomp $line[17];
		my $newLine = '';
		if($verbose)
		  {
		    $newLine = $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[14]."\t".$line[17]."\n";
		  }
		else
		 {
                     $newLine = $line[0]."\t\t".$line[1]."\t\t".$line[2]."\t\t".$line[14]."\t\t".$line[17]."\n";
		 }
		if($newLine =~ m/$accession/)
		  {
		     print $newLine;
		  }
		 elsif(($newLine =~ m/$highConsequenceStr/) && ($advanced))
	          {
		    push @HighConsequenceArray, $newLine;
		 }
	     }
	}   
      elsif($virusName =~ /^[A-Z][a-z]+\s+([A-Za-z][A-Za-z0-9]?[0-9]?\s)?([A-Za-z0-9]+(-|_)?[A-Za-z0-9]+|sp|str|strain)\.?\s+([A-Z]\s)?[A-Za-z0-9]+(-|_)?[A-Za-z0-9]+/)
        {
	    @line = split(/\t/, $_);
	    if($line[16] =~ /virus|10239/)
	    {
		chomp $line[17];
		my $newLine = '';
		if($verbose)
		 {
		    $newLine = $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[14]."\t".$line[17]."\n";
		 }
		else
		 {
                     $newLine = "'".$line[0]."'\t'".$line[1]."'\t'".$line[2]."'\t'".$line[14]."'\t'".$line[17]."'\n";
		 }
		if($newLine =~ m/$virusName/)
		{
                   print $newLine;     
		}
		elsif(($newLine =~ m/$highConsequenceStr/) && ($advanced))
	        {
                   push @HighConsequenceArray, $newLine;
		}
	    }
     	}
      elsif(($advanced) && ($virusName !~ /^[A-Z][a-z]+\s+[A-Za-z0-9]+(-|_)?[A-Za-z0-9]+\.?\s+[A-Za-z0-9]+(-|_)?[A-Za-z0-9]+/) && ($accession !~ /^[A-Z][A-Z][A-Z0-9]+\.?[0-9]$/) && ($accession !~ /^N(C|T|X|Z)_[0-9]+\.?[0-9]$/))
        {
	    @line = split(/\t/, $_);
      	    if($line[16] =~ /virus|10239/)
	    {
		chomp $line[17];
      	my $newLine = "'".$line[0]."'\t'".$line[1]."'\t'".$line[2]."'\t'".$line[14]."'\t'".$line[17]."'\n";
      		if($newLine =~ m/$highConsequenceStr/)
      		{
                    push @HighConsequenceArray, $newLine;
      		}
		elsif($sortIt)
		{
                    push @SortingArray, $newLine;
		}
      	    }
	}
      #else
      #{
      #	  die "Search string too short, lacks sufficient information, or --highConsequence option not provided...$!";
      #}
  }

if($advanced)
  {
      if($verbose)
        {
           print "\n<------- Search for High-Consequence Viral Pathogen ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      my $i = 0;
      foreach(@HighConsequenceArray)
      {
	  print $HighConsequenceArray[$i];
	  $i++;
      }
  }

