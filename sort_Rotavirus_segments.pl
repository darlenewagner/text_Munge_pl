use strict;
use Getopt::Long;
use warnings;

## Expects a multifasta of Rotavirus (or Influenza) segments
## chomp nucleotide sequence lines and sort fasta records by links

my $reference = '';
my $help = '';

## Say '--reference' to create a spliced Rotavirus genome reference
GetOptions(
    'reference|r|' => \$reference,
    'help|h|'      => \$help,
    ) or die "Error in command line option argument\n";

if($help)
{
    print "\nsort_Rotavirus_segments.pl is a script for sorting a multifasta file of Rotavirus or Influenza segments descending order of length.\n";
    print "INPUT is a single multifasta file as a commmand line argument. OUTPUT is a multifasta or single fasta file to STDOUT\n";
    print "##########\n";
    print "USAGE 1: Create a sorted multifasta from input multifasta -> \nsort_Rotavirus_segments.pl test_Rotavirus/Rotavirus_A_1976_G2P.4.fasta\n";
    print "### or ###\n";
    print "USAGE 2: Create a single fasta reference with 150-bp spacers separating segments ->\nsort_Rotavirus_segments.pl --reference test_Rotavirus/Rotavirus_A_1976_G2P.4.fasta\n";
    print "##########\n\n";

    die;
}

my $fasta = $ARGV[0];

open(FASTA, $fasta) || die "Can't find fasta-formatted input file, $fasta $!";

my @GENOME = <FASTA>;

my %Segments = {};

#my %sortedSegments {};

my @Titles = ();

my $i = 0;
my $prev = 0;

foreach(@GENOME)
{
     if($GENOME[$i] =~ /^>/)
     {
	 $prev = $i;
	 $Segments{$GENOME[$prev]} = '';
     }
     elsif($GENOME[$i] =~ /^[ACGNT]+/)
     {
	 chomp $GENOME[$i];
	 $Segments{$GENOME[$prev]} = $Segments{$GENOME[$prev]}.$GENOME[$i];
     }

     $i++;
 }

my @sortedKeys = sort {length($Segments{$b}) <=> length($Segments{$a})} keys(%Segments);


my $spacer = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";


if($reference)
{
    my $header = '';
    my $sequence = '';
    
    for(my $count = 0; $count < scalar @sortedKeys; $count++)
      {
	  if($count == 0)
	    {
	      $header = $sortedKeys[$count];
  	      $sequence = $Segments{$sortedKeys[$count]}.$spacer;
	    }
	  elsif($count < 10)
	  {
	      $sequence = $sequence.$Segments{$sortedKeys[$count]}.$spacer;
	  }
	  else
	  {
             $sequence = $sequence.$Segments{$sortedKeys[$count]}."\n";
	  }
      }
    #print $header, $sequence;

    print $header;
    for(my $line = 0; $line < length($sequence); $line = $line + 60)
      {
        print substr($sequence, $line, 60), "\n";
      }
    
  }
else
{
    ## output sorted segments, longest to shortest, without concatenating or deleting fasta headers
  my $k = 0;    
  foreach(@sortedKeys)
   {
     print $sortedKeys[$k];
     for(my $line = 0; $line < length($Segments{$sortedKeys[$k]}); $line = $line + 60)
      {
        print substr($Segments{$sortedKeys[$k]}, $line, 60), "\n";
      }
     #print $Segments{$sortedKeys[$k]}, "\n";
     $k++;
    }
   print "\n";
 }
