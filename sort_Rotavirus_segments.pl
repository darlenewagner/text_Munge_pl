use strict;
use Getopt::Long;
use warnings;

## Expects a multifasta of Rotavirus (or Influenza) segments
## chomp nucleotide sequence lines and sort fasta records by links

my $reference = '';

## Say '--reference Y' to create a spliced Rotavirus genome reference
GetOptions(
    'reference|r|' => \$reference,
    ) or die "Error in command line option argument\n";

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
    print $header, $sequence;
  }
else
  {
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
