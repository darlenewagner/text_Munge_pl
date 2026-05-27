use strict;
use Statistics::Lite;

## Read results of 'samtools depth my_reads.sorted.bam' into STDIN to compute
## Lower Quartile, Median, and Upper Quartile of read mapping depths

my @array = ();
my ($lower, $median, $upper) = 0;

while(<STDIN>)
  {
     if($_ =~ /[A-Z0-9_\.]+\s+[0-9]+\s+[0-9]+/)
     {
   	 my @line = split(/\s+/, $_);
	 push @array, $line[2];
     }
     else
     {
   	 die "Unexpected input format!\n";
     }
  }

## to guarantee numerical sort with numerical output
my @sorted_numbers = sort { $a <=> $b } @array;

my $arrLength = scalar @sorted_numbers;

 my @lowerArray = ();
 my @upperArray = ();

 if(($arrLength % 2) == 1){
    my $ii = int($arrLength / 2);  
    $median = $sorted_numbers[$ii];

    for(my $y = 0; $y < $ii; $y++){
	push @lowerArray, $sorted_numbers[$y];
    }

    for(my $z = $ii + 1; $z < scalar @sorted_numbers; $z++){
	push @upperArray, $sorted_numbers[$z];
    }
    
   }
 else{
     my $jj = int($arrLength / 2);
     my $halfway = ($sorted_numbers[$jj - 1] + $sorted_numbers[$jj])*0.5;
     $median = $halfway;

    for(my $y = 0; $y < $jj - 1; $y++){
	push @lowerArray, $sorted_numbers[$y];
    }

    for(my $z = $jj + 1; $z < scalar @sorted_numbers; $z++){
	push @upperArray, $sorted_numbers[$z];
    }

     
    }

  $lower = Statistics::Lite::median(@lowerArray);
  $upper = Statistics::Lite::median(@upperArray);

## Simple output.  Will add verbose output later.

  print $lower, "\t", $median, "\t", $upper, "\n";

 

