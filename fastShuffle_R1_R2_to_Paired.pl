#!/usr/bin/perl

## Bioperl-driven method for shuffling reads using filenames as positional arguments
## WARNING: This script is slow for .fastq files exceeding 500,000 reads

use Bio::SeqIO::fastq;
use strict;
use warnings;

# Shuffle reads from forward and reverse reads in .fastq
# Validate that there are 2 input files and 1 output file
if(@ARGV != 3) {
    die "Usage: $0 Reads_R1_001.fastq Reads_R2_001.fastq Reads_Pairs.fastq\n";
}

my ( $R1_File_Name, $R2_File_Name, $Paired_File_Name ) = @ARGV;

# Validate forward reads .fastq file existence and/or nonzero length
  if(!-e $R1_File_Name){
     die "Error: Input file '$R1_File_Name' does not exist.\n";
  }

  if(!-e $R2_File_Name){
     die "Error: Input file '$R2_File_Name' does not exist.\n";
  }

  if(($R1_File_Name =~ /\.gz$/) || ($R2_File_Name =~ /\.gz$/))
  {
      my $new_R1 = `gunzip $R1_File_Name`;
      my $new_R2 = `gunzip $R2_File_Name`;
      $R1_File_Name = $new_R1;
      $R2_File_Name = $new_R2;
  } 


my $fileLength = -s $R1_File_Name;

  if($fileLength == 0){
     die "Error: Input file '$R1_File_Name' is 0-length and does not contain data.\n";
  }

my $fileLength2 = -s $R2_File_Name;

 if($fileLength2 == 0){
     die "Error: Input file '$R2_File_Name' is 0-length and does not contain data.\n";
  }

 if(($fileLength + $fileLength2) > 1500000000){
      print "Total size of files exceeds 1.5 GB.\n";
      print "Are you sure you wish to continue? ";
      my $response = <STDIN>;
      chomp $response;
      if($response ne 'Y'){
	  die "Aborting read shuffle operation.\n";
      }

 }
    


my $R1 = Bio::SeqIO->new(
     -format  => 'fastq',
     -file    => $R1_File_Name,
     -verbose => -1,
     -noclose => 1,
     -streambased => 1);

my $R2 = Bio::SeqIO->new(
    -format  => 'fastq',
    -file    => $R2_File_Name,
    -verbose => -1,
    -noclose => 1,
    -streambased => 1);

open(OUT, '>', $Paired_File_Name) || die "Output file name $Paired_File_Name not provided $!";

my $reads_in_R1 = `wc -l < "$R1_File_Name"`;

if($? != 0) {
    warn "WARNING: Failed to run wc command on '$R1_File_Name'\n";
}

my $header_R1 = `head -1 "$R1_File_Name"`;

my $strand_R1 = "1";

if(scalar split(" ", $header_R1) > 1)
 {
     my @id_line = split(" ", $header_R1);
     $strand_R1 = $id_line[1];
 }


my $reads_in_R2 = `wc -l < "$R2_File_Name"`;

if($? != 0) {
    warn "WARNING: Failed to run wc command on '$R2_File_Name'\n";
}

my $header_R2 = `head -1 "$R2_File_Name"`;

my $strand_R2 = "2";

if(scalar split(" ", $header_R2) > 1)
 {
     my @id_line = split(" ", $header_R2);
     $strand_R2 = $id_line[1];
 }


$reads_in_R1 = $reads_in_R1 / 4;
$reads_in_R2 = $reads_in_R2 / 4; 

printf "R1 read count: %i\n", $reads_in_R1;
printf "R2 read count: %i\n", $reads_in_R2;




if($reads_in_R1 != $reads_in_R2)
  {
      print "There are unpaired reads indicated by unequal read counts.\n";
      print "Do you wish to continue anyway? (Y/n)";
      my $response = <STDIN>;
      chomp $response;
       
      if($response ne 'Y'){
	  die "Aborting read shuffle operation.\n";
      }
  }

while(my $seqR1 = $R1->next_seq){
    my $seqR2 = $R2->next_seq;

    my $id = $seqR1->id;
    my $strand = $strand_R1;
    my $dna = $seqR1->seq;
    my $qual = join("", map{chr($_ + 33)} @{$seqR1->qual});

    print OUT $id, " ", $strand, "\n", $dna, "\n+\n", $qual, "\n";

    $id = $seqR2->id;
    $strand = $strand_R2;
    $dna = $seqR2->seq;
    $qual = join("", map{chr($_ + 33)} @{$seqR2->qual});

    print OUT $id, " ", $strand, "\n", $dna, "\n+\n", $qual, "\n";

}


close(OUT);

exit;
