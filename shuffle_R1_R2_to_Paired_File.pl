#!/usr/bin/perl

## Bioperl-driven method for shuffling reads using filenames as positional arguments
## WARNING: This script is slow for .fastq files exceeding 500,000 reads

use Bio::SeqIO;
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
     -variant => 'sanger',
     -file    => $R1_File_Name);

my $R2 = Bio::SeqIO->new(
    -format  => 'fastq',
    -variant => 'sanger',
    -file    => $R2_File_Name);

my $out = Bio::SeqIO->new(
    -format  => 'fastq',
    -file    => ">$Paired_File_Name");

  if(!-e $Paired_File_Name){
     die "Error: Input file '$Paired_File_Name' does not exist.\n";
  }


#  $fileLength = -s $Paired_File_Name;
#  if($fileLength == 0){
#     die "Error: Input file '$Paired_File_Name' is 0-length and does not contain data.\n";
#  }

my $reads_in_R1 = `wc -l < "$R1_File_Name"`;

if($? != 0) {
    warn "WARNING: Failed to run wc command on '$R1_File_Name'\n";
}

my $reads_in_R2 = `wc -l < "$R2_File_Name"`;

if($? != 0) {
    warn "WARNING: Failed to run wc command on '$R2_File_Name'\n";
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
    $out->write_seq($seqR1);
   $out->write_seq($seqR2);
  }

exit;
