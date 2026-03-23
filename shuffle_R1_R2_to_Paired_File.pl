#!/usr/bin/perl

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

 $fileLength = -s $R2_File_Name;

 if($fileLength == 0){
     die "Error: Input file '$R2_File_Name' is 0-length and does not contain data.\n";
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


while(my $seqR1 = $R1->next_seq){
    my $seqR2 = $R2->next_seq;
    $out->write_seq($seqR1);
    $out->write_seq($seqR2);
  }

exit;
