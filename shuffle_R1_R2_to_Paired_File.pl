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




my $R1 = Bio::SeqIO->new(
     -format  => 'fastq',
     -variant => 'illumina',
     -file    => $R1_File_Name);

my $R2 = Bio::SeqIO->new(
    -format  => 'fastq',
    -variant => 'illumina',
    -file    => $R2_File_Name);

my $out = Bio::SeqIO->new(
    -format  => 'fastq',
    -variant => 'illumina',
    -file    => ">$Paired_File_Name");

while(my $seqR1 = $R1->next_seq){
    $out->write_seq($seqR1);
  }

exit;
