#!/usr/bin/perl

use Bio::SeqIO;
use strict;
my $R1_File_Name = @ARGV[1];
my $R2_File_Name = @ARGV[2];

my @baseFileName = split(/(_R1|_1\.)/, $R1_File_Name);


my $R1 = Bio::SeqIO->new(
     -format  => 'fastq-illumina',
     -file    => $R1_File_Name);

my $out = Bio::SeqIO->new(
    -format  => 'fastq-illumina',
    -file    => $baseFileName[0]);



