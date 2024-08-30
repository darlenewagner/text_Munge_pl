#!/usr/bin/perl
use strict;
# Prepare metagenomic profile data for presentation in .html via Krona/2.7
# Create krona-readable metagenomic profiles in seven-column:
#name    taxID   taxRank genomeSize      numReads        numUniqueReads  abundance
# or four-column format:
#name    taxID    taxRank    abundance

$accession = $ARGV[0];

my $help = '';
my $verbose = '';
my $advanced = '';

GetOptions(
    'help|h|'  => \$help,
    'verbose|v|'  => \$verbose,
    'sevenColumn|s|' => \$advanced,
    ) or die "Error in command line option arguments\n";

if($help)
   {
       print "\nUsage: CAMI_format_to_Krona.pl --verbose --sevenColumn  MetAlgin_or_MetaPhlAn.format.file\n";
       print "\n<------------------------------------------- Options -------------------------------------------->\n";
       print "--verbose,                Show original column names in STDERR\n";
       print "--sevenColumn,            Output in seven-column format (four-column default)\n";
       print "--help,                   Display this help message\n\n";
       exit;
   }


