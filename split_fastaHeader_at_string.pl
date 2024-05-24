use strict;

## split fasta header at input string, $trimmer, to shorten header in fasta input file, $fasta

my $trimmer = $ARGV[0];

chomp $trimmer;

if($trimmer eq '-s')
  {
      $trimmer = '\s';
  }

my $fasta = $ARGV[1];

open(FASTA, $fasta) || die "Can't find fasta-formatted input file, $fasta $!";

while(<FASTA>)
{
    if($_ =~ /^>/)
    {
	my @newHeader = split $trimmer, $_;
	print $newHeader[0], "\n";
    }
    else
    {
       print $_;
    }
}

close(FASTA)

