use strict;

## convert fasta files with single-line sequence data into multiline wrapped to 60 characters

while(<STDIN>)
{
    if($_ =~ /^[A-Z]/i)
    {
	chomp;
	for(my $i = 0; $i < length($_); $i = $i + 60)
	  {
	    print substr($_, $i, 60), "\n";
	  }
    }
    elsif($_ =~ /^>/)
    {
	print;
    }
}

