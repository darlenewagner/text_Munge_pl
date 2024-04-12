use strict;

my @firstAndMiddle = ();
my @surname = ();
my $i = 0;
my @line = ();

while(<STDIN>)  ## Read names, assume they are all on one line
  {
      if($_ =~ /^[A-Z]/i)
        {
	    @line = split(/,\s/, $_);
	}
  }

foreach(@line)
  {
      my @innerLine = split(/\s/, $line[$i]);
      my $last = pop @innerLine;
      my $the_rest = join ' ', @innerLine;
      $the_rest =~ s/^and\s//g;
      print $last, ", ", $the_rest, "\n";
      $i++;
  }
