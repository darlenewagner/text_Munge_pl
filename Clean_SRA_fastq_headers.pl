#!/usr/bin/perl
use strict;
use Getopt::Long;
use warnings;
use Data::Dumper;
## Input BlastN file is tab-delimited with 'virus' or '10239' in 17th column
## Output format option from original BlastN command had 18 columns:
## -outfmt '6 qseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore sgi sacc staxids sblastnames stitle'

## Alphanumeric Genbank or RefSeq expected
my $accession = '';

## Two or more words expected, the first one strictly alphabetical
## Is empty string until assigned after $accession fails format test
my $virusName = '';

## For searching when $accession or $virusName are not supplied
my $rotavirusStr = qr/(Human\srotavirus\s(A|B|C|D|E|F|G)|Rotavirus\s(A|B|C|D|E|F|G))/;
my @HighConsequenceArray = ();
my %Sortable;

## A Hash of Hashes for organizing rotavirus segment matches
my %RotavirusSegments = ();

$accession = $ARGV[0];


## Help section and usage
my $help = '';
my $verbose = '';
my $sortIt = '';
my $rotavirus = '';

GetOptions(
    'help|h|'  => \$help,
    'verbose'  => \$verbose,
    'sort'     => \$sortIt,
    'rotavirus|r|' => \$rotavirus,
    ) or die "Error in command line option arguments\n";

if($help)
   {
       print "\nUsage: BlastN_Contig_Match.pl --rotavirus --sort < BlastN-Output.tsv\n";
       print "\n<------------------------------------------------- Options -------------------------------------------------->\n";
       print "ACCESSION or VIRUS NAME STR#, Search key (required)\n";
       print "--verbose,                    Output column headers (optional)\n";
       print "--sort,                       Sort output by MatchID (NCBI Accession) and AlnLength\n";
       print "--rotavirus,                  Search Human rotavirus or rotavirus (no search key needed)\n";
       print "--help,                       Display this help message\n\n";
       exit;
   } 


## Determine input string '$accession' type

if($accession =~ /^[A-Z][a-z]+\s+([A-Za-z][A-Za-z0-9]?[0-9]?\s)?([A-Za-z0-9]+(-|_)?[A-Za-z0-9]+|sp|str|strain)\.?\s+([A-Z]\s)?[A-Za-z0-9]+(-|_)?[A-Za-z0-9]+/)
  {
    $virusName = $ARGV[0];
     if($verbose)
        {
	    print "Search by species/strain name:\n";
	    print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
  }
elsif($accession =~ /^Human\srotavirus\s(A|B|C|D|E|F|G)/)
  {
      $virusName = $ARGV[0];
      if($verbose)
      {
	   print "Search Human rotavirus without categorization of segments:\n";
	   print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
      }
  }
elsif($accession =~ /^[A-Z][A-Z][A-Z0-9]+\.?[0-9]$/)
  {
    if($verbose)
      {
	  print "Search by GenBank accession:\n";
	  print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
      }
  }
elsif($accession =~ /^N(C|T|X|Z)_[0-9]+\.?[0-9]$/)
  {
     if($verbose)
       {
	  print "Search by RefSeq accession:\n";
	  print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
       }
  }
elsif(!$rotavirus)
  {
      die "Search string too short and/or lacks sufficient information --$!";
  }


## Iterate through blastN output
## Must be tabular format with:
## -outfmt '6 qseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore sgi sacc staxids sblastnames stitle'


my $i = 0;
my @line = ();

while(<STDIN>)  ## Read BlastN file, use mode according to which string is nonempty
  {
      if(($accession =~ /^[A-Z][A-Z][A-Z0-9]+\.?[0-9]$/) || ($accession =~ /^N(C|T|X|Z)_[0-9]+\.?[0-9]$/))
        { 
	    @line = split(/\t/, $_);
	    if($line[16] =~ /virus|10239/)
	    {
		chomp $line[17];
		my $newLine = '';
		if($verbose)
		  {
		    $newLine = $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[14]."\t".$line[17]."\n";
		  }
		else
		 {
                     $newLine = $line[0]."\t\t".$line[1]."\t\t".$line[2]."\t\t".$line[14]."\t\t".$line[17]."\n";
		 }
		if($newLine =~ m/$accession/)
		  {
               	    if($sortIt)
		      {
			  $Sortable{$line[14]}{$line[2]} = $newLine;
		      }
		    else
		      {
			  print $newLine;
		      }
		  }
		 elsif(($newLine =~ m/$rotavirusStr/) && ($rotavirus))
	          {
		    push @HighConsequenceArray, $newLine;
		 }
	     }
	}   
      elsif($virusName =~ /^[A-Z][a-z]+\s+([A-Za-z][A-Za-z0-9]?[0-9]?\s)?([A-Za-z0-9]+(-|_)?[A-Za-z0-9]+|sp|str|strain)\.?\s+([A-Z]\s)?[A-Za-z0-9]+(-|_)?[A-Za-z0-9]+/)  
        {
	    @line = split(/\t/, $_);
	    if($line[16] =~ /virus|10239/)
	    {
		chomp $line[17];
		my $newLine = '';
		if($verbose)
		 {
		    $newLine = $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[14]."\t".$line[17]."\n";
		 }
		else
		 {
                     $newLine = "'".$line[0]."'\t'".$line[1]."'\t'".$line[2]."'\t'".$line[14]."'\t'".$line[17]."'\n";
		 }
		if($newLine =~ m/$virusName/)
		{
		    #print $newLine;
		   if($sortIt)
		      {
			  $Sortable{$line[14]}{$line[2]} = $newLine;
		      }
		    else
		      {
                         print $newLine
		      }

		}
		elsif(($newLine =~ m/$rotavirusStr/) && ($rotavirus))
	        {
                   push @HighConsequenceArray, $newLine;
		}
	    }
     	}
      elsif($rotavirus)
        {
	    @line = split(/\t/, $_);
      	    if($line[16] =~ /virus|10239/)
	    {
		chomp $line[17];
		my $newLine = '';

		if($verbose)
		 {
                    $newLine = $line[0]."\t".$line[1]."\t".$line[2]."\t".$line[14]."\t".$line[17]."\n";
		 }
		else
		 {
 		    $newLine = "'".$line[0]."'\t'".$line[1]."'\t'".$line[2]."'\t'".$line[14]."'\t'".$line[17]."'\n";
                 }
		
      		if($newLine =~ m/$rotavirusStr/)
      		{
                    #push @HighConsequenceArray, $newLine;
                    
		    if($sortIt && ( ($newLine =~ /segment\s1(\s|,)/i) || ($newLine =~ /RNA\-dependent\sRNA\spolymerase\sVP1\s\(VP1\)\sgene/) || ($newLine =~ /VP1\sgene\sfor\sVP1\sprotein/)))
		      {
			  $RotavirusSegments{'segment_1'}{$line[14]}{$line[2]} = $newLine;
			  
		      }
		    elsif($sortIt && ( ($newLine =~ /segment\s2(\s|,)/i) || ($newLine =~ /(core\scapsid\sprotein|RNA\sviral\sgenome\sbinding)\sVP2\s\(VP2\)\sgene/)))
		    {
                        $RotavirusSegments{'segment_2'}{$line[14]}{$line[2]} = $newLine;
                    }
		    elsif($sortIt && ( ($newLine =~ /segment\s3(\s|,)/i) || ($newLine =~ /RNA\scapping\s(protein|enzyme)\sVP3\s\(VP3\)\sgene/)))
		    {
                        $RotavirusSegments{'segment_3'}{$line[14]}{$line[2]} = $newLine;
                    }
		    elsif($sortIt && (($newLine =~ /segment\s4(\s|,)/i) || ($newLine =~ /outer\scapsid\sspike\sprotein\sVP4\s\(VP4\)\sgene/)))
		    {
                        $RotavirusSegments{'segment_4'}{$line[14]}{$line[2]} = $newLine;
                    }
		    elsif($sortIt && (($newLine =~ /segment\s5(\s|,)/i) || ($newLine =~ /non\-structural\sprotein\s1\s\(NSP1\)\sgene/)))
		    {
                        $RotavirusSegments{'segment_5'}{$line[14]}{$line[2]} = $newLine;
                    }
		    elsif($sortIt && ( ($newLine =~ /segment\s6(\s|,)/i) || ($newLine =~ /inner\scapsid\sprotein\sVP6\s\(VP6\)\sgene/)))
		    {
                        $RotavirusSegments{'segment_6'}{$line[14]}{$line[2]} = $newLine;
                    }
		    elsif($sortIt && ( ($newLine =~ /segment\s7(\s|,)/i) || ($newLine =~ /non\-structural\sprotein\s3\s\(NSP3\)\sgene/)))
		    {
                        $RotavirusSegments{'segment_7'}{$line[14]}{$line[2]} = $newLine;
                    }
		    elsif($sortIt && ( ($newLine =~ /segment\s8(\s|,)/i) || ($newLine =~ /non\-structural\sprotein\s2\s\(NSP2\)\sgene/)))
		    {
                        $RotavirusSegments{'segment_8'}{$line[14]}{$line[2]} = $newLine;
                    }
		    elsif($sortIt && ( ($newLine =~ /segment\s9(\s|,)/i) || ($newLine =~ /capsid\sglycoprotein\sVP7\s\(VP7\)\sgene/) || ($newLine =~ /outer\scapsid\sglycoprotein\sVP7\s\(VP7\)\sgene/) || ($newLine =~ /NSP3\sgene,/)))
		    {
                        $RotavirusSegments{'segment_9'}{$line[14]}{$line[2]} = $newLine;
                    }
		    elsif($sortIt && ( ($newLine =~ /segment\s10(\s|,)/i) || ($newLine =~ /non\-structural\sprotein\s4\s\(NSP4\)\sgene/)))
		      {
			  $RotavirusSegments{'segment_10'}{$line[14]}{$line[2]} = $newLine;
		      }
		    elsif($sortIt && ( ($newLine =~ /segment\s11(\s|,)/i) || ($newLine =~ /non\-structural\sprotein\s5\s\(NSP5\)\sand\snon-structural\sprotein\s6\s\(NSP6\)\sgenes/) ))
		      {
			  $RotavirusSegments{'segment_11'}{$line[14]}{$line[2]} = $newLine;
    	              }

      		}

      	    }
	}
      #else
      #{
      #	  die "Search string too short, lacks sufficient information, or --highConsequence option not provided...$!";
      #}
  }

## Option '--sort' with search string included

if($sortIt && !$rotavirus)
{
    foreach my $acc (sort keys %Sortable)
    {
	foreach my $aln (sort { $b <=> $a } keys %{ $Sortable{$acc}} )
	{
	    print $Sortable{$acc}{$aln};
	}
    }
    
}


## When option '--highConsequence' is used, no search string is needed

if($rotavirus)
{
    if($sortIt)
    {
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 1 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      my $i = 0;
      if( exists $RotavirusSegments{'segment_1'} )
         {
	    foreach my $acc (sort keys %{ $RotavirusSegments{'segment_1'} })
	    {
		foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_1'}{$acc}} )
		{
		    print $RotavirusSegments{'segment_1'}{$acc}{$aln};
		}
	    }
	 }
      else
      {
    	    print "Segment 1 not found.\n";	  
      }
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 2 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      if( exists $RotavirusSegments{'segment_2'} )
       {		
	    foreach my $acc (sort keys %{ $RotavirusSegments{'segment_2'} })
	    {
		foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_2'}{$acc}} )
		{
		    print $RotavirusSegments{'segment_2'}{$acc}{$aln};
		}
	    }
       }
      else
	{
    	    print "Segment 2 not found.\n";          
	}
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 3 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      if( exists $RotavirusSegments{'segment_3'} )
	    {
          foreach my $acc (sort keys %{ $RotavirusSegments{'segment_3'} })
	    {
		foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_3'}{$acc}} )
		{
		    print $RotavirusSegments{'segment_3'}{$acc}{$aln};
		}
	    }
	    }
      else
	    {
                 print "Segment 3 not found.\n"; 
	    }
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 4 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      if( exists $RotavirusSegments{'segment_4'} )
	    {
	    foreach my $acc (sort keys %{ $RotavirusSegments{'segment_4'} })
	    {
		foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_4'}{$acc}} )
		{
		    print $RotavirusSegments{'segment_4'}{$acc}{$aln};
		}
	    }
	    }
      else
	    {
       	       print "Segment 4 not found.\n";		
	    }
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 5 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      if( exists $RotavirusSegments{'segment_5'} )
         {
	    foreach my $acc (sort keys %{ $RotavirusSegments{'segment_5'} })
	    {
		foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_5'}{$acc}} )
		{
		    print $RotavirusSegments{'segment_5'}{$acc}{$aln};
		}
	    }
         }
      else
        {
    	    print "Segment 5 not found.\n";
        }
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 6 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      if( exists $RotavirusSegments{'segment_6'} )
        {      
	    foreach my $acc (sort keys %{ $RotavirusSegments{'segment_6'} })
	    {
		foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_6'}{$acc}} )
		{
		    print $RotavirusSegments{'segment_6'}{$acc}{$aln};
		}
	    }
        }
      else
	    {
              	print "Segment 6 not found.\n";
	    }
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 7 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      if( exists $RotavirusSegments{'segment_7'} )
	{	
	    foreach my $acc (sort keys %{ $RotavirusSegments{'segment_7'} })
	    {
		foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_7'}{$acc}} )
		{
		    print $RotavirusSegments{'segment_7'}{$acc}{$aln};
		}
	    }
	}
      else
	  {
              	print "Segment 7 not found.\n";
	  }
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 8 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      if( exists $RotavirusSegments{'segment_8'} )
	    {
	    foreach my $acc (sort keys %{ $RotavirusSegments{'segment_8'} })
	    {
		foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_8'}{$acc}} )
		{
		    print $RotavirusSegments{'segment_8'}{$acc}{$aln};
		}
	    }
	    }
      else
	    {
		print "Segment 8 not found.\n";
	    }
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 9 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      if( exists $RotavirusSegments{'segment_9'} )
	    {
	    foreach my $acc (sort keys %{ $RotavirusSegments{'segment_9'} })
	    {
		foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_9'}{$acc}} )
		{
		    print $RotavirusSegments{'segment_9'}{$acc}{$aln};
		}
	    }
	    }
      else
	    {
		print "Segment 9 not found.\n";
	    }
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 10 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}
      if( exists $RotavirusSegments{'segment_10'} )
	    {
   	        foreach my $acc (sort keys %{ $RotavirusSegments{'segment_10'} })
	        {
		    foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_10'}{$acc}} )
		    {
		        print $RotavirusSegments{'segment_10'}{$acc}{$aln};
		    }
	        }
	    }
      else
	    {
		 print "Segment 10 not found.\n";		
	    }
      
      if($verbose)
        {
           print "\n<------- Search for Rotavirus Segment 11 ------->\n";
           print "QueryID\t%Ident.\tAlnLength\tMatchID\tMatchName\n";
	}

            if( exists $RotavirusSegments{'segment_11'} )
              {
	       foreach my $acc (sort keys %{ $RotavirusSegments{'segment_11'} })
	         { 
		   foreach my $aln (sort { $b <=> $a } keys %{ $RotavirusSegments{'segment_11'}{$acc}} )
		    {
		       print $RotavirusSegments{'segment_11'}{$acc}{$aln};
		    }
		
	         }
	      }
            else
             {
		 print "Segment 11 not found.\n";
	     }
	


    }
      else
      {
         die "Sort option required for rotavirus search.\n";
      }
  }

