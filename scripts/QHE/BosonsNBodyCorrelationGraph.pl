#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $FractionDescriptionFile = "";
my $NbrNBody = 3;

my $Result = GetOptions ("fraction=s" => \$FractionDescriptionFile); 


my @Fractions;
my @FractionNames;
unless (open (INFILE, $FractionDescriptionFile))
  {
    die ("can't open ".$FractionDescriptionFile."\n");
  }
my $TmpLine;
foreach $TmpLine (<INFILE>)
  {
    chomp ($TmpLine);
    my @TmpArray = split (/\|/, $TmpLine);
    push (@Fractions, $TmpArray[0] / $TmpArray[1]);
    if ($TmpArray[1] != 1)
      {
	push (@FractionNames, $TmpArray[0]."/".$TmpArray[1]." ".$TmpArray[4]);
      }
    else
      {
	push (@FractionNames, $TmpArray[0]." ".$TmpArray[4]);
      }
  }
close (INFILE);

my $TmpFileName;
my $Index = 0;
while ($Index <= $#FractionNames)
  {
    my $CurrentFraction = $Fractions[$Index];
    my %Correlations;
    my $Flag = 0;
    foreach $TmpFileName (<correlations*body.dat>)
      {
	$TmpFileName =~ /correlations(\d*)body\.dat/;
	my $NbrNbody =  $1;
	my $Values = "";
	unless (open (INFILE, $TmpFileName))
	  {
	    die ("can't open ".$TmpFileName."\n");
	  }
	while (defined($TmpLine = <INFILE>))
	  {
	    chomp ($TmpLine);
	    my @TmpArray = split (/ /, $TmpLine);
	    if (int($TmpArray[1] * 1000) == int($CurrentFraction * 1000))
	      {
		$Values .= (1.0 / $TmpArray[0])." ".($TmpArray[2] / ($TmpArray[0] ** $NbrNbody))."\n";
	      }
	  }
	close (INFILE);
	if ($Values ne "")
	  {
	    $Flag = 1;
	    $Correlations{$NbrNbody } = $Values;
	  }
      }
    if ($Flag != 0)
      { 
	my $FractionName = $FractionNames[$Index];
	$FractionName =~ s/\//\_/g;
	$FractionName =~ s/ /\_/g;
	$FractionName = "correlation_".$FractionName.".ps";
	&PlotCorrelations (\%Correlations, $FractionNames[$Index], $FractionName);
      }
    $Index++;
  }


# plot correlations for a given fraction
#
# $_[0] = reference on the hash table containing the size and corresponding correlation for each nbody value
# $_[1] = fraction name

sub PlotCorrelations
  {
    my $Correlations = $_[0];
    my $FractionName = $_[1];
    my $OutputFile = $_[2];
    my $NbrNbody;
    my $LastNbrNbody = 0;
    my $Values;
    my $TmpFilename = CreateTemporaryFileName()."nbody";
    while (($NbrNbody, $Values) = each (%$Correlations))
      {
	$LastNbrNbody = $NbrNbody;
	unless (open (OUTFILE, ">".$TmpFilename.$NbrNbody))
	  {
	    die ("can't create file ".$TmpFilename.$NbrNbody."\n");
	  }
	print OUTFILE $Values;
	close (OUTFILE);
      }
    unless (open (OUTFILE, ">".$TmpFilename.".p"))
      {
	die ("can't create file ".$TmpFilename.".p\n");
      }
    print OUTFILE ("set xrange [0:0.3]
set xlabel \"1/N\"
set ylabel \"corr\"
set size 1.5, 0.6
set terminal postscript portrait enhanced \"Helvetica\" 14
set output \"".$OutputFile."\"
");
    my @SortedNbrNbody = keys(%$Correlations);
    @SortedNbrNbody = sort {$a <=> $b} (@SortedNbrNbody);
    foreach $NbrNbody (@SortedNbrNbody)
      {
	print OUTFILE ("f".$NbrNbody."(x)= a".$NbrNbody."*x+b".$NbrNbody."
fit f".$NbrNbody."(x) \"".$TmpFilename.$NbrNbody."\" using 1:2 via a".$NbrNbody.",b".$NbrNbody."\n");
      }
    print OUTFILE ("plot ");
    my $Index = 5;
    my $Index2 = 1;
    foreach $NbrNbody (@SortedNbrNbody)
      {
	print OUTFILE "\"".$TmpFilename.$NbrNbody."\" using 1:2 title \"".$FractionName." m=".$NbrNbody."\" with points ".$Index.", f".$NbrNbody."(x) notitle with lines 1";
	$Index++;
	$Index2++;
	if ($NbrNbody != $SortedNbrNbody[$#SortedNbrNbody])
	  {
	    print OUTFILE ", ";
	  }

      }
    close (OUTFILE);
    my $Command = "gnuplot ".$TmpFilename.".p";
    `$Command`;
    while (($NbrNbody, $Values) = each (%$Correlations))
      {
	unlink ($TmpFilename.$NbrNbody);
      }    
#    unlink ($TmpFilename.".p");
  }

# create a temporary file name
#
# return value =temporary file name
 
sub CreateTemporaryFileName
  {
    return "tmp".time();
  }

