#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $PathToDiagonalizationProgram = "/home/regnault/development/DMRG/DiagHam/build/src/Programs/QHE/QHEOnSphere/QHEBosonsTwoBodyGeneric";
my $PathToOverlapProgram = "/home/regnault/development/DMRG/DiagHam/build/src/Programs/QHE/QHEBosonsDeltaOverlap";
my $Step = 0;
my $MinValue = 0.0;
my $NbrValues = 0;
my $NbrBosons = 0;
my $LzMax = 0;
my $ReferenceVector = "";
my $PlotFlag = 0;
my $DiagonalizationProgramOptions = "";
my $GridFile = "";
my $KeepFlag = 0;


my $Result = GetOptions ("progdiag:s" => \$PathToDiagonalizationProgram, "progoverlap:s" => \$PathToOverlapProgram, 
			 "diagoption:s" => \$DiagonalizationProgramOptions,
			 "ref=s" => \$ReferenceVector, "nbrbosons:i" => \$NbrBosons, "lzmax:i" => \$LzMax,
			 "step:f" => \$Step, "min:f" => \$MinValue, "nbrvalues:i" => \$NbrValues, "plot" => \$PlotFlag,
			 "keep" => \$KeepFlag, "grid:s" => \$GridFile);



if ($GridFile eq "")
  {
    if ($NbrValues <= 0)
      {
	die ("invalid value for nbrvalues option\n");
      }
    if ($Step <= 0)
      {
	die ("invalid value for step option\n");
      }
  }

if (($ReferenceVector eq "") || (!(-e $ReferenceVector)))
  {
    die ("invalid value for ref option, or file ".$ReferenceVector." does not exist\n");
  }
if ($NbrBosons == 0)
  {
    $ReferenceVector =~ /\_n\_(\d+)/;
    $NbrBosons = $1;
    if (!defined($NbrBosons))
      {
	die ("can't find number of bosons from file name ".$ReferenceVector.", use the nbrbosons option\n");
      }
  }
if ($LzMax == 0)
  {
    $ReferenceVector =~ /\_2s\_(\d+)/;
    $LzMax = $1;
    if (!defined($LzMax))
      {
	die ("can't find lzmax from file name ".$ReferenceVector.", use the lzmax option\n");
      }
  }


my $DataFileName = "bosons_v2v0_n_".$NbrBosons."_2s_".$LzMax;
unless (open (OUTFILE, ">".$DataFileName.".overlap"))
  {
    die ("can't open ".$DataFileName.".overlap\n");
  }
close(OUTFILE);


if ($GridFile eq "")
  {
    while ($NbrValues > 0)
      {
	unless (open (OUTFILE, ">tmppseudopotential.dat"))
	  {
	    die ("can't open tmppseudopotential.dat\n");
	  }
	print OUTFILE "Pseudopotentials = 1 0 ".$MinValue;
	my $TmpLz = 3;
	while ($TmpLz <= $LzMax)
	  {
	    print OUTFILE " 0";
	    $TmpLz++;
	  }
	print OUTFILE "\n";
	close(OUTFILE);
	
	my $DiagOutputFileName = $MinValue;
	$DiagOutputFileName =~ s/\./\_/;
	$DiagOutputFileName = "v2v0".$DiagOutputFileName;
	my $Overlap = &EvalauteOverlap($PathToDiagonalizationProgram, $PathToOverlapProgram, $NbrBosons, $LzMax, $ReferenceVector, $DiagOutputFileName, $KeepFlag);
	unless (open (OUTFILE, ">>".$DataFileName.".overlap"))
	  {
	    die ("can't open ".$DataFileName.".overlap\n");
	  }
	print OUTFILE $MinValue." ".$Overlap."\n";
	close(OUTFILE);

	unlink("tmppseudopotential.dat");
	$NbrValues--;
	$MinValue += $Step;
      }
    
    if ($PlotFlag == 1)
      {
	$MinValue += $Step;
	unless (open (OUTFILE, ">".$DataFileName.".p"))
	  {
	    die ("can't create file ".$DataFileName.".p\n");
	  }
	print OUTFILE ("set xrange [0:".$MinValue."]
set yrange [0:1.1]
set xlabel \"V2/V0\"
set ylabel \"overlap\"
set size 1, 0.9
set terminal postscript landscape enhanced \"Helvetica\" 14
set output \"".$DataFileName.".ps\"
plot \"".$DataFileName.".overlap\" using 1:2 title \"N=".$NbrBosons." 2S=".$LzMax."\", \"".$DataFileName.".overlap\" using 1:2 notitle with lines
");  
	my $Command = "gnuplot ".$DataFileName.".p";
	`$Command`;
	unlink ($DataFileName.".p");
      }
  }
else
  {
    my @GridVertices;
    &ParseGridDefinition ($GridFile, \@GridVertices);
    if ($#GridVertices == 0)
      {
	die ($GridFile." is an invalid grid file\n");
      }
    my $VertexCount = 0;
    my $VertexPosition;
    foreach $VertexPosition (@GridVertices)
      {
	unless (open (OUTFILE, ">tmppseudopotential.dat"))
	  {
	    die ("can't open tmppseudopotential.dat\n");
	  }
	print OUTFILE "Pseudopotentials = ".$VertexPosition."\n";
	close(OUTFILE);

	my $DiagOutputFileName = "v2v0".$VertexCount;
	my $Overlap = 0.0;
#&EvalauteOverlap($PathToDiagonalizationProgram, $PathToOverlapProgram, $NbrBosons, $LzMax, $ReferenceVector, $DiagOutputFileName, $KeepFlag);

	unless (open (OUTFILE, ">>".$DataFileName.".overlap"))
	  {
	    die ("can't open ".$DataFileName.".overlap\n");
	  }
	print OUTFILE $VertexPosition."|".$Overlap."\n";
	close(OUTFILE);

	unlink("tmppseudopotential.dat");
	$VertexCount++;
      }
  }


# evaluate overlap betwwen maximally symmetric ground state evalauted for a given pseudopotential configuration, and a given eigenstate
#
# $_[0] = diagonalization program (with full path if any)
# $_[1] = overlap program (with full path if any)
# $_[2] = number of bosons
# $_[3] = maximum Lz value
# $_[4] = file that contains the reference vector (with path)
# $_[5] = interaction name
# $_[6] = flag to indicate if partial datas have to be kept (0 if false)
# return value = squqre overlap

sub EvalauteOverlap()
  {
    my $PathToDiagonalizationProgram = $_[0];
    my $PathToOverlapProgram = $_[1];
    my $NbrBosons = $_[2];
    my $LzMax = $_[3];
    my $ReferenceVector = $_[4];
    my $DiagOutputFileName = $_[5];
    my $KeepFlag = $_[6];


    my $ParityValue = 0;
    if ((($NbrBosons % 2) == 1) && (($LzMax % 2) == 1))
      {
	$ParityValue = 1;
      }

    my $Command = $PathToDiagonalizationProgram." -p ".$NbrBosons." -l ".$LzMax." --nbr-lz 1 -n 1 --force-reorthogonalize --eigenstate --interaction-name ".$DiagOutputFileName." --interaction-file tmppseudopotential.dat ".$DiagonalizationProgramOptions;
    `$Command`;

    $DiagOutputFileName = "bosons_".$DiagOutputFileName."_n_".$NbrBosons."_2s_".$LzMax."_lz";

    my $Overlap = -1.0;
    if (-e $DiagOutputFileName."_".$ParityValue.".0.vec")
      {
	$Command = $PathToOverlapProgram." -p ".$NbrBosons." -l ".$LzMax." --exact-state ".$ReferenceVector." --use-exact ".$DiagOutputFileName."_".$ParityValue.".0.vec";
	$Overlap = `$Command`;
	chomp ($Overlap);
	$Overlap =~ s/^overlap \= //; 
	$Overlap *= $Overlap;
      }

    if ($KeepFlag == 0)
      {
	if (-e $DiagOutputFileName."_".$ParityValue.".0.vec")
	  {
	    unlink($DiagOutputFileName."_".$ParityValue.".0.vec")
	  }
	if (-e $DiagOutputFileName.".dat")
	  {
	    unlink ($DiagOutputFileName.".dat");
	  }
      }
    return $Overlap;
  }


# parse a grid definition file
#
# $_[0] = name file that contains the grid definition
# $_[1] = reference on the array where vertex positions will be stored (each position is a string that gives pseudo-potential definition)
# $_[2] = maximum Lz value

sub ParseGridDefinition()
  {
    my $GridFile = $_[0];
    my $GridVertices = $_[1];
    my $LzMax = $_[2];

    unless (open (INFILE, $GridFile))
      {
	die ("can't open ".$GridFile."\n");
      }
    my @Pseudopotentials;
    my $TmpLine;
    while (defined($TmpLine = <INFILE>))
      {
	chomp ($TmpLine);
	$TmpLine =~ s/^\s+//;
	$TmpLine =~ s/\s+$//;
	$TmpLine =~ s/^\#.*//;
	if ($TmpLine ne "")
	  {
	    my @TmpArray = split (/\s+/, $TmpLine);
	    if (($#TmpArray != 4) || (!($TmpArray[0] =~ /^\d+$/)) || (!($TmpArray[1] =~ /^\-?\d*\.?\d+$/)) || (!($TmpArray[2] =~ /^\-?\d*\.?\d+$/))
		|| (!($TmpArray[3] =~ /^\d+$/)))
	      {
		close (INFILE);
		return;
	      }
	    my $Index = shift (@TmpArray);
	    $Pseudopotentials[$Index] = \@TmpArray;
	  }
      }
    close (INFILE);
  }
