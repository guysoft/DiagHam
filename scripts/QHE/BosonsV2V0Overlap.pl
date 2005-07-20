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

my $ParityValue = 0;
if ((($NbrBosons % 2) == 1) && (($LzMax % 2) == 1))
  {
    $ParityValue = 1;
  }

if ($GridFile eq "")
  {
    while ($NbrValues > 0)
      {
	unless (open (OUTFILE, ">tmppseudopotential.dat"))
	  {
	    die ("can't open tmppseudopotential.dat\n");
	  }
	print OUTFILE "Pseudopotentials = 1 0 0.158 0 0 0 ".$MinValue;
	my $TmpLz = 7;
	while ($TmpLz <= $LzMax)
	  {
	    print OUTFILE " 0";
	    $TmpLz++;
	  }
	print OUTFILE "\n";
	close(OUTFILE);
	
	my $DiagOutputFileName = $MinValue;
	$DiagOutputFileName =~ s/\./\_/;
	my $Command = $PathToDiagonalizationProgram." -p ".$NbrBosons." -l ".$LzMax." --nbr-lz 1 -n 1 --force-reorthogonalize --eigenstate --interaction-name v2v0".$DiagOutputFileName." --interaction-file tmppseudopotential.dat ".$DiagonalizationProgramOptions;
	`$Command`;
	$DiagOutputFileName = "bosons_v2v0".$DiagOutputFileName."_n_".$NbrBosons."_2s_".$LzMax."_lz";
	if (-e $DiagOutputFileName."_".$ParityValue.".0.vec")
	  {
	    $Command = $PathToOverlapProgram." -p ".$NbrBosons." -l ".$LzMax." --exact-state ".$ReferenceVector." --use-exact ".$DiagOutputFileName."_".$ParityValue.".0.vec";
	    my $Overlap = `$Command`;
	    chomp ($Overlap);
	    $Overlap =~ s/^overlap \= //; 
	    unless (open (OUTFILE, ">>".$DataFileName.".overlap"))
	      {
		die ("can't open ".$DataFileName.".overlap\n");
	      }
	    print OUTFILE $MinValue." ".($Overlap * $Overlap)."\n";
	    close(OUTFILE);
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
      }
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

