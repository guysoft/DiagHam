#!/usr/bin/perl -w

use strict 'vars';

use Getopt::Long;


my $SpectrumFile = "";
my $Eigenvalue = 0.0;
my $Error = 1.0e-12;
my $LFlag = 0;
my $RowFlag = 0;
my $Result = GetOptions ("spectrum=s" => \$SpectrumFile, "eigenvalue:s" => \$Eigenvalue, "error:s" => \$Error, "lsort" => \$LFlag,
			"row" => \$RowFlag); 

#if (($SpectrumFile eq "") || (!(-e $SpectrumFile)) || (!($Eigenvalue =~ /^[\+\-]?\d*\.?\d*e?\d+$/)) || (!()))
if ($SpectrumFile eq "")
  {
    die ("usage: SphereSpectrumDegenracy.pl --spectrum file_name [--eigenvalue 0.0 --error 1.0e-12 --lsort]\n");
  }

if (abs($Eigenvalue) < $Error)
  {
    $Eigenvalue = 0.0;
  }
my %Degeneracy;

unless (open (INFILE, $SpectrumFile))
  {
    die ("can't open ".$SpectrumFile."\n");
  }
my $TmpLine;
while (defined($TmpLine = <INFILE>))
  {
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    if ((abs ($TmpArray[1] - $Eigenvalue) < abs($Error * $Eigenvalue)) || (($Eigenvalue == 0.0) && (abs ($TmpArray[1]) < $Error)))
      {
	if (defined($Degeneracy{$TmpArray[0]}))
	  {
	    $Degeneracy{$TmpArray[0]}++;
	  }
	else
	  {
	    $Degeneracy{$TmpArray[0]} = 1;
	  }
      }
  }
close (INFILE);

if ($LFlag != 0)
  {
    my $CurrentNbr = 0;
    foreach $TmpLine (sort {$b <=> $a} (keys(%Degeneracy)))
      {	
	$Degeneracy{$TmpLine} -= $CurrentNbr;
	$CurrentNbr += $Degeneracy{$TmpLine};
      }
    
  }

if ($RowFlag == 0)
  {
    foreach $TmpLine (sort {$a <=> $b} (keys(%Degeneracy)))
      {
	print $TmpLine." ".$Degeneracy{$TmpLine}."\n";
      }
  }
else
  {
    foreach $TmpLine (sort {$a <=> $b} (keys(%Degeneracy)))
      {
	print $Degeneracy{$TmpLine}." ";
      }
    print "\n";
  }
