#!/usr/bin/perl -w

use strict 'vars';


if (!(defined($ARGV[0])))
  {
    die ("usage: MagneticTranslation2Translation filename_mt [filename]\n");
  }

unless (open (INFILE, $ARGV[0]))
  {
    die ("can't open ".$ARGV[0]."\n");
  }
my $TmpLine;
my %Spectrum;
foreach $TmpLine (<INFILE>)
  {
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    my $TmpKey = $TmpArray[1];
    if (!(defined($Spectrum{$TmpKey})))
      {
	my @Array = ($TmpArray[2]);
	$Spectrum{$TmpKey} = \@Array;
      }
    else
      {
	my $Array = $Spectrum{$TmpKey};
	push (@$Array, $TmpArray[2]);
      }
  }
close (INFILE);

if (defined($ARGV[1]))
  {
    my %Spectrum2;

   unless (open (INFILE, $ARGV[1]))
      {
	die ("can't open ".$ARGV[1]."\n");
      }
    foreach $TmpLine (<INFILE>)
      {
	chomp ($TmpLine);
	my @TmpArray = split (/ /, $TmpLine);
	my $TmpKey = $TmpArray[0];
	if (!(defined($Spectrum2{$TmpKey})))
	  {
	    my @Array = ($TmpArray[1]);
	    $Spectrum2{$TmpKey} = \@Array;
	  }
	else
	  {
	    my $Array = $Spectrum2{$TmpKey};
	    push (@$Array, $TmpArray[1]);
	  }
      }
    close (INFILE);

    foreach $TmpLine (sort {$a <=> $b} (keys(%Spectrum)))
      {
	my $Array = $Spectrum{$TmpLine};
	my $Array3 = $Spectrum2{$TmpLine};
	my @Array2 = sort {$a <=> $b} (@$Array3);
	my $Value;
	my $Index = 0;
	foreach $Value (sort {$a <=> $b} (@$Array))
	  {
	    print $TmpLine." ".$Value." ".$Array2[$Index]." ";
	    my $Error = abs($Array2[$Index] - $Value);
	    if ($Error > (1e-13 * abs($Value)))
	      {
		print $Error;
		print STDERR $TmpLine." ".$Value." ".$Array2[$Index]." ".$Error."\n";
	      }
	    else
	      {
		print "0";
	      }
	    print "\n";
	    $Index++;
	  }
      }

  }
else
  {
    foreach $TmpLine (sort {$a <=> $b} (keys(%Spectrum)))
      {
	my $Array = $Spectrum{$TmpLine};
	my $Value;
	foreach $Value (sort {$a <=> $b} (@$Array))
	  {
	    print $TmpLine." ".$Value."\n";
	  }
      }
  }

